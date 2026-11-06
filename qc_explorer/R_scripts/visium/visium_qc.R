#!/usr/bin/env Rscript
###############################################################################
# visium_qc.R — sample-sheet-driven quality control for 10x Visium (human data)
#
# Usage
#   Rscript visium_qc.R <samplesheet.csv> [out_dir]
#   ...or open bladder cancer.Rproj in RStudio, set SAMPLE_SHEET_PATH below, and click Source.
#   (After sourcing once, you can also call main("path/to/sheet.csv") yourself.)
#
# The sample sheet has ONE ROW PER SAMPLE (one capture area / slide section).
# Required columns: sample_id, path. Every other column is optional; a blank
# cell means "use the default" in DEFAULTS below. Each sample is QC'd and saved
# as its OWN Seurat object (with its tissue image); nothing is merged.
#
# Per sample, the script:
#   1. reads the counts + spatial files (Space Ranger output, including GEO-style
#      prefixed/flattened files and .gz images), or a Seurat .rds
#   2. keeps in-tissue spots only
#   3. computes QC metrics per spot (UMIs, genes, % mito; % ribo and % hemoglobin reported only)
#   4. removes outlier spots with per-sample MAD thresholds (same method as scrna_qc.R)
#   5. drops rare genes and saves <sample_id>_qc.rds (raw counts, not normalized)
# Per dataset folder it also writes summary tables, QC plots, and spatial maps
# showing WHERE the removed spots are on the tissue.
#
# Not done here, on purpose: doublet detection and ambient-RNA correction
# (a Visium spot is a mixture of cells by design) and normalization (do it
# downstream; RCTD deconvolution needs raw counts).
###############################################################################

## =============================================================================
## Settings (used when no command-line arguments are given)
## =============================================================================
# Relative to the repo root (where "bladder cancer.Rproj" is).
SAMPLE_SHEET_PATH <- "datasets/VISIUM_DATASET/samplesheet_visium.csv"   # <- change to your dataset
OUT_DIR           <- NULL   # NULL = "<folder of the sample sheet>/qc_output"

## =============================================================================
## Defaults for every sample-sheet column
## =============================================================================
DEFAULTS <- list(
  # what the sample is
  dataset        = "dataset1", # output folder name; samples are never merged
  format         = "auto",     # auto | spaceranger | rds
  gene_column    = 2,          # features file: 1 = Ensembl ID, 2 = gene symbol
  mito_pattern   = "^MT-",     # regex for mitochondrial genes (human gene symbols)
  # spot filters (per-sample MAD cutoffs, with hard floors/caps and a mito clamp)
  nmads_low      = 3,          # MADs below the median (log scale) for UMIs/genes
  nmads_high     = 5,          # MADs above the median (log scale) for UMIs/genes
  nmads_mito     = 3,          # MADs above the median for % mito
  min_counts     = 250,        # hard floor on UMIs per spot (always applied)
  min_features   = 200,        # hard floor on genes per spot (always applied)
  max_counts     = NA,         # optional hard cap on UMIs
  max_features   = NA,         # optional hard cap on genes
  mito_floor     = 10,         # MAD mito cutoff is never stricter than this
  mito_ceiling   = 25,         # MAD mito cutoff is never looser than this
  # genes and reproducibility
  min_spots_gene = 3,          # keep genes seen in >= this many of the sample's kept spots
  seed           = 1234
)
LOGICAL_PARAMS   <- character(0)
CHARACTER_PARAMS <- c("dataset", "format", "mito_pattern")
NUMERIC_PARAMS   <- setdiff(names(DEFAULTS), c(LOGICAL_PARAMS, CHARACTER_PARAMS))
KNOWN_FORMATS    <- c("spaceranger", "rds")

# Human gene-symbol patterns for ribosomal and hemoglobin genes (reported only; never used to filter)
HUMAN_GENES <- list(
  ribo = "^RP[SL][0-9]+[AXYL]?[0-9]*$|^RPLP[0-2]$|^RPSA$",
  hb   = "^HB[ABDEGMQZ][0-9]?$"
)

REASON_COLOURS <- c(pass = "grey75", low_counts = "#E69F00", low_features = "#D55E00",
                    high_counts = "#56B4E9", high_features = "#0072B2", mito = "#CC79A7")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

## =============================================================================
## Logging: every message goes to the console and to <out_dir>/qc_log.txt
## =============================================================================
.qc_log <- new.env()

start_log <- function(out_dir) {
  .qc_log$file <- file.path(out_dir, "qc_log.txt")
  cat("", file = .qc_log$file)
}

log_msg <- function(...) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(...))
  message(line)
  if (!is.null(.qc_log$file)) cat(line, "\n", file = .qc_log$file, append = TRUE, sep = "")
}

## =============================================================================
## Sample sheet: read, fill in defaults, check everything before running
## =============================================================================
is_absolute <- function(p) grepl("^([A-Za-z]:|/|~)", p)

# Relative paths in the sheet are relative to the sheet's own folder
resolve_path <- function(p, sheet_dir) {
  if (is.na(p)) return(NA_character_)
  p <- gsub("\\\\", "/", p)
  if (!is_absolute(p)) p <- file.path(sheet_dir, p)
  normalizePath(p, winslash = "/", mustWork = FALSE)
}

guess_format <- function(p) {
  if (dir.exists(p)) return("spaceranger")
  if (grepl("\\.rds$", p, ignore.case = TRUE)) return("rds")
  NA_character_
}

# Turn one sheet row into a full parameter list. Returns list(params, problems).
resolve_sample <- function(row, i, sheet_dir) {
  problems <- character(0)
  add <- function(...) problems <<- c(problems, sprintf(...))
  sid <- row$sample_id
  p <- DEFAULTS
  for (nm in names(DEFAULTS)) {
    if (!nm %in% names(row) || is.na(row[[nm]])) next
    v <- row[[nm]]
    if (nm %in% NUMERIC_PARAMS) {
      v <- suppressWarnings(as.numeric(v))
      if (is.na(v)) add("[%s] %s must be a number, got '%s'", sid, nm, row[[nm]])
    } else if (nm == "format") {
      v <- tolower(v)
    }
    p[[nm]] <- v
  }
  if (is.na(sid)) add("Row %d has no sample_id", i)
  p$sample_id <- sid
  p$path      <- resolve_path(row$path, sheet_dir)

  if (is.na(p$path) || !file.exists(p$path)) {
    add("[%s] path not found: %s", sid, row$path)
  } else if (p$format == "auto") {
    p$format <- guess_format(p$path)
    if (is.na(p$format)) add("[%s] can't guess the format of %s; fill in the format column", sid, p$path)
  }
  if (!is.na(p$format) && p$format != "auto" && !p$format %in% KNOWN_FORMATS)
    add("[%s] unknown format '%s' (use %s)", sid, p$format, paste(KNOWN_FORMATS, collapse = ", "))
  if (isTRUE(p$mito_floor > p$mito_ceiling))
    add("[%s] mito_floor is larger than mito_ceiling", sid)
  for (nm in c("nmads_low", "nmads_high", "nmads_mito"))
    if (isTRUE(p[[nm]] <= 0)) add("[%s] %s must be > 0", sid, nm)
  list(params = p, problems = problems)
}

# Read the whole sheet; stop with ALL problems listed if anything is wrong.
read_sample_sheet <- function(sample_sheet_path, sheet_dir) {
  sheet <- read.csv(sample_sheet_path, colClasses = "character", na.strings = c("", "NA"),
                    check.names = FALSE, strip.white = TRUE)
  names(sheet) <- tolower(trimws(names(sheet)))
  sheet <- sheet[rowSums(!is.na(sheet)) > 0, , drop = FALSE]

  missing_cols <- setdiff(c("sample_id", "path"), names(sheet))
  if (length(missing_cols)) stop("Sample sheet is missing required column(s): ", paste(missing_cols, collapse = ", "))

  ignored_cols <- setdiff(names(sheet), c("sample_id", "path", names(DEFAULTS)))
  if (length(ignored_cols)) log_msg("Ignoring extra column(s): %s", paste(ignored_cols, collapse = ", "))

  problems <- character(0)
  dups <- unique(sheet$sample_id[duplicated(sheet$sample_id)])
  if (length(dups)) problems <- c(problems, sprintf("Duplicate sample_id: %s", paste(dups, collapse = ", ")))

  resolved <- lapply(seq_len(nrow(sheet)), function(i) resolve_sample(sheet[i, , drop = FALSE], i, sheet_dir))
  problems <- c(problems, unlist(lapply(resolved, `[[`, "problems")))
  if (length(problems)) stop(paste(c("Sample sheet problems (nothing was run):", problems), collapse = "\n  - "))

  samples <- lapply(resolved, `[[`, "params")
  names(samples) <- vapply(samples, function(p) p$sample_id, "")
  samples
}

write_resolved_params <- function(samples, out_dir) {
  resolved <- do.call(rbind, lapply(samples, function(p) as.data.frame(
    p[c("sample_id", "path", names(DEFAULTS))], stringsAsFactors = FALSE)))
  write.csv(resolved, file.path(out_dir, "samplesheet_resolved.csv"), row.names = FALSE)
  log_msg("Sample sheet OK: %d sample(s) in %d dataset(s). Resolved parameters -> samplesheet_resolved.csv",
          length(samples), length(unique(resolved$dataset)))
}

## =============================================================================
## Readers: Space Ranger output (any layout) or a Seurat .rds
## =============================================================================
# Find exactly one file matching a pattern anywhere under `dir`. When several
# match (e.g. filtered_ and raw_ matrices), prefer paths containing "filtered".
find_file <- function(dir, pattern, what, required = TRUE) {
  f <- list.files(dir, pattern = pattern, full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  if (length(f) > 1 && any(grepl("filtered", f, ignore.case = TRUE))) f <- f[grepl("filtered", f, ignore.case = TRUE)]
  if (length(f) > 1) stop("More than one ", what, " file under ", dir, ": ", paste(basename(f), collapse = ", "),
                          ". Put each sample in its own folder.")
  if (length(f) == 0) {
    if (required) stop("No ", what, " file found under ", dir)
    return(NA_character_)
  }
  f
}

# Copy a file to `dest`, un-gzipping it if needed (GEO often gzips images and JSON)
copy_plain <- function(src, dest) {
  if (!grepl("\\.gz$", src, ignore.case = TRUE)) { file.copy(src, dest, overwrite = TRUE); return(invisible(dest)) }
  inp <- gzfile(src, "rb"); out <- file(dest, "wb")
  on.exit({ close(inp); close(out) })
  repeat { b <- readBin(inp, "raw", 1e7); if (!length(b)) break; writeBin(b, out) }
  invisible(dest)
}

# tissue_positions(.csv / _list.csv): Space Ranger v1 has no header, v2 has one.
read_positions <- function(path) {
  first <- readLines(if (grepl("\\.gz$", path)) gzfile(path) else path, n = 1)
  has_header <- grepl("barcode", first, ignore.case = TRUE)
  pos <- read.csv(path, header = has_header, stringsAsFactors = FALSE)
  pos <- pos[, 1:6]
  colnames(pos) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col")
  list(table = pos, has_header = has_header)
}

read_counts_spaceranger <- function(dir, gene_column) {
  h5 <- find_file(dir, "feature_bc_matrix\\.h5$", "count matrix .h5", required = FALSE)
  if (!is.na(h5)) {
    m <- Read10X_h5(h5, use.names = gene_column == 2)
    if (is.list(m)) m <- m[["Gene Expression"]]
    return(m)
  }
  mtx  <- find_file(dir, "matrix\\.mtx(\\.gz)?$", "matrix.mtx")
  bcs  <- find_file(dir, "barcodes\\.tsv(\\.gz)?$", "barcodes.tsv")
  feat <- find_file(dir, "(features|genes)\\.tsv(\\.gz)?$", "features.tsv", required = FALSE)
  if (is.na(feat)) {   # some GEO series share one features file across samples, one folder up
    feat <- find_file(dirname(dir), "(features|genes)\\.tsv(\\.gz)?$", "features.tsv (sample or parent folder)")
  }
  ft <- read.delim(feat, header = FALSE, stringsAsFactors = FALSE)
  m <- ReadMtx(mtx = mtx, cells = bcs, features = feat, feature.column = min(gene_column, ncol(ft)))
  if (ncol(ft) >= 3 && any(ft[[3]] == "Gene Expression") && any(ft[[3]] != "Gene Expression"))
    m <- m[ft[[3]] == "Gene Expression", , drop = FALSE]
  m
}

# Build a Seurat object with counts, the tissue image, and spot positions in the metadata
load_spaceranger <- function(dir, p) {
  counts <- read_counts_spaceranger(dir, p$gene_column)
  rownames(counts) <- make.unique(gsub("_", "-", rownames(counts)))   # Seurat forbids "_" in gene names
  pos_f  <- find_file(dir, "tissue_positions(_list)?\\.csv(\\.gz)?$", "tissue_positions.csv")
  sf_f   <- find_file(dir, "scalefactors_json\\.json(\\.gz)?$", "scalefactors_json.json")
  img_f  <- find_file(dir, "tissue_lowres_image\\.png(\\.gz)?$", "tissue_lowres_image.png", required = FALSE)
  img_name <- "tissue_lowres_image.png"
  if (is.na(img_f)) {
    img_f <- find_file(dir, "tissue_hires_image\\.png(\\.gz)?$", "tissue image (lowres or hires .png)")
    img_name <- "tissue_hires_image.png"
    log_msg("[%s] no lowres image; using the hires image", p$sample_id)
  }

  # Seurat's image reader wants standard file names in one folder, so stage them in a temp folder
  pos <- read_positions(pos_f)
  tmp <- file.path(tempdir(), paste0("visium_", p$sample_id))
  dir.create(tmp, showWarnings = FALSE, recursive = TRUE)
  copy_plain(pos_f, file.path(tmp, if (pos$has_header) "tissue_positions.csv" else "tissue_positions_list.csv"))
  copy_plain(sf_f,  file.path(tmp, "scalefactors_json.json"))
  copy_plain(img_f, file.path(tmp, img_name))
  image <- Read10X_Image(image.dir = tmp, image.name = img_name, filter.matrix = TRUE)

  in_tissue <- pos$table$barcode[pos$table$in_tissue == 1]
  keep <- intersect(colnames(counts), in_tissue)
  if (!length(keep)) stop("No barcodes shared between the count matrix and the in-tissue spots in ", basename(pos_f))
  if (ncol(counts) > length(keep))
    log_msg("[%s] %d of %d barcodes are in tissue; keeping those", p$sample_id, length(keep), ncol(counts))
  counts <- counts[, keep, drop = FALSE]

  obj <- CreateSeuratObject(counts = counts, assay = "Spatial")
  image <- image[Cells(obj)]
  DefaultAssay(image) <- "Spatial"
  obj[["slice1"]] <- image
  posm <- pos$table[match(Cells(obj), pos$table$barcode), c("array_row", "array_col", "pxl_row", "pxl_col")]
  rownames(posm) <- Cells(obj)
  AddMetaData(obj, posm)
}

load_rds <- function(path) {
  obj <- readRDS(path)
  if (!inherits(obj, "Seurat")) stop("The .rds holds a ", class(obj)[1], "; expected a Seurat object")
  if (!"Spatial" %in% names(obj@assays)) stop("The Seurat object has no 'Spatial' assay")
  obj[["Spatial"]] <- JoinLayers(obj[["Spatial"]])
  obj
}

# Spot positions for plotting and deconvolution: x = column, y = row (pixels, full-res image)
spot_xy <- function(obj) {
  md <- obj[[]]
  if (all(c("pxl_row", "pxl_col") %in% names(md))) {
    return(data.frame(cell = Cells(obj), x = md$pxl_col, y = md$pxl_row))
  }
  tc <- GetTissueCoordinates(obj)
  if (all(c("imagerow", "imagecol") %in% names(tc))) return(data.frame(cell = rownames(tc), x = tc$imagecol, y = tc$imagerow))
  data.frame(cell = rownames(tc), x = tc$y, y = tc$x)   # Seurat v5 VisiumV2 stores x = row, y = column
}

## =============================================================================
## QC steps
## =============================================================================
mad_bounds <- function(x, nlow, nhigh) {
  med <- median(x); dev <- mad(x)
  c(low = med - nlow * dev, high = med + nhigh * dev)
}
na_inf <- function(x) if (is.na(x)) Inf else x

# Step 1-2: read the sample, in-tissue spots only, spot names prefixed with the sample ID
load_sample <- function(p) {
  log_msg("[%s] reading %s (%s)", p$sample_id, p$path, p$format)
  obj <- if (p$format == "rds") load_rds(p$path) else load_spaceranger(p$path, p)
  obj <- RenameCells(obj, add.cell.id = p$sample_id)
  obj$sample  <- p$sample_id
  obj$dataset <- p$dataset
  log_msg("[%s] %d genes x %d in-tissue spots", p$sample_id, nrow(obj), ncol(obj))
  obj
}

# Step 3: per-spot QC metrics (added to the object's metadata)
compute_qc_metrics <- function(obj, p) {
  m <- LayerData(obj, assay = "Spatial", layer = "counts")
  mt_genes <- grep(p$mito_pattern, rownames(m), value = TRUE)
  if (length(mt_genes) == 0) log_msg("[%s] WARNING: no genes match the mito pattern '%s' (Ensembl IDs instead of gene symbols? check gene_column)", p$sample_id, p$mito_pattern)
  n_counts <- Matrix::colSums(m)
  pct <- function(genes) 100 * Matrix::colSums(m[genes, , drop = FALSE]) / pmax(n_counts, 1)
  obj$percent.mt   <- pct(mt_genes)
  obj$percent.ribo <- pct(grep(HUMAN_GENES$ribo, rownames(m), value = TRUE))
  obj$percent.hb   <- pct(grep(HUMAN_GENES$hb,   rownames(m), value = TRUE))
  attr(obj, "n_mito_genes") <- length(mt_genes)
  obj
}

# Step 4a: this sample's MAD-based cutoffs (same rules as scrna_qc.R)
compute_thresholds <- function(md, p, n_mito_genes) {
  lc <- mad_bounds(log10(pmax(md$nCount_Spatial, 1)),   p$nmads_low, p$nmads_high)
  lf <- mad_bounds(log10(pmax(md$nFeature_Spatial, 1)), p$nmads_low, p$nmads_high)
  mt_mad <- median(md$percent.mt) + p$nmads_mito * mad(md$percent.mt)
  data.frame(sample = p$sample_id, dataset = p$dataset, spots_in = nrow(md),
             min_counts   = max(p$min_counts,   floor(10^lc[["low"]])),
             max_counts   = min(na_inf(p$max_counts),   ceiling(10^lc[["high"]])),
             min_features = max(p$min_features, floor(10^lf[["low"]])),
             max_features = min(na_inf(p$max_features), ceiling(10^lf[["high"]])),
             max_mito_mad = round(mt_mad, 2),
             max_mito     = round(min(max(mt_mad, p$mito_floor), p$mito_ceiling), 2),
             median_counts = median(md$nCount_Spatial), median_features = median(md$nFeature_Spatial),
             median_mito = round(median(md$percent.mt), 2), n_mito_genes = n_mito_genes)
}

# Step 4b: flag spots outside the cutoffs; returns the metadata with qc_reason and keep
apply_spot_filters <- function(md, th, p) {
  checks <- list(
    low_counts    = md$nCount_Spatial   <  th$min_counts,
    high_counts   = md$nCount_Spatial   >  th$max_counts,
    low_features  = md$nFeature_Spatial <  th$min_features,
    high_features = md$nFeature_Spatial >  th$max_features,
    mito          = md$percent.mt       >  th$max_mito)
  md$qc_reason <- "pass"
  for (nm in rev(names(checks))) md$qc_reason[checks[[nm]]] <- nm   # first failing check wins
  md$keep <- md$qc_reason == "pass"
  log_msg("[%s] spot filters: %d of %d spots pass (mito cutoff %.1f%%, min genes %d, min UMIs %d)",
          p$sample_id, sum(md$keep), nrow(md), th$max_mito, as.integer(th$min_features), as.integer(th$min_counts))
  md
}

# Runs steps 1-4 for one sample; returns the full object (for plots), the kept object, and tables
qc_one_sample <- function(p) {
  set.seed(p$seed)
  obj <- load_sample(p)                                         # 1-2. read, in-tissue spots
  obj <- compute_qc_metrics(obj, p)                             # 3. metrics
  md  <- obj[[]]
  th  <- compute_thresholds(md, p, attr(obj, "n_mito_genes"))  # 4a. cutoffs
  md  <- apply_spot_filters(md, th, p)                          # 4b. filter
  obj$qc_reason <- md$qc_reason
  xy <- spot_xy(obj)
  md$x <- xy$x; md$y <- xy$y
  md$cell <- rownames(md)
  list(sample_id = p$sample_id, dataset = p$dataset, obj = obj, meta = md, thresholds = th,
       min_spots_gene = p$min_spots_gene)
}

## =============================================================================
## Outputs
## =============================================================================
# Step 5: drop rare genes and save the kept spots (with the image) for this sample
save_sample <- function(r, ds_dir) {
  dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)
  kept <- r$meta$cell[r$meta$keep]
  sub <- subset(r$obj, cells = kept)
  m <- LayerData(sub, assay = "Spatial", layer = "counts")
  genes_keep <- rownames(m)[Matrix::rowSums(m > 0) >= r$min_spots_gene]
  sub <- subset(sub, features = genes_keep)
  sub$qc_reason <- NULL
  Idents(sub) <- "sample"
  f <- file.path(ds_dir, paste0(r$sample_id, "_qc.rds"))
  saveRDS(sub, f)
  log_msg("[%s] saved %s: %d spots, %d genes (in >= %d spots)", r$sample_id, f, ncol(sub), nrow(sub), r$min_spots_gene)
}

summarise_sample <- function(r) {
  m <- r$meta
  n <- function(reason) sum(m$qc_reason == reason)
  data.frame(dataset = r$dataset, sample = r$sample_id, status = "ok", spots_in = nrow(m),
             low_counts = n("low_counts"), high_counts = n("high_counts"),
             low_features = n("low_features"), high_features = n("high_features"),
             high_mito = n("mito"), spots_kept = sum(m$keep), pct_kept = round(100 * mean(m$keep), 1))
}

# One spatial map per sample, combined into one figure
spatial_panels <- function(meta, fill_fun, title) {
  plots <- lapply(split(meta, meta$sample), function(d) {
    fill_fun(ggplot(d, aes(x, y))) + scale_y_reverse() + coord_fixed() +
      labs(title = d$sample[1], x = NULL, y = NULL) +
      theme_void(base_size = 10) + theme(legend.position = "right", plot.title = element_text(hjust = 0.5))
  })
  wrap_plots(plots, ncol = min(length(plots), 3)) + plot_annotation(title = title)
}

make_plots <- function(meta, thresholds, plot_dir) {
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  theme_set(theme_bw(base_size = 11))
  ns <- length(unique(meta$sample)); nc <- min(ns, 3); nr <- ceiling(ns / nc)
  w <- max(6, 3.2 * nc)

  # 01: violins before vs after, cutoffs dashed
  vdf <- rbind(transform(meta, stage = "before"), transform(meta[meta$keep, ], stage = "after"))
  vdf$stage <- factor(vdf$stage, levels = c("before", "after"))
  lines <- function(cols) do.call(rbind, lapply(cols, function(cl) data.frame(sample = thresholds$sample, y = thresholds[[cl]])))
  vln <- function(var, cols, log_y, title) {
    ld <- lines(cols); ld <- ld[is.finite(ld$y), ]
    g <- ggplot(vdf, aes(stage, .data[[var]], fill = stage)) +
      geom_violin(scale = "width", linewidth = 0.3) + geom_jitter(width = 0.15, size = 0.15, alpha = 0.25) +
      facet_wrap(~ sample, ncol = nc) +
      scale_fill_manual(values = c(before = "grey75", after = "#4C9BE8"), guide = "none") +
      labs(title = title, x = NULL, y = NULL)
    if (nrow(ld)) g <- g + geom_hline(data = ld, aes(yintercept = y), linetype = "dashed", colour = "red3")
    if (log_y) g <- g + scale_y_log10()
    g
  }
  p1 <- vln("nCount_Spatial", c("min_counts", "max_counts"), TRUE, "UMIs per spot (log)") /
        vln("nFeature_Spatial", c("min_features", "max_features"), TRUE, "Genes per spot (log)") /
        vln("percent.mt", "max_mito", FALSE, "% mitochondrial")
  ggsave(file.path(plot_dir, "01_violin_before_after.png"), p1, width = w, height = 3 * nr * 3, dpi = 150, limitsize = FALSE)

  # 02: UMIs vs genes, coloured by QC result
  meta$qc_reason <- factor(meta$qc_reason, levels = names(REASON_COLOURS))
  p2 <- ggplot(meta, aes(nCount_Spatial, nFeature_Spatial, colour = qc_reason)) +
    geom_point(size = 0.5, alpha = 0.7) + scale_x_log10() + scale_y_log10() +
    facet_wrap(~ sample, ncol = nc) + scale_colour_manual(values = REASON_COLOURS, drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    labs(title = "UMIs vs genes per spot, coloured by QC result", x = "UMIs (log)", y = "Genes (log)", colour = "QC result")
  ggsave(file.path(plot_dir, "02_counts_vs_genes_by_reason.png"), p2, width = w + 2, height = 3.5 * nr + 0.5, dpi = 150, limitsize = FALSE)

  # 03: WHERE the removed spots are on the tissue
  p3 <- spatial_panels(meta, function(g) g + geom_point(aes(colour = qc_reason), size = 0.6) +
                         scale_colour_manual(values = REASON_COLOURS, drop = FALSE, name = "QC result"),
                       "Removed spots on the tissue (grey = kept)")
  ggsave(file.path(plot_dir, "03_spatial_removed_spots.png"), p3, width = 4.5 * nc + 1.5, height = 4.2 * nr + 0.6, dpi = 150, limitsize = FALSE)

  # 04-05: UMIs and % mito across the tissue (all in-tissue spots, before filtering)
  p4 <- spatial_panels(meta, function(g) g + geom_point(aes(colour = log10(pmax(nCount_Spatial, 1))), size = 0.6) +
                         scale_colour_viridis_c(name = "log10 UMIs"), "UMIs per spot across the tissue")
  ggsave(file.path(plot_dir, "04_spatial_umis.png"), p4, width = 4.5 * nc + 1.5, height = 4.2 * nr + 0.6, dpi = 150, limitsize = FALSE)
  p5 <- spatial_panels(meta, function(g) g + geom_point(aes(colour = percent.mt), size = 0.6) +
                         scale_colour_viridis_c(option = "magma", name = "% mito"), "% mitochondrial per spot across the tissue")
  ggsave(file.path(plot_dir, "05_spatial_mito.png"), p5, width = 4.5 * nc + 1.5, height = 4.2 * nr + 0.6, dpi = 150, limitsize = FALSE)
}

write_dataset_outputs <- function(rs, ds_dir) {
  dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)
  keep_cols <- c("cell", "sample", "dataset", "nCount_Spatial", "nFeature_Spatial", "percent.mt", "percent.ribo",
                 "percent.hb", "array_row", "array_col", "pxl_row", "pxl_col", "x", "y", "qc_reason", "keep")
  meta_all <- do.call(rbind, lapply(rs, function(r) r$meta[, intersect(keep_cols, names(r$meta))]))
  th_all   <- do.call(rbind, lapply(rs, `[[`, "thresholds"))
  sm       <- do.call(rbind, lapply(rs, summarise_sample))
  write.csv(th_all,   file.path(ds_dir, "qc_thresholds.csv"), row.names = FALSE)
  write.csv(meta_all, file.path(ds_dir, "spot_qc_metadata.csv"), row.names = FALSE)
  tryCatch(make_plots(meta_all, th_all, file.path(ds_dir, "plots")),
           error = function(e) log_msg("[%s] WARNING: plotting failed (%s)", basename(ds_dir), conditionMessage(e)))
  sm
}

write_overview <- function(summaries, failures, samples, out_dir) {
  for (sid in names(failures)) {
    summaries[[sid]] <- data.frame(dataset = samples[[sid]]$dataset, sample = sid,
                                   status = paste("failed:", failures[[sid]]))
  }
  if (!length(summaries)) return(invisible(NULL))
  cols <- unique(unlist(lapply(summaries, names)))
  fill <- function(d) { for (cl in setdiff(cols, names(d))) d[[cl]] <- NA; d[cols] }
  overview <- do.call(rbind, lapply(summaries, fill))
  rownames(overview) <- NULL
  write.csv(overview, file.path(out_dir, "all_samples_summary.csv"), row.names = FALSE)
  print(overview[, intersect(c("dataset", "sample", "status", "spots_in", "spots_kept", "pct_kept"), cols)], row.names = FALSE)
  invisible(overview)
}

## =============================================================================
## main(): the whole pipeline, in order
## =============================================================================
main <- function(sample_sheet_path = SAMPLE_SHEET_PATH, out_dir = OUT_DIR) {
  if (!file.exists(sample_sheet_path)) stop("Sample sheet not found: ", sample_sheet_path,
                                            "\n(current working directory: ", getwd(), ")")
  sheet_dir <- dirname(normalizePath(sample_sheet_path, winslash = "/"))
  if (is.null(out_dir)) out_dir <- file.path(sheet_dir, "qc_output")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  start_log(out_dir)

  samples <- read_sample_sheet(sample_sheet_path, sheet_dir)
  write_resolved_params(samples, out_dir)

  results  <- list()
  failures <- list()
  for (p in samples) {
    r <- tryCatch({
      r <- qc_one_sample(p)
      save_sample(r, file.path(out_dir, p$dataset))   # save now, then drop the big object
      r$obj <- NULL
      r
    }, error = function(e) {
      log_msg("[%s] FAILED: %s", p$sample_id, conditionMessage(e)); conditionMessage(e) })
    if (is.character(r)) failures[[p$sample_id]] <- r else results[[p$sample_id]] <- r
  }

  summaries <- list()
  for (ds in unique(vapply(results, function(r) r$dataset, ""))) {
    rs <- Filter(function(r) r$dataset == ds, results)
    summaries[[ds]] <- write_dataset_outputs(rs, file.path(out_dir, ds))
  }

  write_overview(summaries, failures, samples, out_dir)
  log_msg("Done: %d sample(s) OK, %d failed. Outputs in %s", length(results), length(failures), out_dir)
  invisible(list(results = results, failures = failures, out_dir = out_dir))
}

## =============================================================================
## Run
## =============================================================================
args <- commandArgs(trailingOnly = TRUE)
run <- main(sample_sheet_path = if (length(args) >= 1) args[1] else SAMPLE_SHEET_PATH,
            out_dir           = if (length(args) >= 2) args[2] else OUT_DIR)
if (length(run$failures) && !interactive()) quit(status = 1)
