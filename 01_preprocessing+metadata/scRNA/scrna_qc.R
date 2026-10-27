#!/usr/bin/env Rscript
###############################################################################
# scrna_qc.R — sample-sheet-driven scRNA-seq quality control (human data)
#
# Usage
#   Rscript scrna_qc.R <samplesheet.csv> [out_dir]
#   ...or open in RStudio, set SAMPLE_SHEET_PATH below, and click Source.
#   (After sourcing once, you can also call main("path/to/sheet.csv") yourself.)
#
# The sample sheet has ONE ROW PER SAMPLE. Required columns: sample_id, path.
# Every other column is optional. A blank cell means "use the default" listed
# in DEFAULTS below. Every sample is QC'd and saved as its OWN Seurat object;
# nothing is merged. `dataset` only decides which output folder a sample's
# files go in (and which samples share a summary table and plot panel).
#
# Per sample, the pipeline:
#   1. reads the counts (10x folder, 10x .h5, dense text table, or .rds)
#   2. calls cells with emptyDrops if the matrix is unfiltered (is_filtered = FALSE)
#   3. computes QC metrics (UMIs, genes, % mito; % ribo and % hemoglobin are reported only)
#   4. removes outliers using per-sample MAD thresholds
#   5. removes doublets with scDblFinder
#   6. estimates ambient RNA with DecontX (uses raw droplets as background if given)
#   7. drops rare genes and saves <sample_id>_qc.rds
# Per dataset folder, it also writes summary tables and QC plots (one panel per sample).
#
# Layout of this file
#   Settings and defaults  -> what you normally edit
#   Helper functions       -> logging, sample sheet, readers, QC steps, outputs
#   main()                 -> runs the whole pipeline, in order
#   Last lines             -> call main()
###############################################################################

## =============================================================================
## Settings (used when no command-line arguments are given)
## =============================================================================
SAMPLE_SHEET_PATH <- "C:/ongkeko lab works/bladder cancer/data/GSE130001/samplesheet_GSE130001.csv"
OUT_DIR           <- NULL   # NULL = "<folder of the sample sheet>/qc_output"

## =============================================================================
## Defaults for every sample-sheet column
## =============================================================================
DEFAULTS <- list(
  # what the sample is
  dataset           = "dataset1", # output folder name; samples are never merged
  format            = "auto",     # auto | 10x | 10x_h5 | dense | dense_t | rds
  gene_column       = 2,          # 10x features file: 1 = Ensembl ID, 2 = gene symbol
  mito_pattern      = "^MT-",     # regex for mitochondrial genes (human gene symbols)
  is_filtered       = TRUE,       # FALSE = raw droplet matrix; cells called with emptyDrops
  raw_path          = NA,         # optional raw (unfiltered) matrix -> DecontX background
  # cell filters
  nmads_low         = 3,          # MADs below the median (log scale) for UMIs/genes
  nmads_high        = 5,          # MADs above the median (log scale) for UMIs/genes
  nmads_mito        = 3,          # MADs above the median for % mito
  min_counts        = 500,        # hard floor on UMIs (always applied)
  min_features      = 200,        # hard floor on genes (always applied)
  max_counts        = NA,         # optional hard cap on UMIs
  max_features      = NA,         # optional hard cap on genes
  mito_floor        = 10,         # MAD mito cutoff is never stricter than this
  mito_ceiling      = 25,         # MAD mito cutoff is never looser than this
  # doublets and ambient RNA
  run_doublets      = TRUE,
  doublet_rate      = NA,         # blank = scDblFinder default (~1% per 1,000 cells)
  run_decontx       = TRUE,
  max_contamination = NA,         # e.g. 0.5 removes cells >50% ambient; blank = keep all
  # cell calling (only when is_filtered = FALSE)
  emptydrops_lower  = 100,        # droplets with <= this many UMIs define "empty"
  emptydrops_fdr    = 0.001,
  # genes and reproducibility
  min_cells_gene    = 3,          # keep genes seen in >= this many of the sample's kept cells
  seed              = 1234
)
LOGICAL_PARAMS   <- c("is_filtered", "run_doublets", "run_decontx")
CHARACTER_PARAMS <- c("dataset", "format", "mito_pattern", "raw_path")
NUMERIC_PARAMS   <- setdiff(names(DEFAULTS), c(LOGICAL_PARAMS, CHARACTER_PARAMS))
KNOWN_FORMATS    <- c("10x", "10x_h5", "dense", "dense_t", "rds")

# Human gene-symbol patterns for ribosomal and hemoglobin genes (reported only; never used to filter)
HUMAN_GENES <- list(
  ribo = "^RP[SL][0-9]+[AXYL]?[0-9]*$|^RPLP[0-2]$|^RPSA$",   # RPS/RPL proteins, not e.g. RPS6KA1
  hb   = "^HB[ABDEGMQZ][0-9]?$"                              # HBA1/2, HBB, HBD, ...; not HBP1, HBEGF
)

# Plot colour for each reason a cell can be removed
REASON_COLOURS <- c(pass = "grey70", low_counts = "#E69F00", low_features = "#D55E00",
                    high_counts = "#56B4E9", high_features = "#0072B2", mito = "#CC79A7",
                    doublet = "#009E73", contamination = "#7F3C8D")

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
# convert T/F strings to bool
parse_logical <- function(x) {
  v <- toupper(trimws(x))
  if (v %in% c("TRUE", "T", "YES", "Y", "1")) return(TRUE)
  if (v %in% c("FALSE", "F", "NO", "N", "0")) return(FALSE)
  NA
}
# if selected auto mode, detect mode based on which files exist 
guess_format <- function(p) {
  if (dir.exists(p)) return("10x")
  if (grepl("\\.h5$", p, ignore.case = TRUE)) return("10x_h5")
  if (grepl("\\.rds$", p, ignore.case = TRUE)) return("rds")
  if (grepl("\\.(csv|tsv|txt)(\\.gz)?$", p, ignore.case = TRUE)) return("dense")
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
    if (nm %in% LOGICAL_PARAMS) {
      v <- parse_logical(v)
      if (is.na(v)) add("[%s] %s must be TRUE/FALSE, got '%s'", sid, nm, row[[nm]])
    } else if (nm %in% NUMERIC_PARAMS) {
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
  p$raw_path  <- resolve_path(p$raw_path, sheet_dir)

  if (is.na(p$path) || !file.exists(p$path)) {
    add("[%s] path not found: %s", sid, row$path)
  } else if (p$format == "auto") {
    p$format <- guess_format(p$path)
    if (is.na(p$format)) add("[%s] can't guess the format of %s; fill in the format column", sid, p$path)
  }
  if (!is.na(p$format) && p$format != "auto" && !p$format %in% KNOWN_FORMATS)
    add("[%s] unknown format '%s' (use %s)", sid, p$format, paste(KNOWN_FORMATS, collapse = ", "))
  if (!is.na(p$raw_path) && !file.exists(p$raw_path))
    add("[%s] raw_path not found: %s", sid, p$raw_path)
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
  sheet <- sheet[rowSums(!is.na(sheet)) > 0, , drop = FALSE]   # drop fully blank rows

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

# Stop early if the sheet asks for a step whose package isn't installed
check_packages <- function(samples) {
  need <- c()
  if (any(vapply(samples, function(p) p$run_doublets, TRUE))) need <- c(need, "SingleCellExperiment", "scDblFinder")
  if (any(vapply(samples, function(p) p$run_decontx, TRUE)))  need <- c(need, "celda")
  if (any(!vapply(samples, function(p) p$is_filtered, TRUE))) need <- c(need, "DropletUtils")
  if (any(vapply(samples, function(p) p$format == "10x_h5", TRUE))) need <- c(need, "hdf5r")
  missing_pkgs <- need[!vapply(need, requireNamespace, TRUE, quietly = TRUE)]
  if (length(missing_pkgs)) {
    stop("Missing packages: ", paste(missing_pkgs, collapse = ", "),
         "\nInstall with: BiocManager::install(c(", paste0('"', missing_pkgs, '"', collapse = ", "), "))",
         "\nor set run_doublets / run_decontx to FALSE in the sample sheet.")
  }
}

# Record the final value of every parameter for every sample (for your methods section)
write_resolved_params <- function(samples, out_dir) {
  resolved <- do.call(rbind, lapply(samples, function(p) as.data.frame(
    p[c("sample_id", "path", names(DEFAULTS))], stringsAsFactors = FALSE)))
  write.csv(resolved, file.path(out_dir, "samplesheet_resolved.csv"), row.names = FALSE)
  log_msg("Sample sheet OK: %d sample(s) in %d dataset(s). Resolved parameters -> samplesheet_resolved.csv",
          length(samples), length(unique(resolved$dataset)))
}

## =============================================================================
## Readers: every format ends up as a genes x cells sparse count matrix
## =============================================================================
as_dgc <- function(m) as(as(m, "CsparseMatrix"), "dgCMatrix")

read_10x_dir <- function(dir, gene_column) {
  find1 <- function(pattern, what) {
    f <- list.files(dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
    if (length(f) == 0) stop("No ", what, " file in ", dir)
    if (length(f) > 1) stop("More than one ", what, " file in ", dir,
                            " (put each sample in its own folder): ", paste(basename(f), collapse = ", "))
    f
  }
  mtx   <- find1("matrix\\.mtx(\\.gz)?$", "matrix.mtx")
  cells <- find1("barcodes\\.tsv(\\.gz)?$", "barcodes.tsv")
  feats <- find1("(features|genes)\\.tsv(\\.gz)?$", "features/genes.tsv")
  ft <- read.delim(feats, header = FALSE, stringsAsFactors = FALSE)
  gene_column <- min(gene_column, ncol(ft))
  m <- ReadMtx(mtx = mtx, cells = cells, features = feats, feature.column = gene_column)
  if (ncol(ft) >= 3 && any(ft[[3]] == "Gene Expression") && any(ft[[3]] != "Gene Expression")) {
    m <- m[ft[[3]] == "Gene Expression", , drop = FALSE]   # drop antibody / CRISPR features
  }
  m
}

read_10x_h5 <- function(path, gene_column) {
  m <- Read10X_h5(path, use.names = gene_column == 2)
  if (is.list(m)) m <- m[["Gene Expression"]]
  m
}

read_dense <- function(path, cells_in_rows = FALSE) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    df <- as.data.frame(data.table::fread(path, data.table = FALSE))
  } else {
    df <- read.table(path, header = TRUE, sep = if (grepl("\\.csv", path)) "," else "\t",
                     check.names = FALSE, stringsAsFactors = FALSE)
  }
  if (ncol(df) < 2 || nrow(df) < 2) stop("dense table needs a gene-name column plus at least one cell column: ", path)
  num <- vapply(df[, -1, drop = FALSE], is.numeric, TRUE)
  if (!all(num)) stop(sum(!num), " non-numeric column(s) after the first (e.g. '", names(num)[!num][1],
                      "'). Expected: first column = IDs, the rest = counts. Wrong orientation? Try dense_t.")
  ids <- as.character(df[[1]])
  m <- as.matrix(df[, -1, drop = FALSE])
  storage.mode(m) <- "double"
  rownames(m) <- make.unique(ids)
  if (cells_in_rows) m <- t(m)
  as_dgc(m)
}

read_rds <- function(path) {
  obj <- readRDS(path)
  if (inherits(obj, "Seurat")) {
    assay <- if ("RNA" %in% names(obj@assays)) "RNA" else SeuratObject::DefaultAssay(obj)
    a <- obj[[assay]]
    if (inherits(a, "Assay5")) {
      if (length(SeuratObject::Layers(a, search = "counts")) > 1) a <- SeuratObject::JoinLayers(a)
      return(SeuratObject::LayerData(a, layer = "counts"))
    }
    return(a@counts)
  }
  if (inherits(obj, "SingleCellExperiment")) return(SummarizedExperiment::assay(obj, "counts"))
  if (inherits(obj, c("matrix", "Matrix"))) return(obj)
  stop("The .rds holds a ", class(obj)[1], "; expected a Seurat object, SingleCellExperiment or count matrix")
}

read_counts <- function(path, format, gene_column) {
  m <- switch(format,
    "10x"     = read_10x_dir(path, gene_column),
    "10x_h5"  = read_10x_h5(path, gene_column),
    "dense"   = read_dense(path, cells_in_rows = FALSE),
    "dense_t" = read_dense(path, cells_in_rows = TRUE),
    "rds"     = read_rds(path))
  m <- as_dgc(m)
  rownames(m) <- make.unique(gsub("_", "-", rownames(m)))   # Seurat forbids "_" in gene names
  m
}

# Put a matrix's rows into a given gene order, filling absent genes with zeros
align_rows <- function(m, genes) {
  if (identical(rownames(m), genes)) return(m)
  idx <- match(genes, rownames(m))
  ok  <- which(!is.na(idx))
  mt  <- as(m[idx[ok], , drop = FALSE], "TsparseMatrix")
  sparseMatrix(i = ok[mt@i + 1], j = mt@j + 1, x = mt@x,
               dims = c(length(genes), ncol(m)), dimnames = list(genes, colnames(m)))
}

## =============================================================================
## QC steps (each one takes and returns plain objects, so you can run them alone)
## =============================================================================
mad_bounds <- function(x, nlow, nhigh) {
  med <- median(x); dev <- mad(x)
  c(low = med - nlow * dev, high = med + nhigh * dev)
}
na_inf <- function(x) if (is.na(x)) Inf else x

# Step 1: read the counts and prefix cell names with the sample ID
load_sample <- function(p) {
  log_msg("[%s] reading %s (%s)", p$sample_id, p$path, p$format)
  m <- read_counts(p$path, p$format, p$gene_column)
  colnames(m) <- paste0(p$sample_id, "_", colnames(m))
  log_msg("[%s] %d genes x %d barcodes", p$sample_id, nrow(m), ncol(m))
  xs <- head(m@x, 1e5)
  if (length(xs) && any(xs != round(xs))) {
    log_msg("[%s] WARNING: counts are not whole numbers; this looks like normalized data, and QC expects raw counts", p$sample_id)
  }
  m
}

# Step 2: call cells on unfiltered data, and collect empty droplets as DecontX background.
# Returns list(counts = cells only, background = empty droplets or NULL, n_barcodes).
call_cells <- function(m, p) {
  background <- NULL
  n_barcodes <- ncol(m)
  if (!p$is_filtered) {
    ed <- DropletUtils::emptyDrops(m, lower = p$emptydrops_lower)
    is_cell <- !is.na(ed$FDR) & ed$FDR <= p$emptydrops_fdr
    background <- m
    m <- m[, is_cell, drop = FALSE]
    log_msg("[%s] emptyDrops: %d cells called from %d barcodes (FDR <= %g)",
            p$sample_id, ncol(m), n_barcodes, p$emptydrops_fdr)
  } else if (!is.na(p$raw_path)) {
    background <- read_counts(p$raw_path, guess_format(p$raw_path), p$gene_column)
    colnames(background) <- paste0(p$sample_id, "_", colnames(background))
  }
  if (!is.null(background)) {
    background <- background[, !colnames(background) %in% colnames(m), drop = FALSE]
    background <- background[, Matrix::colSums(background) > 0, drop = FALSE]
    background <- align_rows(background, rownames(m))
  }
  list(counts = m, background = background, n_barcodes = n_barcodes)
}

# Step 3: per-cell QC metrics (one row per cell)
compute_qc_metrics <- function(m, p) {
  mt_genes <- grep(p$mito_pattern, rownames(m), value = TRUE)
  if (length(mt_genes) == 0) log_msg("[%s] WARNING: no genes match the mito pattern '%s' (Ensembl IDs instead of gene symbols? check gene_column)", p$sample_id, p$mito_pattern)
  n_counts   <- Matrix::colSums(m)
  n_features <- Matrix::colSums(m > 0)
  pct <- function(genes) 100 * Matrix::colSums(m[genes, , drop = FALSE]) / pmax(n_counts, 1)
  meta <- data.frame(
    cell = colnames(m), sample = p$sample_id, dataset = p$dataset,
    nCount_RNA = n_counts, nFeature_RNA = n_features,
    percent.mt   = pct(mt_genes),
    percent.ribo = pct(grep(HUMAN_GENES$ribo, rownames(m), value = TRUE)),
    percent.hb   = pct(grep(HUMAN_GENES$hb,   rownames(m), value = TRUE)),
    log10GenesPerUMI = log10(pmax(n_features, 1)) / log10(pmax(n_counts, 2)),
    row.names = colnames(m), stringsAsFactors = FALSE)
  attr(meta, "n_mito_genes") <- length(mt_genes)
  meta
}

# Step 4a: work out this sample's MAD-based cutoffs.
# UMIs and genes: log10 scale, median -/+ nmads_low/nmads_high MADs, never past the hard floors/caps.
# % mito: median + nmads_mito MADs, clamped to [mito_floor, mito_ceiling].
compute_thresholds <- function(meta, p, n_barcodes) {
  lc <- mad_bounds(log10(pmax(meta$nCount_RNA, 1)),   p$nmads_low, p$nmads_high)
  lf <- mad_bounds(log10(pmax(meta$nFeature_RNA, 1)), p$nmads_low, p$nmads_high)
  mt_mad <- median(meta$percent.mt) + p$nmads_mito * mad(meta$percent.mt)
  th <- list(min_counts   = max(p$min_counts,   floor(10^lc[["low"]])),
             max_counts   = min(na_inf(p$max_counts),   ceiling(10^lc[["high"]])),
             min_features = max(p$min_features, floor(10^lf[["low"]])),
             max_features = min(na_inf(p$max_features), ceiling(10^lf[["high"]])),
             max_mito_mad = round(mt_mad, 2),                                         # before clamping
             max_mito     = round(min(max(mt_mad, p$mito_floor), p$mito_ceiling), 2))  # cutoff used
  data.frame(sample = p$sample_id, dataset = p$dataset,
             barcodes_in = n_barcodes, cells_in = nrow(meta), th,
             median_counts = median(meta$nCount_RNA), median_features = median(meta$nFeature_RNA),
             median_mito = round(median(meta$percent.mt), 2), n_mito_genes = attr(meta, "n_mito_genes"))
}

# Step 4b: flag cells outside the cutoffs; adds qc_reason and keep columns
apply_cell_filters <- function(meta, th, p) {
  checks <- list(
    low_counts    = meta$nCount_RNA   <  th$min_counts,
    high_counts   = meta$nCount_RNA   >  th$max_counts,
    low_features  = meta$nFeature_RNA <  th$min_features,
    high_features = meta$nFeature_RNA >  th$max_features,
    mito          = meta$percent.mt   >  th$max_mito)
  meta$qc_reason <- "pass"
  for (nm in rev(names(checks))) meta$qc_reason[checks[[nm]]] <- nm   # first failing check wins
  meta$keep <- meta$qc_reason == "pass"
  log_msg("[%s] cell filters: %d of %d cells pass (mito cutoff %.1f%%, min genes %d, min UMIs %d)",
          p$sample_id, sum(meta$keep), nrow(meta), th$max_mito, as.integer(th$min_features), as.integer(th$min_counts))
  meta
}

# Step 5: doublets with scDblFinder, on raw counts of cells still kept
detect_doublets <- function(meta, m, p) {
  meta$scDblFinder.score <- NA_real_
  meta$scDblFinder.class <- NA_character_
  if (!p$run_doublets) return(meta)
  keep <- meta$cell[meta$keep]
  if (length(keep) < 50) {
    log_msg("[%s] WARNING: only %d cells; skipping doublet detection", p$sample_id, length(keep))
    return(meta)
  }
  res <- tryCatch({
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = m[, keep]))
    set.seed(p$seed)
    if (is.na(p$doublet_rate)) scDblFinder::scDblFinder(sce) else scDblFinder::scDblFinder(sce, dbr = p$doublet_rate)
  }, error = function(e) { log_msg("[%s] WARNING: scDblFinder failed (%s); no doublets removed", p$sample_id, conditionMessage(e)); NULL })
  if (is.null(res)) return(meta)
  meta[keep, "scDblFinder.score"] <- res$scDblFinder.score
  meta[keep, "scDblFinder.class"] <- as.character(res$scDblFinder.class)
  dbl <- keep[res$scDblFinder.class == "doublet"]
  meta[dbl, "qc_reason"] <- "doublet"
  meta[dbl, "keep"] <- FALSE
  log_msg("[%s] scDblFinder: %d doublets among %d cells", p$sample_id, length(dbl), length(keep))
  meta
}

# Step 6: ambient RNA with DecontX, on cells still kept. Returns list(meta, decontx_counts).
correct_ambient <- function(meta, m, background, p) {
  meta$decontX_contamination <- NA_real_
  if (!p$run_decontx) return(list(meta = meta, decontx = NULL))
  keep <- meta$cell[meta$keep]
  res <- tryCatch({
    if (is.null(background)) {
      celda::decontX(x = m[, keep], seed = p$seed, verbose = FALSE)
    } else {
      celda::decontX(x = m[, keep], background = background, seed = p$seed, verbose = FALSE)
    }
  }, error = function(e) { log_msg("[%s] WARNING: DecontX failed (%s); no ambient correction", p$sample_id, conditionMessage(e)); NULL })
  if (is.null(res)) return(list(meta = meta, decontx = NULL))
  meta[keep, "decontX_contamination"] <- res$contamination
  log_msg("[%s] DecontX: median contamination %.1f%%%s", p$sample_id, 100 * median(res$contamination),
          if (is.null(background)) "" else sprintf(" (background: %d empty droplets)", ncol(background)))
  if (!is.na(p$max_contamination)) {
    bad <- keep[res$contamination > p$max_contamination]
    meta[bad, "qc_reason"] <- "contamination"
    meta[bad, "keep"] <- FALSE
    log_msg("[%s] removed %d cells with contamination > %.2f", p$sample_id, length(bad), p$max_contamination)
  }
  list(meta = meta, decontx = res$decontXcounts)
}

# Runs steps 1-6 for one sample and returns everything needed to save it
qc_one_sample <- function(p) {
  set.seed(p$seed)
  raw     <- load_sample(p)                                     # 1. read
  called  <- call_cells(raw, p)                                 # 2. cell calling / background
  m       <- called$counts
  meta    <- compute_qc_metrics(m, p)                           # 3. metrics
  th      <- compute_thresholds(meta, p, called$n_barcodes)     # 4a. cutoffs
  meta    <- apply_cell_filters(meta, th, p)                    # 4b. filter
  meta    <- detect_doublets(meta, m, p)                        # 5. doublets
  amb     <- correct_ambient(meta, m, called$background, p)     # 6. ambient RNA
  meta    <- amb$meta

  kept <- meta$cell[meta$keep]
  list(sample_id = p$sample_id, dataset = p$dataset, meta = meta, thresholds = th,
       counts  = m[, kept, drop = FALSE],
       decontx = if (is.null(amb$decontx)) NULL else amb$decontx[, kept, drop = FALSE],
       min_cells_gene = p$min_cells_gene)
}

## =============================================================================
## Outputs: Seurat objects, summary tables, plots
## =============================================================================
# Step 7: drop rare genes and save one Seurat object for this sample
save_sample <- function(r, ds_dir) {
  sid <- r$sample_id
  genes_keep <- rownames(r$counts)[Matrix::rowSums(r$counts > 0) >= r$min_cells_gene]
  md <- r$meta[colnames(r$counts), c("sample", "dataset", "percent.mt", "percent.ribo", "percent.hb",
                                      "log10GenesPerUMI", "scDblFinder.score", "decontX_contamination")]
  seu <- CreateSeuratObject(counts = r$counts[genes_keep, , drop = FALSE], project = sid,
                            meta.data = md, min.cells = 0, min.features = 0)
  seu$orig.ident <- sid
  Idents(seu) <- "sample"
  if (!is.null(r$decontx)) seu[["decontX"]] <- CreateAssayObject(counts = r$decontx[genes_keep, colnames(seu), drop = FALSE])
  f <- file.path(ds_dir, paste0(sid, "_qc.rds"))
  saveRDS(seu, f)
  log_msg("[%s] saved %s: %d cells, %d genes (in >= %d cells)", sid, f, ncol(seu), nrow(seu), r$min_cells_gene)
}

# One summary row per sample: cells in, removed by each reason, kept
summarise_sample <- function(r) {
  m <- r$meta
  n <- function(reason) sum(m$qc_reason == reason)
  data.frame(dataset = r$dataset, sample = r$sample_id, status = "ok",
             barcodes_in = r$thresholds$barcodes_in, cells_in = nrow(m),
             low_counts = n("low_counts"), high_counts = n("high_counts"),
             low_features = n("low_features"), high_features = n("high_features"),
             high_mito = n("mito"), doublets = n("doublet"),
             high_contamination = n("contamination"), cells_kept = sum(m$keep),
             pct_kept = round(100 * mean(m$keep), 1),
             median_decontX_pct = round(100 * median(m$decontX_contamination, na.rm = TRUE), 1))
}

make_plots <- function(meta, thresholds, plot_dir) {
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  theme_set(theme_bw(base_size = 11))
  ns <- length(unique(meta$sample)); nc <- min(ns, 4); nr <- ceiling(ns / nc)
  w <- max(6, 2.6 * nc)

  vdf <- rbind(transform(meta, stage = "before"), transform(meta[meta$keep, ], stage = "after"))
  vdf$stage <- factor(vdf$stage, levels = c("before", "after"))
  lines <- function(cols) do.call(rbind, lapply(cols, function(cl)
    data.frame(sample = thresholds$sample, y = thresholds[[cl]])))
  vln <- function(var, cols, log_y, title) {
    ld <- lines(cols); ld <- ld[is.finite(ld$y), ]
    g <- ggplot(vdf, aes(stage, .data[[var]], fill = stage)) +
      geom_violin(scale = "width", linewidth = 0.3) +
      geom_jitter(width = 0.15, size = 0.15, alpha = 0.25) +
      facet_wrap(~ sample, ncol = nc) +
      scale_fill_manual(values = c(before = "grey75", after = "#4C9BE8"), guide = "none") +
      labs(title = title, x = NULL, y = NULL)
    if (nrow(ld)) g <- g + geom_hline(data = ld, aes(yintercept = y), linetype = "dashed", colour = "red3")
    if (log_y) g <- g + scale_y_log10()
    g
  }
  p1 <- vln("nCount_RNA", c("min_counts", "max_counts"), TRUE, "UMIs per cell (log)") /
        vln("nFeature_RNA", c("min_features", "max_features"), TRUE, "Genes per cell (log)") /
        vln("percent.mt", "max_mito", FALSE, "% mitochondrial")
  ggsave(file.path(plot_dir, "01_violin_before_after.png"), p1, width = w, height = 3 * nr * 3, dpi = 150, limitsize = FALSE)

  meta$qc_reason <- factor(meta$qc_reason, levels = names(REASON_COLOURS))
  p2 <- ggplot(meta, aes(nCount_RNA, nFeature_RNA, colour = qc_reason)) +
    geom_point(size = 0.4, alpha = 0.7) + scale_x_log10() + scale_y_log10() +
    facet_wrap(~ sample, ncol = nc) +
    scale_colour_manual(values = REASON_COLOURS, drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    labs(title = "UMIs vs genes per cell, coloured by QC result", x = "UMIs (log)", y = "Genes (log)", colour = "QC result")
  ggsave(file.path(plot_dir, "02_counts_vs_genes_by_reason.png"), p2, width = w + 2, height = 3.5 * nr + 0.5, dpi = 150, limitsize = FALSE)

  p3 <- ggplot(meta, aes(nCount_RNA, percent.mt)) + geom_point(size = 0.4, alpha = 0.5) +
    geom_hline(data = thresholds, aes(yintercept = max_mito), linetype = "dashed", colour = "red3") +
    geom_vline(data = thresholds, aes(xintercept = min_counts), linetype = "dashed", colour = "red3") +
    scale_x_log10() + facet_wrap(~ sample, ncol = nc) +
    labs(title = "UMIs vs % mitochondrial (dashed = cutoffs)", x = "UMIs (log)", y = "% mito")
  ggsave(file.path(plot_dir, "03_counts_vs_mito.png"), p3, width = w, height = 3.5 * nr + 0.5, dpi = 150, limitsize = FALSE)

  d <- meta[!is.na(meta$scDblFinder.score), ]
  if (nrow(d)) {
    p4 <- ggplot(d, aes(scDblFinder.score, fill = scDblFinder.class)) + geom_histogram(bins = 50) +
      facet_wrap(~ sample, ncol = nc, scales = "free_y") +
      scale_fill_manual(values = c(singlet = "grey70", doublet = "#009E73")) +
      labs(title = "scDblFinder doublet scores", x = "doublet score", y = "cells", fill = NULL)
    ggsave(file.path(plot_dir, "04_doublet_scores.png"), p4, width = w, height = 3 * nr + 0.5, dpi = 150, limitsize = FALSE)
  }
  d <- meta[!is.na(meta$decontX_contamination), ]
  if (nrow(d)) {
    p5 <- ggplot(d, aes(decontX_contamination)) + geom_histogram(bins = 50, fill = "#4C9BE8") +
      facet_wrap(~ sample, ncol = nc, scales = "free_y") +
      labs(title = "DecontX estimated ambient contamination", x = "fraction of counts from ambient RNA", y = "cells")
    ggsave(file.path(plot_dir, "05_decontX_contamination.png"), p5, width = w, height = 3 * nr + 0.5, dpi = 150, limitsize = FALSE)
  }
}

# Everything for one dataset folder: tables, one .rds per sample, plots. Returns the summary rows.
write_dataset_outputs <- function(rs, ds_dir) {
  dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)
  meta_all <- do.call(rbind, lapply(rs, `[[`, "meta"))
  th_all   <- do.call(rbind, lapply(rs, `[[`, "thresholds"))
  sm       <- do.call(rbind, lapply(rs, summarise_sample))
  rownames(meta_all) <- meta_all$cell
  write.csv(sm,       file.path(ds_dir, "qc_summary.csv"), row.names = FALSE)
  write.csv(th_all,   file.path(ds_dir, "qc_thresholds.csv"), row.names = FALSE)
  write.csv(meta_all, file.path(ds_dir, "cell_qc_metadata.csv"), row.names = FALSE)

  for (r in rs) {
    tryCatch(save_sample(r, ds_dir),
             error = function(e) log_msg("[%s] WARNING: saving failed (%s)", r$sample_id, conditionMessage(e)))
  }
  tryCatch(make_plots(meta_all, th_all, file.path(ds_dir, "plots")),
           error = function(e) log_msg("[%s] WARNING: plotting failed (%s)", basename(ds_dir), conditionMessage(e)))
  sm
}

# Top-level table covering every sample, including ones that failed
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
  print(overview[, intersect(c("dataset", "sample", "status", "cells_in", "cells_kept", "pct_kept"), cols)], row.names = FALSE)
  invisible(overview)
}

## =============================================================================
## main(): the whole pipeline, in order
## =============================================================================
main <- function(sample_sheet_path = SAMPLE_SHEET_PATH, out_dir = OUT_DIR) {
  # 1. Set up the output folder and log
  if (!file.exists(sample_sheet_path)) stop("Sample sheet not found: ", sample_sheet_path)
  sheet_dir <- dirname(normalizePath(sample_sheet_path, winslash = "/"))
  if (is.null(out_dir)) out_dir <- file.path(sheet_dir, "qc_output")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  start_log(out_dir)

  # 2. Read and check the sample sheet (stops here if anything is wrong)
  samples <- read_sample_sheet(sample_sheet_path, sheet_dir)
  check_packages(samples)
  write_resolved_params(samples, out_dir)

  # 3. QC each sample; a failure in one sample doesn't stop the others
  results  <- list()
  failures <- list()
  for (p in samples) {
    r <- tryCatch(qc_one_sample(p), error = function(e) {
      log_msg("[%s] FAILED: %s", p$sample_id, conditionMessage(e)); conditionMessage(e) })
    if (is.character(r)) failures[[p$sample_id]] <- r else results[[p$sample_id]] <- r
  }

  # 4. Save objects, tables and plots, grouped into one folder per dataset
  summaries <- list()
  for (ds in unique(vapply(results, function(r) r$dataset, ""))) {
    rs <- Filter(function(r) r$dataset == ds, results)
    summaries[[ds]] <- write_dataset_outputs(rs, file.path(out_dir, ds))
  }

  # 5. Run-level summary and package versions
  write_overview(summaries, failures, samples, out_dir)
  writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
  log_msg("Done: %d sample(s) OK, %d failed. Outputs in %s", length(results), length(failures), out_dir)

  invisible(list(results = results, failures = failures, out_dir = out_dir))
}

## =============================================================================
## Run
## =============================================================================
args <- commandArgs(trailingOnly = TRUE)
run <- main(sample_sheet_path = if (length(args) >= 1) args[1] else SAMPLE_SHEET_PATH,
            out_dir           = if (length(args) >= 2) args[2] else OUT_DIR)
if (length(run$failures) && !interactive()) quit(status = 1)   # non-zero exit for scripts/SLURM
