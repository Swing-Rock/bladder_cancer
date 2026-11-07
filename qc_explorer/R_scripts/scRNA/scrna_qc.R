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
# nothing is merged. `dataset` only decides which output folder a sample goes in.
# Optional extra columns: patient (added to every cell as `patient`), cutoff_<name>
# (fixes one cutoff for that sample instead of the MAD rule; see CUTOFFS), annot_<name>
# (annotation settings, read by annotate_cells.R); any other column is kept as it is.
#
# Per sample, the pipeline:
#   1. reads the counts (10x folder, 10x .h5, dense text table, or .rds)
#   2. detects whether the matrix still contains empty droplets and, if so, calls cells with emptyDrops
#   3. computes QC metrics (UMIs, genes, % mito; % ribo and % hemoglobin are reported only)
#   4. removes outliers using per-sample MAD thresholds
#   5. removes doublets with scDblFinder
#   6. estimates ambient RNA with DecontX (uses raw droplets as background if given)
#   7. drops rare genes and saves the sample
#
# Output folder (the same layout as the other QC scripts, annotate_cells.R and the QC app)
#   samplesheet.csv   one row per sample with every setting used: give it back to this script
#                     to rerun with exactly the same settings (it then writes to the same folder)
#   summary.csv       one row per sample: cells in, removed for each reason, kept, cutoffs used
#   log.txt           every run adds to it
#   <dataset>/<sample_id>/
#     <sample_id>.rds        the QC'd Seurat object (raw counts; + decontX assay); annotation adds its labels
#     <sample_id>_cells.csv  one row per cell, kept or not: QC metrics, why it was removed, doublet / ambient scores
#     qc/                    QC plots
#     doublets_ambient/      scDblFinder and DecontX plots
#     annotation/            written by annotate_cells.R
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
# Relative to the repo root (where "bladder cancer.Rproj" is), so it works on anyone's copy.
# Open the .Rproj in RStudio (or run Rscript from the repo root) and this path resolves.
SAMPLE_SHEET_PATH <- "datasets/GSE130001/samplesheet_GSE130001.csv"
OUT_DIR           <- NULL   # NULL = "<folder of the sample sheet>/qc_output" (or, for an output folder's own
                            #        samplesheet.csv, that output folder: the rerun replaces its results)

## =============================================================================
## Defaults for every sample-sheet column
## =============================================================================
DEFAULTS <- list(
  # what the sample is
  dataset           = "dataset1", # output folder name; samples are never merged
  format            = "auto",     # auto | 10x | 10x_h5 | dense | dense_t | rds
  gene_column       = 2,          # 10x features file: 1 = Ensembl ID, 2 = gene symbol
  mito_pattern      = "^MT-",     # regex for mitochondrial genes (human gene symbols)
  is_filtered       = "auto",     # auto = detect from the matrix; TRUE = cells only; FALSE = raw droplets (emptyDrops)
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
  # cell calling (only for raw matrices); emptydrops_lower is also the "empty" UMI level used by auto-detection
  emptydrops_lower  = 100,        # droplets with <= this many UMIs define "empty"
  emptydrops_fdr    = 0.001,
  # genes and reproducibility
  min_cells_gene    = 3,          # keep genes seen in >= this many of the sample's kept cells
  seed              = 1234
)
# Cutoffs a sample-sheet column can fix for one sample: cutoff_<name>, e.g. cutoff_max_mito = 15
# (Inf = no maximum). Blank = the MAD rule above. The QC app fills these in when you set a cutoff
# by hand, so its samplesheet.csv gives the same cutoffs here.
CUTOFFS <- c("min_counts", "max_counts", "min_features", "max_features", "max_mito")
LOGICAL_PARAMS   <- c("run_doublets", "run_decontx")
RAW_FRACTION     <- 0.5                         # matrix counts as raw if more than this share of barcodes are empty-looking
CHARACTER_PARAMS <- c("dataset", "format", "mito_pattern", "raw_path", "is_filtered")
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
## Logging: every message goes to the console and to <out_dir>/log.txt (added to, never wiped)
## =============================================================================
.qc_log <- new.env()
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

start_log <- function(out_dir, what = "") {
  .qc_log$file <- file.path(normalizePath(out_dir, winslash = "/"), "log.txt")
  cat(sprintf("\n==== %s  scrna_qc.R %s ====\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), what),
      file = .qc_log$file, append = TRUE)
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

parse_logical <- function(x) {
  v <- toupper(trimws(x))
  if (v %in% c("TRUE", "T", "YES", "Y", "1")) return(TRUE)
  if (v %in% c("FALSE", "F", "NO", "N", "0")) return(FALSE)
  NA
}

guess_format <- function(p) {
  if (dir.exists(p)) return("10x")
  if (grepl("\\.h5$", p, ignore.case = TRUE)) return("10x_h5")
  if (grepl("\\.rds$", p, ignore.case = TRUE)) return("rds")
  if (grepl("\\.(csv|tsv|txt)(\\.gz)?$", p, ignore.case = TRUE)) return("dense")
  NA_character_
}

# The sheet columns that aren't QC parameters: fixed cutoffs (cutoff_<name>), the patient,
# and any other column (kept as is, e.g. notes or annot_* settings for annotate_cells.R)
sheet_extras <- function(row, sid, add) {
  out <- list()
  for (nm in CUTOFFS) {
    col <- paste0("cutoff_", nm)
    if (!col %in% names(row) || is.na(row[[col]])) next
    v <- suppressWarnings(as.numeric(row[[col]]))
    if (is.na(v)) add("[%s] %s must be a number (Inf = no maximum), got '%s'", sid, col, row[[col]]) else out[[col]] <- v
  }
  pc <- intersect(c("patient", "patient_id"), names(row))
  out$patient <- if (length(pc) && !is.na(row[[pc[1]]])) as.character(row[[pc[1]]]) else NA_character_
  other <- setdiff(names(row), c("sample_id", "path", "patient", "patient_id", names(DEFAULTS), paste0("cutoff_", CUTOFFS)))
  out$extra_columns <- lapply(setNames(other, other), function(nm) as.character(row[[nm]]))
  out
}
fixed_cutoffs <- function(p) CUTOFFS[vapply(CUTOFFS, function(nm) { v <- p[[paste0("cutoff_", nm)]]; length(v) == 1 && !is.na(v) }, TRUE)]

# Turn one sheet row into a full parameter list. Returns list(params, problems).
resolve_sample <- function(row, i, sheet_dir) {
  problems <- character(0)
  add <- function(...) problems <<- c(problems, sprintf(...))
  sid <- row$sample_id
  p <- DEFAULTS
  for (nm in names(DEFAULTS)) {
    if (!nm %in% names(row) || is.na(row[[nm]])) next
    v <- row[[nm]]
    if (nm == "is_filtered") {
      v <- if (tolower(trimws(v)) == "auto") "auto" else as.character(parse_logical(v))
      if (is.na(v)) add("[%s] is_filtered must be auto, TRUE or FALSE, got '%s'", sid, row[[nm]])
    } else if (nm %in% LOGICAL_PARAMS) {
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
  p <- c(p, sheet_extras(row, sid, add))

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

  other_cols <- setdiff(names(sheet), c("sample_id", "path", "patient", "patient_id", names(DEFAULTS), paste0("cutoff_", CUTOFFS)))
  if (length(other_cols)) log_msg("Not QC settings, copied as they are to samplesheet.csv: %s", paste(other_cols, collapse = ", "))

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
  if (any(vapply(samples, function(p) p$is_filtered == "FALSE", TRUE))) need <- c(need, "DropletUtils")
  if (any(vapply(samples, function(p) p$format == "10x_h5", TRUE))) need <- c(need, "hdf5r")
  missing_pkgs <- need[!vapply(need, requireNamespace, TRUE, quietly = TRUE)]
  if (length(missing_pkgs)) {
    stop("Missing packages: ", paste(missing_pkgs, collapse = ", "),
         "\nInstall with: BiocManager::install(c(", paste0('"', missing_pkgs, '"', collapse = ", "), "))",
         "\nor set run_doublets / run_decontx to FALSE in the sample sheet.")
  }
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

# Is this a raw droplet matrix? Raw matrices are mostly empty barcodes (a few UMIs each);
# filtered matrices have almost none. Returns list(raw = TRUE/FALSE, frac_empty, n_empty).
detect_raw <- function(m, lower) {
  n_empty <- sum(Matrix::colSums(m) <= lower)
  frac <- n_empty / max(ncol(m), 1)
  list(raw = frac > RAW_FRACTION, frac_empty = frac, n_empty = n_empty)
}

# Decide whether to call cells, checking the sheet's is_filtered against what the matrix looks like
resolve_filter_state <- function(m, p) {
  d <- detect_raw(m, p$emptydrops_lower)
  looks <- if (d$raw) "raw (mostly empty droplets)" else "filtered (cells only)"
  evidence <- sprintf("%s of %s barcodes have <= %d UMIs (%.1f%%)", format(d$n_empty, big.mark = ","),
                      format(ncol(m), big.mark = ","), as.integer(p$emptydrops_lower), 100 * d$frac_empty)
  if (p$is_filtered == "auto") {
    log_msg("[%s] matrix looks %s: %s", p$sample_id, looks, evidence)
    return(!d$raw)
  }
  said_filtered <- p$is_filtered == "TRUE"
  if (said_filtered == d$raw) {
    stop(sprintf("is_filtered = %s in the sample sheet, but the matrix looks %s: %s. Check the path (filtered_ vs raw_ matrix), or set is_filtered to auto.",
                 p$is_filtered, looks, evidence))
  }
  log_msg("[%s] is_filtered = %s, consistent with the matrix: %s", p$sample_id, p$is_filtered, evidence)
  said_filtered
}

# Step 2: call cells on unfiltered data, and collect empty droplets as DecontX background.
# Returns list(counts = cells only, background = empty droplets or NULL, n_barcodes, matrix_type).
call_cells <- function(m, p) {
  background <- NULL
  n_barcodes <- ncol(m)
  is_filtered <- resolve_filter_state(m, p)
  if (!is_filtered) {
    if (!requireNamespace("DropletUtils", quietly = TRUE))
      stop("this is a raw matrix, so cells must be called with emptyDrops, but DropletUtils is not installed: BiocManager::install(\"DropletUtils\")")
    if (!is.na(p$raw_path)) log_msg("[%s] raw_path ignored: the main matrix is already raw and supplies its own empty droplets", p$sample_id)
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
  list(counts = m, background = background, n_barcodes = n_barcodes,
       matrix_type = if (is_filtered) "filtered" else "raw (emptyDrops)")
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
  out <- data.frame(sample = p$sample_id, dataset = p$dataset,
                    barcodes_in = n_barcodes, cells_in = nrow(meta), th,
                    median_counts = median(meta$nCount_RNA), median_features = median(meta$nFeature_RNA),
                    median_mito = round(median(meta$percent.mt), 2), n_mito_genes = attr(meta, "n_mito_genes"))
  fix_cutoffs(out, p)
}

# Cutoffs fixed in the sample sheet (cutoff_<name> columns) replace the MAD ones
fix_cutoffs <- function(th, p) {
  for (nm in fixed_cutoffs(p)) th[[nm]] <- p[[paste0("cutoff_", nm)]]
  th
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

# The score above which scDblFinder called cells doublets. scDblFinder stores it in the
# object's metadata; if it's missing, use the lowest score among cells it called doublets.
doublet_cutoff <- function(res) {
  thr <- tryCatch(S4Vectors::metadata(res)$scDblFinder.threshold, error = function(e) NULL)
  if (is.numeric(thr) && length(thr) == 1 && is.finite(thr)) return(thr)
  s <- res$scDblFinder.score[res$scDblFinder.class == "doublet"]
  if (length(s)) min(s) else NA_real_
}

# Step 5: doublets with scDblFinder, on raw counts of cells still kept.
# Returns list(meta, threshold): threshold is the doublet score cutoff (NA if not run).
detect_doublets <- function(meta, m, p) {
  meta$scDblFinder.score <- NA_real_
  meta$scDblFinder.class <- NA_character_
  skipped <- list(meta = meta, threshold = NA_real_)
  if (!p$run_doublets) return(skipped)
  keep <- meta$cell[meta$keep]
  if (length(keep) < 50) {
    log_msg("[%s] WARNING: only %d cells; skipping doublet detection", p$sample_id, length(keep))
    return(skipped)
  }
  res <- tryCatch({
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = m[, keep]))
    set.seed(p$seed)
    if (is.na(p$doublet_rate)) scDblFinder::scDblFinder(sce) else scDblFinder::scDblFinder(sce, dbr = p$doublet_rate)
  }, error = function(e) { log_msg("[%s] WARNING: scDblFinder failed (%s); no doublets removed", p$sample_id, conditionMessage(e)); NULL })
  if (is.null(res)) return(skipped)
  meta[keep, "scDblFinder.score"] <- res$scDblFinder.score
  meta[keep, "scDblFinder.class"] <- as.character(res$scDblFinder.class)
  dbl <- keep[res$scDblFinder.class == "doublet"]
  meta[dbl, "qc_reason"] <- "doublet"
  meta[dbl, "keep"] <- FALSE
  thr <- doublet_cutoff(res)
  log_msg("[%s] scDblFinder: %d doublets among %d cells (%.1f%%; score cutoff %.3f)",
          p$sample_id, length(dbl), length(keep), 100 * length(dbl) / length(keep), thr)
  list(meta = meta, threshold = thr)
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
  th$input_matrix <- called$matrix_type
  meta    <- apply_cell_filters(meta, th, p)                    # 4b. filter
  dbl     <- detect_doublets(meta, m, p)                        # 5. doublets
  meta    <- dbl$meta
  th$doublet_threshold <- round(dbl$threshold, 4)
  scored  <- !is.na(meta$scDblFinder.class)
  th$doublet_pct       <- if (any(scored)) round(100 * mean(meta$qc_reason[scored] == "doublet"), 2) else NA_real_
  amb     <- correct_ambient(meta, m, called$background, p)     # 6. ambient RNA
  meta    <- amb$meta

  kept <- meta$cell[meta$keep]
  list(sample_id = p$sample_id, dataset = p$dataset, meta = meta, thresholds = th,
       counts  = m[, kept, drop = FALSE],
       decontx = if (is.null(amb$decontx)) NULL else amb$decontx[, kept, drop = FALSE],
       min_cells_gene = p$min_cells_gene, sample_meta = list(patient = p$patient))
}

## =============================================================================
## Per-sample outputs: Seurat object, summary numbers, plots
## =============================================================================
# Step 7: drop rare genes and save one Seurat object for this sample
save_sample <- function(r, sdir) {
  dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
  sid <- r$sample_id
  genes_keep <- rownames(r$counts)[Matrix::rowSums(r$counts > 0) >= r$min_cells_gene]
  md <- r$meta[colnames(r$counts), c("sample", "dataset", "percent.mt", "percent.ribo", "percent.hb",
                                      "log10GenesPerUMI", "scDblFinder.score", "decontX_contamination")]
  seu <- CreateSeuratObject(counts = r$counts[genes_keep, , drop = FALSE], project = sid,
                            meta.data = md, min.cells = 0, min.features = 0)
  seu$orig.ident <- sid
  Idents(seu) <- "sample"
  if (!is.null(r$decontx)) seu[["decontX"]] <- CreateAssayObject(counts = r$decontx[genes_keep, colnames(seu), drop = FALSE])
  for (nm in names(r$sample_meta)) if (!is.na(r$sample_meta[[nm]])) seu@meta.data[[nm]] <- r$sample_meta[[nm]]   # patient, ...
  f <- file.path(sdir, paste0(sid, ".rds"))
  saveRDS(seu, f)
  log_msg("[%s] saved %s: %d cells, %d genes (in >= %d cells)", sid, f, ncol(seu), nrow(seu), r$min_cells_gene)
  invisible(f)
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

# QC plots in plot_dir; the doublet and ambient-RNA plots (when those steps ran) in extra_dir
make_plots <- function(meta, thresholds, plot_dir, extra_dir = plot_dir) {
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
    cut_df <- thresholds[is.finite(thresholds$doublet_threshold), c("sample", "doublet_threshold", "doublet_pct")]
    p4 <- ggplot(d, aes(scDblFinder.score, fill = scDblFinder.class)) + geom_histogram(bins = 50) +
      facet_wrap(~ sample, ncol = nc, scales = "free_y") +
      scale_fill_manual(values = c(singlet = "grey70", doublet = "#009E73")) +
      labs(title = "scDblFinder doublet scores",
           subtitle = "Dashed line = this sample's cutoff; cells to the right are called doublets",
           x = "doublet score", y = "cells", fill = NULL)
    if (nrow(cut_df)) {
      p4 <- p4 +
        geom_vline(data = cut_df, aes(xintercept = doublet_threshold), linetype = "dashed", colour = "red3") +
        geom_text(data = cut_df, inherit.aes = FALSE, size = 3, hjust = 1.05, vjust = 1.5, colour = "red3",
                  aes(x = doublet_threshold, y = Inf,
                      label = sprintf("cutoff %.3f\n%.1f%% doublets", doublet_threshold, doublet_pct)))
    }
    dir.create(extra_dir, recursive = TRUE, showWarnings = FALSE)
    ggsave(file.path(extra_dir, "doublet_scores.png"), p4, width = w, height = 3 * nr + 0.5, dpi = 150, limitsize = FALSE)
  }
  d <- meta[!is.na(meta$decontX_contamination), ]
  if (nrow(d)) {
    p5 <- ggplot(d, aes(decontX_contamination)) + geom_histogram(bins = 50, fill = "#4C9BE8") +
      facet_wrap(~ sample, ncol = nc, scales = "free_y") +
      labs(title = "DecontX estimated ambient contamination", x = "fraction of counts from ambient RNA", y = "cells")
    dir.create(extra_dir, recursive = TRUE, showWarnings = FALSE)
    ggsave(file.path(extra_dir, "decontX_contamination.png"), p5, width = w, height = 3 * nr + 0.5, dpi = 150, limitsize = FALSE)
  }
}

SCRIPT_NAME   <- "scrna_qc.R"
UNIT_FILE     <- "cells"     # <sample>_cells.csv
TABLE_COLUMNS <- NULL        # columns of the per-cell table (NULL = all)
SUMMARY_PRINT <- c("dataset", "sample_id", "status", "cells_in", "cells_kept", "pct_kept")

## =============================================================================
## Output folder: samplesheet.csv, summary.csv, log.txt and one folder per sample.
## The same in scrna_qc.R, visium_qc.R and xenium_qc.R; annotate_cells.R adds to it.
## =============================================================================
sample_dir <- function(out_dir, dataset, sample_id) file.path(out_dir, dataset, sample_id)

# Default output folder: <sheet folder>/qc_output, or, when the sheet is an output folder's own
# samplesheet.csv, that folder (the rerun replaces its results)
default_out_dir <- function(sheet) {
  dir <- dirname(normalizePath(sheet, winslash = "/"))
  if (tolower(basename(sheet)) == "samplesheet.csv" && file.exists(file.path(dir, "summary.csv"))) dir else file.path(dir, "qc_output")
}

# A path as written to samplesheet.csv: relative to the output folder when on the same drive, so
# the sheet keeps working when the project folder moves (relative paths are read from the sheet's folder)
relative_path <- function(path, from) {
  if (length(path) != 1 || is.na(path) || !nzchar(path)) return(NA_character_)
  a <- strsplit(normalizePath(path, winslash = "/", mustWork = FALSE), "/")[[1]]
  b <- strsplit(normalizePath(from, winslash = "/", mustWork = FALSE), "/")[[1]]
  same <- function(x, y) if (.Platform$OS.type == "windows") tolower(x) == tolower(y) else x == y
  if (!same(a[1], b[1])) return(paste(a, collapse = "/"))
  n <- 0
  while (n < min(length(a), length(b)) && same(a[n + 1], b[n + 1])) n <- n + 1
  rel <- c(rep("..", length(b) - n), a[-seq_len(n)])
  if (length(rel)) paste(rel, collapse = "/") else "."
}

# Replace one sample's row in a CSV (or add it); the other rows stay as they are
upsert_csv <- function(file, row, key = "sample_id") {
  row <- as.data.frame(lapply(row, function(v) if (is.null(v) || !length(v) || is.na(v[1])) NA_character_ else as.character(v[1])),
                       stringsAsFactors = FALSE, check.names = FALSE)
  out <- row
  if (file.exists(file)) {
    old <- read.csv(file, colClasses = "character", na.strings = "", check.names = FALSE)
    if (key %in% names(old)) {
      cols <- union(names(row), names(old))
      for (cl in setdiff(cols, names(old))) old[[cl]] <- if (nrow(old)) NA_character_ else character(0)
      for (cl in setdiff(cols, names(row))) row[[cl]] <- NA_character_
      out <- old[cols]
      i <- match(row[[key]], out[[key]])
      if (is.na(i)) out <- rbind(out, row[cols]) else out[i, ] <- row[cols]
    }
  }
  write.csv(out, file, row.names = FALSE, na = "")
  invisible(file)
}

# Change some columns of one sample's row (adding the row if it isn't there)
update_csv_row <- function(file, id, values, key = "sample_id") {
  old <- if (file.exists(file)) read.csv(file, colClasses = "character", na.strings = "", check.names = FALSE)
  row <- if (!is.null(old) && id %in% old[[key]]) as.list(old[match(id, old[[key]]), , drop = FALSE]) else setNames(list(id), key)
  for (nm in names(values)) row[[nm]] <- values[[nm]]
  upsert_csv(file, row, key)
}

# One samplesheet.csv row: every setting the sample was run with. Given back to this script,
# the sheet reproduces the run (cutoffs set by hand included, as cutoff_* columns).
sheet_row <- function(p, out_dir) {
  row <- list(sample_id = p$sample_id, path = relative_path(p$path, out_dir), dataset = p$dataset, patient = p$patient %||% NA)
  for (nm in setdiff(names(DEFAULTS), "dataset")) row[[nm]] <- p[[nm]] %||% NA
  if ("raw_path" %in% names(row)) row$raw_path <- relative_path(row$raw_path, out_dir)
  for (nm in CUTOFFS) row[[paste0("cutoff_", nm)]] <- p[[paste0("cutoff_", nm)]] %||% NA
  for (nm in names(p$extra_columns)) row[[nm]] <- p$extra_columns[[nm]]
  row
}

# One summary.csv row: counts per QC result, the cutoffs used, and (filled in by annotate_cells.R) the labels
summary_row <- function(r, p, manual, saved_by) {
  sm <- summarise_sample(r)
  th <- r$thresholds
  row <- c(list(sample_id = r$sample_id, dataset = r$dataset, patient = p$patient %||% NA, status = "ok"),
           as.list(sm[setdiff(names(sm), c("sample", "dataset", "status"))]))
  for (nm in CUTOFFS) row[[paste0("cutoff_", nm)]] <- th[[nm]]
  row$fixed_cutoffs <- if (length(manual)) paste(manual, collapse = ";") else "none"
  for (nm in setdiff(names(th), c("sample", "dataset", CUTOFFS, names(row)))) row[[nm]] <- th[[nm]]
  row[c("annotated_at", "cell_labels", "malignant")] <- NA            # cleared: the cells may have changed
  if (identical(UNIT_FILE, "spots")) row[c("deconvolved_at", "rctd_mean_proportions")] <- NA
  c(row, list(qc_saved_by = saved_by, qc_saved_at = format(Sys.time(), "%Y-%m-%d %H:%M")))
}

# Everything for one sample: <dataset>/<sample>/ (object, per-cell table, plots) and its rows in
# samplesheet.csv and summary.csv. Saving a sample again replaces all of it, an earlier annotation
# included (it was made on the old cells). The QC app calls this too when you save.
# manual: the cutoffs fixed instead of the MAD rule (cutoff_* columns here; set by hand in the app)
write_sample_outputs <- function(r, p, out_dir, manual = fixed_cutoffs(p), saved_by = SCRIPT_NAME) {
  sdir <- sample_dir(out_dir, r$dataset, r$sample_id)
  old <- file.path(sdir, c("qc", "doublets_ambient", "annotation", "deconvolution"))
  unlink(old[dir.exists(old)], recursive = TRUE)
  dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
  rds <- save_sample(r, sdir)
  meta <- if (is.null(TABLE_COLUMNS)) r$meta else r$meta[, intersect(TABLE_COLUMNS, names(r$meta)), drop = FALSE]
  for (nm in names(r$sample_meta)) if (!is.na(r$sample_meta[[nm]])) meta[[nm]] <- r$sample_meta[[nm]]
  table_file <- file.path(sdir, sprintf("%s_%s.csv", r$sample_id, UNIT_FILE))
  write.csv(meta, table_file, row.names = FALSE)
  tryCatch(make_plots(r$meta, r$thresholds, file.path(sdir, "qc"), file.path(sdir, "doublets_ambient")),
           error = function(e) log_msg("[%s] WARNING: plotting failed (%s)", r$sample_id, conditionMessage(e)))
  upsert_csv(file.path(out_dir, "samplesheet.csv"), sheet_row(p, out_dir))
  upsert_csv(file.path(out_dir, "summary.csv"), summary_row(r, p, manual, saved_by))
  log_msg("[%s] saved to %s: %d of %d %s kept%s", r$sample_id, sdir, sum(r$meta$keep), nrow(r$meta), UNIT_FILE,
          if (length(manual)) paste0("; fixed cutoffs (not the MAD rule): ", paste(sprintf("%s = %s", manual, unlist(r$thresholds[manual])), collapse = ", ")) else "")
  invisible(list(dir = sdir, rds = rds, table = table_file))
}

# A sample that failed: its settings still go in samplesheet.csv (fix and rerun), and summary.csv says so
record_failure <- function(p, out_dir, msg) {
  upsert_csv(file.path(out_dir, "samplesheet.csv"), sheet_row(p, out_dir))
  had <- file.exists(file.path(sample_dir(out_dir, p$dataset, p$sample_id), paste0(p$sample_id, ".rds")))
  update_csv_row(file.path(out_dir, "summary.csv"), p$sample_id,
                 list(dataset = p$dataset, status = paste0("failed", if (had) " (the files are from an earlier run)", ": ", msg)))
}

## =============================================================================
## main(): the whole pipeline, in order
## =============================================================================
main <- function(sample_sheet_path = SAMPLE_SHEET_PATH, out_dir = OUT_DIR) {
  # 1. Set up the output folder and log
  if (!file.exists(sample_sheet_path)) stop("Sample sheet not found: ", sample_sheet_path,
                                            "\n(current working directory: ", getwd(), ")")
  sheet_dir <- dirname(normalizePath(sample_sheet_path, winslash = "/"))
  if (is.null(out_dir)) out_dir <- default_out_dir(sample_sheet_path)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir <- normalizePath(out_dir, winslash = "/")
  start_log(out_dir, sample_sheet_path)

  # 2. Read and check the sample sheet (stops here if anything is wrong)
  samples <- read_sample_sheet(sample_sheet_path, sheet_dir)
  check_packages(samples)
  log_msg("Sample sheet OK: %d sample(s) in %d dataset(s); output folder %s", length(samples),
          length(unique(vapply(samples, `[[`, "", "dataset"))), out_dir)

  # 3. QC and save each sample; a failure in one sample doesn't stop the others
  ok <- character(0); failures <- list()
  for (p in samples) {
    r <- NULL
    res <- tryCatch({ r <- qc_one_sample(p); write_sample_outputs(r, p, out_dir); TRUE },
                    error = function(e) { log_msg("[%s] FAILED: %s", p$sample_id, conditionMessage(e)); conditionMessage(e) })
    if (isTRUE(res)) ok <- c(ok, p$sample_id) else { failures[[p$sample_id]] <- res; record_failure(p, out_dir, res) }
    gc(verbose = FALSE)
  }

  # 4. What happened
  sm <- read.csv(file.path(out_dir, "summary.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  print(sm[sm$sample_id %in% names(samples), intersect(SUMMARY_PRINT, names(sm))], row.names = FALSE)
  log_msg("Done: %d sample(s) OK, %d failed. Outputs in %s", length(ok), length(failures), out_dir)
  invisible(list(ok = ok, failures = failures, out_dir = out_dir))
}

## =============================================================================
## Run
## =============================================================================
args <- commandArgs(trailingOnly = TRUE)
run <- main(sample_sheet_path = if (length(args) >= 1) args[1] else SAMPLE_SHEET_PATH,
            out_dir           = if (length(args) >= 2) args[2] else OUT_DIR)
if (length(run$failures) && !interactive()) quit(status = 1)   # non-zero exit for scripts/SLURM
