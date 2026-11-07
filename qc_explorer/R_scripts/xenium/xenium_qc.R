#!/usr/bin/env Rscript
###############################################################################
# xenium_qc.R — sample-sheet-driven quality control for 10x Xenium (human data)
#
# Written to match scRNA/scrna_qc.R and visium/visium_qc.R: same sample-sheet
# rules, same MAD-based cutoffs with hard floors/caps, same output layout.
# The QC app (app_scripts/) uses the functions in this file; it can also run on its own:
#
# Usage
#   Rscript R_scripts/xenium/xenium_qc.R <samplesheet.csv> [out_dir]
#   ...or open bladder cancer.Rproj in RStudio, set SAMPLE_SHEET_PATH below, and click Source.
#
# The sample sheet has ONE ROW PER SAMPLE (one Xenium region / output folder).
# Required columns: sample_id, path. Every other column is optional; a blank
# cell means "use the default" in DEFAULTS below. Each sample is QC'd and saved
# as its OWN Seurat object (with cell centroids as an FOV); nothing is merged.
#
# Per sample, the script:
#   1. reads the cell x gene matrix (cell_feature_matrix.h5 or the
#      cell_feature_matrix/ folder) and the per-cell table (cells.csv.gz or
#      cells.parquet), or a Seurat .rds made with LoadXenium()
#   2. computes QC metrics per cell: transcripts, genes, cell and nucleus area,
#      and the % of decoded counts that came from negative-control probes /
#      codewords (a false-positive rate for that cell)
#   3. removes outlier cells with per-sample MAD thresholds (log scale), hard
#      floors/caps, a cell-area window and a % control cap
#   4. drops rare genes (off by default: Xenium panels are targeted) and saves the
#      sample (raw counts, not normalized), with QC plots and spatial maps
#
# Output folder (the same layout as scrna_qc.R; see there): samplesheet.csv (every setting
# used, reusable as the input sheet), summary.csv, log.txt, and per sample
# <dataset>/<sample_id>/<sample_id>.rds, <sample_id>_cells.csv, qc/ and annotation/.
# Extra sheet columns: patient, cutoff_<name> (see CUTOFFS), annot_<name>; others are kept.
#
# Not done here, on purpose: doublet detection and ambient-RNA correction
# (droplet-based ideas that don't apply to imaging-based segmentation).
###############################################################################

## =============================================================================
## Settings (used when no command-line arguments are given)
## =============================================================================
SAMPLE_SHEET_PATH <- "datasets/XENIUM_DATASET/samplesheet_xenium.csv"   # <- change to your dataset
OUT_DIR           <- NULL   # NULL = "<folder of the sample sheet>/qc_output" (or, for an output folder's own
                            #        samplesheet.csv, that output folder: the rerun replaces its results)

## =============================================================================
## Defaults for every sample-sheet column
## =============================================================================
DEFAULTS <- list(
  # what the sample is
  dataset         = "dataset1", # output folder name; samples are never merged
  format          = "auto",     # auto | xenium | rds
  # cell filters (per-sample MAD cutoffs on the log scale, with hard floors/caps)
  nmads_low       = 3,          # MADs below the median (log scale) for transcripts/genes
  nmads_high      = 5,          # MADs above the median (log scale) for transcripts/genes
  nmads_area      = 5,          # MADs below/above the median (log scale) for cell area
  min_counts      = 10,         # hard floor on transcripts per cell (always applied)
  min_features    = 5,          # hard floor on genes per cell (always applied)
  max_counts      = NA,         # optional hard cap on transcripts
  max_features    = NA,         # optional hard cap on genes
  min_cell_area   = NA,         # optional hard floor on cell area (um^2)
  max_cell_area   = NA,         # optional hard cap on cell area (um^2)
  max_pct_control = 5,          # max % of a cell's decoded counts from negative controls
  # genes and reproducibility
  min_cells_gene  = 0,          # keep genes seen in >= this many kept cells (0 = keep the whole panel)
  seed            = 1234
)
# Cutoffs a sample-sheet column can fix for one sample: cutoff_<name>, e.g. cutoff_min_cell_area = 20
# (Inf = no maximum). Blank = the rules above. The QC app fills these in when you set a cutoff
# by hand, so its samplesheet.csv gives the same cutoffs here.
CUTOFFS <- c("min_counts", "max_counts", "min_features", "max_features", "min_cell_area", "max_cell_area", "max_pct_control")
LOGICAL_PARAMS   <- character(0)
CHARACTER_PARAMS <- c("dataset", "format")
NUMERIC_PARAMS   <- setdiff(names(DEFAULTS), c(LOGICAL_PARAMS, CHARACTER_PARAMS))
KNOWN_FORMATS    <- c("xenium", "rds")

# Columns of cells.csv.gz that count negative-control signal (not every Xenium version has all)
CONTROL_COLUMNS <- c("control_probe_counts", "control_codeword_counts", "unassigned_codeword_counts",
                     "deprecated_codeword_counts", "genomic_control_counts")

REASON_COLOURS <- c(pass = "grey75", low_counts = "#E69F00", low_features = "#D55E00",
                    high_counts = "#56B4E9", high_features = "#0072B2", small_area = "#8C564B",
                    large_area = "#17BECF", high_control = "#CC79A7")

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
  cat(sprintf("\n==== %s  xenium_qc.R %s ====\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), what),
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

resolve_path <- function(p, sheet_dir) {
  if (is.na(p)) return(NA_character_)
  p <- gsub("\\\\", "/", p)
  if (!is_absolute(p)) p <- file.path(sheet_dir, p)
  normalizePath(p, winslash = "/", mustWork = FALSE)
}

guess_format <- function(p) {
  if (dir.exists(p)) return("xenium")
  if (grepl("\\.rds$", p, ignore.case = TRUE)) return("rds")
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
  p <- c(p, sheet_extras(row, sid, add))

  if (is.na(p$path) || !file.exists(p$path)) {
    add("[%s] path not found: %s", sid, row$path)
  } else if (p$format == "auto") {
    p$format <- guess_format(p$path)
    if (is.na(p$format)) add("[%s] can't guess the format of %s; fill in the format column", sid, p$path)
  }
  if (!is.na(p$format) && p$format != "auto" && !p$format %in% KNOWN_FORMATS)
    add("[%s] unknown format '%s' (use %s)", sid, p$format, paste(KNOWN_FORMATS, collapse = ", "))
  for (nm in c("nmads_low", "nmads_high", "nmads_area"))
    if (isTRUE(p[[nm]] <= 0)) add("[%s] %s must be > 0", sid, nm)
  if (isTRUE(p$min_cell_area > p$max_cell_area))
    add("[%s] min_cell_area is larger than max_cell_area", sid)
  list(params = p, problems = problems)
}

read_sample_sheet <- function(sample_sheet_path, sheet_dir) {
  sheet <- read.csv(sample_sheet_path, colClasses = "character", na.strings = c("", "NA"),
                    check.names = FALSE, strip.white = TRUE)
  names(sheet) <- tolower(trimws(names(sheet)))
  sheet <- sheet[rowSums(!is.na(sheet)) > 0, , drop = FALSE]

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

check_packages <- function(samples) {
  need <- c("hdf5r")
  if (any(vapply(samples, function(p) p$format == "xenium", TRUE))) need <- c(need, "data.table")
  missing_pkgs <- need[!vapply(need, requireNamespace, TRUE, quietly = TRUE)]
  if (length(missing_pkgs)) {
    stop("Missing packages: ", paste(missing_pkgs, collapse = ", "),
         "\nInstall with: install.packages(c(", paste0('"', missing_pkgs, '"', collapse = ", "), "))")
  }
}

## =============================================================================
## Readers: Xenium output folder (any layout, incl. GEO-prefixed files) or a Seurat .rds
## =============================================================================
find_file <- function(dir, pattern, what, required = TRUE) {
  f <- list.files(dir, pattern = pattern, full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  if (length(f) > 1) stop("More than one ", what, " file under ", dir, ": ", paste(basename(f), collapse = ", "),
                          ". Put each sample in its own folder.")
  if (length(f) == 0) {
    if (required) stop("No ", what, " file found under ", dir)
    return(NA_character_)
  }
  f
}

# Returns list(genes = genes x cells counts, controls = per-cell control counts from the matrix, or NULL)
read_xenium_matrix <- function(dir) {
  h5 <- find_file(dir, "cell_feature_matrix\\.h5$", "cell_feature_matrix.h5", required = FALSE)
  if (!is.na(h5)) {
    m <- Read10X_h5(h5, use.names = TRUE)
  } else {
    mtx  <- find_file(dir, "matrix\\.mtx(\\.gz)?$", "matrix.mtx (cell_feature_matrix folder)")
    bcs  <- find_file(dirname(mtx), "barcodes\\.tsv(\\.gz)?$", "barcodes.tsv")
    feat <- find_file(dirname(mtx), "features\\.tsv(\\.gz)?$", "features.tsv")
    ft <- read.delim(feat, header = FALSE, stringsAsFactors = FALSE)
    all <- ReadMtx(mtx = mtx, cells = bcs, features = feat, feature.column = 2, unique.features = TRUE)
    m <- if (ncol(ft) >= 3) {
      types <- ft[[3]]
      lapply(split(seq_along(types), types), function(i) all[i, , drop = FALSE])
    } else all
  }
  if (!is.list(m)) return(list(genes = m, controls = NULL))
  if (!"Gene Expression" %in% names(m)) stop("No 'Gene Expression' features in the Xenium matrix")
  ctrl <- m[setdiff(names(m), "Gene Expression")]
  controls <- if (length(ctrl)) Reduce(`+`, lapply(ctrl, Matrix::colSums)) else NULL
  list(genes = m[["Gene Expression"]], controls = controls)
}

read_cells_table <- function(dir) {
  f <- find_file(dir, "cells\\.csv(\\.gz)?$", "cells.csv.gz", required = FALSE)
  if (!is.na(f)) return(as.data.frame(data.table::fread(f, data.table = FALSE)))
  f <- find_file(dir, "cells\\.parquet$", "cells.csv.gz or cells.parquet")
  if (!requireNamespace("arrow", quietly = TRUE))
    stop("only cells.parquet was found, which needs the arrow package: install.packages(\"arrow\")")
  as.data.frame(arrow::read_parquet(f))
}

# Build list(counts, cells) for a Xenium output folder; cells has one row per matrix column
load_xenium_dir <- function(dir, p) {
  mx <- read_xenium_matrix(dir)
  counts <- mx$genes
  rownames(counts) <- make.unique(gsub("_", "-", rownames(counts)))   # Seurat forbids "_" in gene names
  cells <- read_cells_table(dir)
  if (!"cell_id" %in% names(cells)) stop("cells table has no cell_id column")
  cells$cell_id <- as.character(cells$cell_id)
  idx <- match(colnames(counts), cells$cell_id)
  if (all(is.na(idx))) stop("No cell IDs shared between the matrix and the cells table")
  if (anyNA(idx)) {
    log_msg("[%s] WARNING: %d matrix cells are missing from the cells table; dropping them", p$sample_id, sum(is.na(idx)))
    counts <- counts[, !is.na(idx), drop = FALSE]
    idx <- idx[!is.na(idx)]
  }
  cells <- cells[idx, , drop = FALSE]
  ctrl_cols <- intersect(CONTROL_COLUMNS, names(cells))
  control <- if (length(ctrl_cols)) {
    rowSums(as.matrix(cells[, ctrl_cols, drop = FALSE]))
  } else if (!is.null(mx$controls)) {
    mx$controls[colnames(counts)]
  } else {
    log_msg("[%s] WARNING: no negative-control counts found; %% control will be 0", p$sample_id)
    rep(0, ncol(counts))
  }
  info <- data.frame(
    x = cells$x_centroid, y = cells$y_centroid, control_counts = as.numeric(control),
    cell_area    = if ("cell_area"    %in% names(cells)) cells$cell_area    else NA_real_,
    nucleus_area = if ("nucleus_area" %in% names(cells)) cells$nucleus_area else NA_real_,
    row.names = colnames(counts))
  list(counts = counts, info = info, obj = NULL)
}

# A Seurat object made with LoadXenium(): counts from the Xenium assay, controls from nCount_* columns
load_xenium_rds <- function(path, p) {
  obj <- readRDS(path)
  if (!inherits(obj, "Seurat")) stop("The .rds holds a ", class(obj)[1], "; expected a Seurat object")
  assay <- if ("Xenium" %in% names(obj@assays)) "Xenium" else DefaultAssay(obj)
  a <- obj[[assay]]
  if (inherits(a, "Assay5") && length(SeuratObject::Layers(a, search = "counts")) > 1) a <- SeuratObject::JoinLayers(a)
  counts <- SeuratObject::LayerData(a, layer = "counts")
  md <- obj[[]]
  ctrl_cols <- grep("^nCount_(ControlProbe|ControlCodeword|BlankCodeword|UnassignedCodeword|GenomicControl|DeprecatedCodeword)",
                    names(md), value = TRUE)
  tc <- tryCatch(GetTissueCoordinates(obj), error = function(e) NULL)
  xy <- if (!is.null(tc) && "cell" %in% names(tc)) tc[match(colnames(counts), tc$cell), c("x", "y")] else
    data.frame(x = rep(NA_real_, ncol(counts)), y = NA_real_)
  info <- data.frame(
    x = xy$x, y = xy$y,
    control_counts = if (length(ctrl_cols)) rowSums(as.matrix(md[colnames(counts), ctrl_cols, drop = FALSE])) else 0,
    cell_area    = if ("cell_area"    %in% names(md)) md[colnames(counts), "cell_area"]    else NA_real_,
    nucleus_area = if ("nucleus_area" %in% names(md)) md[colnames(counts), "nucleus_area"] else NA_real_,
    row.names = colnames(counts))
  list(counts = counts, info = info, obj = obj)
}

## =============================================================================
## QC steps
## =============================================================================
mad_bounds <- function(x, nlow, nhigh) {
  med <- median(x); dev <- mad(x)
  c(low = med - nlow * dev, high = med + nhigh * dev)
}
na_inf <- function(x) if (is.na(x)) Inf else x
na_0   <- function(x) if (is.na(x)) 0 else x

# Step 1: read the sample, cell names prefixed with the sample ID
load_sample <- function(p) {
  log_msg("[%s] reading %s (%s)", p$sample_id, p$path, p$format)
  x <- if (p$format == "rds") load_xenium_rds(p$path, p) else load_xenium_dir(p$path, p)
  new_names <- paste0(p$sample_id, "_", colnames(x$counts))
  if (!is.null(x$obj)) x$obj <- RenameCells(x$obj, new.names = new_names)
  colnames(x$counts) <- new_names
  rownames(x$info) <- new_names
  log_msg("[%s] %d genes x %d cells", p$sample_id, nrow(x$counts), ncol(x$counts))
  x
}

# Step 2: per-cell QC metrics (one row per cell)
compute_qc_metrics <- function(x, p) {
  m <- x$counts
  n_counts <- Matrix::colSums(m)
  meta <- data.frame(
    cell = colnames(m), sample = p$sample_id, dataset = p$dataset,
    nCount_Xenium = n_counts, nFeature_Xenium = Matrix::colSums(m > 0),
    control_counts = x$info$control_counts,
    pct_control = 100 * x$info$control_counts / pmax(n_counts + x$info$control_counts, 1),
    cell_area = x$info$cell_area, nucleus_area = x$info$nucleus_area,
    nucleus_ratio = x$info$nucleus_area / x$info$cell_area,
    x = x$info$x, y = x$info$y,
    row.names = colnames(m), stringsAsFactors = FALSE)
  meta
}

# Step 3a: this sample's cutoffs. Transcripts and genes: same rule as scrna_qc.R / visium_qc.R.
# Cell area: log scale, median -/+ nmads_area MADs, within the optional hard floor/cap.
compute_thresholds <- function(meta, p) {
  pos <- function(v) log10(pmax(v, 1))
  lc <- mad_bounds(pos(meta$nCount_Xenium),   p$nmads_low, p$nmads_high)
  lf <- mad_bounds(pos(meta$nFeature_Xenium), p$nmads_low, p$nmads_high)
  has_area <- any(is.finite(meta$cell_area) & meta$cell_area > 0)
  la <- if (has_area) mad_bounds(log10(meta$cell_area[meta$cell_area > 0]), p$nmads_area, p$nmads_area) else c(low = -Inf, high = Inf)
  th <- data.frame(sample = p$sample_id, dataset = p$dataset, cells_in = nrow(meta),
             min_counts    = max(p$min_counts,   floor(10^lc[["low"]])),
             max_counts    = min(na_inf(p$max_counts),   ceiling(10^lc[["high"]])),
             min_features  = max(p$min_features, floor(10^lf[["low"]])),
             max_features  = min(na_inf(p$max_features), ceiling(10^lf[["high"]])),
             min_cell_area = round(max(na_0(p$min_cell_area), 10^la[["low"]]), 2),
             max_cell_area = round(min(na_inf(p$max_cell_area), 10^la[["high"]]), 2),
             max_pct_control = na_inf(p$max_pct_control),
             median_counts = median(meta$nCount_Xenium), median_features = median(meta$nFeature_Xenium),
             median_cell_area = if (has_area) round(median(meta$cell_area, na.rm = TRUE), 1) else NA_real_,
             median_pct_control = round(median(meta$pct_control), 3))
  fix_cutoffs(th, p)
}

# Cutoffs fixed in the sample sheet (cutoff_<name> columns) replace the MAD ones
fix_cutoffs <- function(th, p) {
  for (nm in fixed_cutoffs(p)) th[[nm]] <- p[[paste0("cutoff_", nm)]]
  th
}

# Step 3b: flag cells outside the cutoffs; adds qc_reason and keep columns
apply_cell_filters <- function(meta, th, p) {
  area <- ifelse(is.na(meta$cell_area), NA, meta$cell_area)
  checks <- list(
    low_counts    = meta$nCount_Xenium   <  th$min_counts,
    high_counts   = meta$nCount_Xenium   >  th$max_counts,
    low_features  = meta$nFeature_Xenium <  th$min_features,
    high_features = meta$nFeature_Xenium >  th$max_features,
    small_area    = !is.na(area) & area  <  th$min_cell_area,
    large_area    = !is.na(area) & area  >  th$max_cell_area,
    high_control  = meta$pct_control     >  th$max_pct_control)
  meta$qc_reason <- "pass"
  for (nm in rev(names(checks))) meta$qc_reason[checks[[nm]]] <- nm   # first failing check wins
  meta$keep <- meta$qc_reason == "pass"
  log_msg("[%s] cell filters: %d of %d cells pass (min transcripts %d, min genes %d, area %.1f-%.1f um^2, max %% control %.1f)",
          p$sample_id, sum(meta$keep), nrow(meta), as.integer(th$min_counts), as.integer(th$min_features),
          th$min_cell_area, th$max_cell_area, th$max_pct_control)
  meta
}

qc_one_sample <- function(p) {
  set.seed(p$seed)
  x    <- load_sample(p)                       # 1. read
  meta <- compute_qc_metrics(x, p)             # 2. metrics
  th   <- compute_thresholds(meta, p)          # 3a. cutoffs
  meta <- apply_cell_filters(meta, th, p)      # 3b. filter
  list(sample_id = p$sample_id, dataset = p$dataset, data = x, meta = meta, thresholds = th,
       min_cells_gene = p$min_cells_gene, sample_meta = list(patient = p$patient))
}

## =============================================================================
## Outputs
## =============================================================================
META_COLS <- c("sample", "dataset", "control_counts", "pct_control", "cell_area", "nucleus_area", "nucleus_ratio", "x", "y")

# Step 4: drop rare genes (if asked) and save the kept cells with their centroids
save_sample <- function(r, sdir) {
  dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
  kept <- r$meta$cell[r$meta$keep]
  m <- r$data$counts[, kept, drop = FALSE]
  genes_keep <- rownames(m)[Matrix::rowSums(m > 0) >= r$min_cells_gene]
  if (!is.null(r$data$obj)) {
    seu <- subset(r$data$obj, cells = kept, features = genes_keep)
    seu <- AddMetaData(seu, r$meta[kept, setdiff(META_COLS, c("x", "y"))])
  } else {
    seu <- CreateSeuratObject(counts = m[genes_keep, , drop = FALSE], assay = "Xenium", project = r$sample_id,
                              meta.data = r$meta[kept, META_COLS], min.cells = 0, min.features = 0)
    xy <- r$meta[kept, c("x", "y")]
    if (all(is.finite(xy$x)) && all(is.finite(xy$y))) {
      cents <- SeuratObject::CreateCentroids(data.frame(x = xy$x, y = xy$y, cell = kept))
      seu[["fov"]] <- SeuratObject::CreateFOV(cents, type = "centroids", assay = "Xenium")
    }
  }
  seu$orig.ident <- r$sample_id
  Idents(seu) <- "sample"
  for (nm in names(r$sample_meta)) if (!is.na(r$sample_meta[[nm]])) seu@meta.data[[nm]] <- r$sample_meta[[nm]]   # patient, ...
  f <- file.path(sdir, paste0(r$sample_id, ".rds"))
  saveRDS(seu, f)
  log_msg("[%s] saved %s: %d cells, %d genes (in >= %d cells)", r$sample_id, f, ncol(seu), nrow(seu), r$min_cells_gene)
  invisible(f)
}

summarise_sample <- function(r) {
  m <- r$meta
  n <- function(reason) sum(m$qc_reason == reason)
  data.frame(dataset = r$dataset, sample = r$sample_id, status = "ok", cells_in = nrow(m),
             low_counts = n("low_counts"), high_counts = n("high_counts"),
             low_features = n("low_features"), high_features = n("high_features"),
             small_area = n("small_area"), large_area = n("large_area"), high_control = n("high_control"),
             cells_kept = sum(m$keep), pct_kept = round(100 * mean(m$keep), 1))
}

MAX_SPATIAL_POINTS <- 3e5   # spatial maps draw at most this many cells per sample (random subset)

spatial_panels <- function(meta, fill_fun, title) {
  plots <- lapply(split(meta, meta$sample), function(d) {
    d <- d[is.finite(d$x) & is.finite(d$y), ]
    if (nrow(d) > MAX_SPATIAL_POINTS) { set.seed(1); d <- d[sample(nrow(d), MAX_SPATIAL_POINTS), ] }
    fill_fun(ggplot(d, aes(x, y))) + scale_y_reverse() + coord_fixed() +
      labs(title = d$sample[1], x = NULL, y = NULL) +
      theme_void(base_size = 10) + theme(legend.position = "right", plot.title = element_text(hjust = 0.5))
  })
  wrap_plots(plots, ncol = min(length(plots), 3)) + plot_annotation(title = title)
}

make_plots <- function(meta, thresholds, plot_dir, extra_dir = NULL) {   # extra_dir: unused (no doublet / ambient steps)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  theme_set(theme_bw(base_size = 11))
  ns <- length(unique(meta$sample)); nc <- min(ns, 3); nr <- ceiling(ns / nc)
  w <- max(6, 3.2 * nc)

  # 01: violins before vs after, cutoffs dashed
  vdf <- rbind(transform(meta, stage = "before"), transform(meta[meta$keep, ], stage = "after"))
  vdf$stage <- factor(vdf$stage, levels = c("before", "after"))
  lines <- function(cols) do.call(rbind, lapply(cols, function(cl) data.frame(sample = thresholds$sample, y = thresholds[[cl]])))
  vln <- function(var, cols, log_y, title) {
    ld <- lines(cols); ld <- ld[is.finite(ld$y) & (!log_y | ld$y > 0), ]
    d <- vdf[is.finite(vdf[[var]]) & (!log_y | vdf[[var]] > 0), ]
    g <- ggplot(d, aes(stage, .data[[var]], fill = stage)) +
      geom_violin(scale = "width", linewidth = 0.3) +
      facet_wrap(~ sample, ncol = nc) +
      scale_fill_manual(values = c(before = "grey75", after = "#4C9BE8"), guide = "none") +
      labs(title = title, x = NULL, y = NULL)
    if (nrow(ld)) g <- g + geom_hline(data = ld, aes(yintercept = y), linetype = "dashed", colour = "red3")
    if (log_y) g <- g + scale_y_log10()
    g
  }
  panels <- list(vln("nCount_Xenium", c("min_counts", "max_counts"), TRUE, "Transcripts per cell (log)"),
                 vln("nFeature_Xenium", c("min_features", "max_features"), TRUE, "Genes per cell (log)"))
  if (any(is.finite(meta$cell_area)))
    panels <- c(panels, list(vln("cell_area", c("min_cell_area", "max_cell_area"), TRUE, "Cell area, um^2 (log)")))
  panels <- c(panels, list(vln("pct_control", "max_pct_control", FALSE, "% counts from negative controls")))
  p1 <- wrap_plots(panels, ncol = 1)
  ggsave(file.path(plot_dir, "01_violin_before_after.png"), p1, width = w, height = 3 * nr * length(panels), dpi = 150, limitsize = FALSE)

  # 02: transcripts vs genes, coloured by QC result
  meta$qc_reason <- factor(meta$qc_reason, levels = names(REASON_COLOURS))
  p2 <- ggplot(meta, aes(pmax(nCount_Xenium, 1), pmax(nFeature_Xenium, 1), colour = qc_reason)) +
    geom_point(size = 0.3, alpha = 0.6) + scale_x_log10() + scale_y_log10() +
    facet_wrap(~ sample, ncol = nc) + scale_colour_manual(values = REASON_COLOURS, drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    labs(title = "Transcripts vs genes per cell, coloured by QC result", x = "Transcripts (log)", y = "Genes (log)", colour = "QC result")
  ggsave(file.path(plot_dir, "02_counts_vs_genes_by_reason.png"), p2, width = w + 2, height = 3.5 * nr + 0.5, dpi = 150, limitsize = FALSE)

  if (any(is.finite(meta$x))) {
    # 03: WHERE the removed cells are
    p3 <- spatial_panels(meta, function(g) g +
                           geom_point(data = function(d) d[d$keep, ], colour = "grey80", size = 0.05, shape = 16) +
                           geom_point(data = function(d) d[!d$keep, ], aes(colour = qc_reason), size = 0.15, shape = 16) +
                           scale_colour_manual(values = REASON_COLOURS, drop = FALSE, name = "QC result") +
                           guides(colour = guide_legend(override.aes = list(size = 3))),
                         "Removed cells on the tissue (grey = kept)")
    ggsave(file.path(plot_dir, "03_spatial_removed_cells.png"), p3, width = 4.5 * nc + 1.5, height = 4.2 * nr + 0.6, dpi = 150, limitsize = FALSE)
    # 04: transcripts across the tissue
    p4 <- spatial_panels(meta, function(g) g + geom_point(aes(colour = log10(pmax(nCount_Xenium, 1))), size = 0.05, shape = 16) +
                           scale_colour_viridis_c(name = "log10 transcripts"), "Transcripts per cell across the tissue")
    ggsave(file.path(plot_dir, "04_spatial_transcripts.png"), p4, width = 4.5 * nc + 1.5, height = 4.2 * nr + 0.6, dpi = 150, limitsize = FALSE)
  }
}

SCRIPT_NAME   <- "xenium_qc.R"
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
if (length(run$failures) && !interactive()) quit(status = 1)
