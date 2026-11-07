#!/usr/bin/env Rscript
###############################################################################
# visium_deconvolution.R — cell-type deconvolution of QC'd Visium spots with RCTD
#
# Usage
#   Rscript visium_deconvolution.R <qc_dir> <reference.rds> <celltype_column>
#     qc_dir           the output folder of visium_qc.R (or the QC app): every
#                      <dataset>/<sample>/<sample>.rds in it is deconvolved
#     reference.rds    an annotated single-cell Seurat object (raw counts) to learn
#                      each cell type's expression profile from
#     celltype_column  the metadata column in the reference that holds cell-type labels
#   ...or open bladder cancer.Rproj in RStudio, set the Settings below, and click Source.
#
# A Visium spot (55 um) captures several cells, so its counts are a mixture.
# RCTD (Cable et al. 2022, Nature Biotechnology; package spacexr, either the
# Bioconductor or the GitHub version: the script detects which) learns a profile
# for each cell type from the reference and estimates, for every spot, what
# proportion of the spot each cell type makes up.
#
# For each sample, the script:
#   1. builds the RCTD reference once from the annotated Seurat object
#   2. runs RCTD in "full" mode (any number of cell types per spot)
#   3. adds RCTD_<celltype> proportion columns + RCTD_dominant to the spot metadata
#   4. writes them into the QC output folder: into <sample>.rds and <sample>_spots.csv, spatial
#      maps into <dataset>/<sample>/deconvolution/, mean proportions into summary.csv, and log.txt
# Running visium_qc.R on a sample again removes its deconvolution: run this again afterwards.
###############################################################################

## =============================================================================
## Settings (used when no command-line arguments are given)
## =============================================================================
# Relative to the repo root (where "bladder cancer.Rproj" is).
QC_DIR          <- "datasets/VISIUM_DATASET/qc_output"                  # <- the visium_qc.R output folder
REFERENCE_RDS   <- "datasets/REFERENCE/reference_annotated.rds"         # <- annotated scRNA-seq Seurat object
CELLTYPE_COLUMN <- "celltype"                                           # <- its cell-type column

PARAMS <- list(
  reference_assay    = "RNA",   # assay in the reference holding raw counts
  min_cells_per_type = 25,      # cell types with fewer cells are dropped (RCTD needs >= 25)
  max_cells_per_type = 10000,   # larger cell types are randomly downsampled to this many cells
  doublet_mode       = "full",  # "full" = any number of types per spot (right for Visium)
  n_cores            = 4,       # parallel cores for RCTD
  min_shared_genes   = 5000,    # warn if reference and spots share fewer genes than this
  seed               = 1234
)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

## =============================================================================
## Logging: console + <qc_dir>/log.txt (added to, never wiped)
## =============================================================================
.dc_log <- new.env()
start_log <- function(qc_dir) {
  .dc_log$file <- file.path(normalizePath(qc_dir, winslash = "/"), "log.txt")
  cat(sprintf("\n==== %s  visium_deconvolution.R ====\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), file = .dc_log$file, append = TRUE)
}
log_msg <- function(...) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(...))
  message(line)
  if (!is.null(.dc_log$file)) cat(line, "\n", file = .dc_log$file, append = TRUE, sep = "")
}

## =============================================================================
## Step 1: the RCTD reference, from any annotated Seurat object
## =============================================================================
# RCTD can't handle "/" in cell-type names, and the names become metadata column names
clean_type_names <- function(x) gsub("[/\\\\]", "-", trimws(as.character(x)))

build_reference <- function(reference_rds, celltype_column, params) {
  log_msg("Reading reference %s", reference_rds)
  ref <- readRDS(reference_rds)
  if (!inherits(ref, "Seurat")) stop("The reference .rds must hold a Seurat object, not a ", class(ref)[1])
  if (!celltype_column %in% colnames(ref[[]]))
    stop("Column '", celltype_column, "' not found in the reference metadata. Available: ",
         paste(colnames(ref[[]]), collapse = ", "))
  assay <- if (params$reference_assay %in% names(ref@assays)) params$reference_assay else DefaultAssay(ref)
  ref[[assay]] <- JoinLayers(ref[[assay]])
  counts <- LayerData(ref, assay = assay, layer = "counts")
  xs <- head(counts@x, 1e5)
  if (length(xs) && any(xs != round(xs)))
    stop("The reference '", assay, "' counts are not whole numbers. RCTD needs raw UMI counts, not normalized data.")

  labels <- clean_type_names(ref[[celltype_column, drop = TRUE]])
  names(labels) <- colnames(ref)
  labels <- labels[!is.na(labels) & labels != "" & labels != "NA"]

  tab <- table(labels)
  small <- names(tab)[tab < params$min_cells_per_type]
  if (length(small)) log_msg("Dropping cell types with < %d cells: %s", params$min_cells_per_type,
                             paste(sprintf("%s (%d)", small, tab[small]), collapse = ", "))
  labels <- labels[!labels %in% small]
  if (length(unique(labels)) < 2) stop("Fewer than 2 cell types left in the reference after filtering")

  set.seed(params$seed)
  keep <- unlist(lapply(split(names(labels), labels), function(cells)
    if (length(cells) > params$max_cells_per_type) sample(cells, params$max_cells_per_type) else cells), use.names = FALSE)
  labels <- labels[keep]
  counts <- counts[, keep, drop = FALSE]

  type_table <- as.data.frame(sort(table(labels), decreasing = TRUE), stringsAsFactors = FALSE)
  names(type_table) <- c("cell_type", "n_cells_used")
  log_msg("Reference: %d cells, %d genes, %d cell types: %s", ncol(counts), nrow(counts), nrow(type_table),
          paste(sprintf("%s (%d)", type_table$cell_type, type_table$n_cells_used), collapse = ", "))

  cell_types <- factor(labels)
  names(cell_types) <- names(labels)
  list(counts = counts, cell_types = cell_types)
}

## =============================================================================
## RCTD itself. spacexr comes in two versions with different interfaces:
##   - Bioconductor (BiocManager::install("spacexr"), R >= 4.5): createRctd() / runRctd()
##   - GitHub (dmcable/spacexr, the original):                   Reference() / SpatialRNA() /
##                                                                create.RCTD() / run.RCTD()
## run_rctd() uses whichever is installed and returns a spots x cell-types weight matrix.
## =============================================================================
spacexr_api <- function() {
  if (exists("createRctd", envir = asNamespace("spacexr"), inherits = FALSE)) "bioconductor" else "github"
}

run_rctd <- function(counts, xy, ref, params) {
  if (spacexr_api() == "bioconductor") {
    spe <- SpatialExperiment::SpatialExperiment(assays = list(counts = counts),
                                                spatialCoords = as.matrix(xy[colnames(counts), c("x", "y")]))
    ref_se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = ref$counts),
      colData = S4Vectors::DataFrame(cell_type = ref$cell_types, row.names = colnames(ref$counts)))
    rctd <- spacexr::createRctd(spe, ref_se, cell_type_col = "cell_type")
    res  <- spacexr::runRctd(rctd, rctd_mode = params$doublet_mode, max_cores = params$n_cores)
    w <- t(as.matrix(SummarizedExperiment::assay(res, "weights")))   # stored cell types x spots
  } else {
    reference <- spacexr::Reference(counts = ref$counts, cell_types = ref$cell_types,
                                    nUMI = Matrix::colSums(ref$counts))
    puck <- spacexr::SpatialRNA(coords = xy[colnames(counts), c("x", "y")], counts = counts,
                                nUMI = Matrix::colSums(counts))
    rctd <- spacexr::create.RCTD(puck, reference, max_cores = params$n_cores)
    rctd <- spacexr::run.RCTD(rctd, doublet_mode = params$doublet_mode)
    w <- as.matrix(rctd@results$weights)                                 # stored spots x cell types
  }
  w
}

## =============================================================================
## Step 2-3: RCTD on one sample
## =============================================================================
# Spot positions: full-resolution pixel coordinates saved by visium_qc.R, or the image's own
spot_xy <- function(obj) {
  md <- obj[[]]
  if (all(c("pxl_row", "pxl_col") %in% names(md))) {
    return(data.frame(x = md$pxl_col, y = md$pxl_row, row.names = Cells(obj)))
  }
  tc <- GetTissueCoordinates(obj)
  if (all(c("imagerow", "imagecol") %in% names(tc))) return(data.frame(x = tc$imagecol, y = tc$imagerow, row.names = rownames(tc)))
  data.frame(x = tc$y, y = tc$x, row.names = rownames(tc))   # Seurat v5 VisiumV2 stores x = row, y = column
}

deconvolve_sample <- function(obj, ref, ref_genes, params, sample_id) {
  obj[["Spatial"]] <- JoinLayers(obj[["Spatial"]])
  counts <- LayerData(obj, assay = "Spatial", layer = "counts")
  shared <- length(intersect(rownames(counts), ref_genes))
  log_msg("[%s] %d spots; %d genes shared with the reference", sample_id, ncol(counts), shared)
  if (shared < params$min_shared_genes)
    log_msg("[%s] WARNING: only %d shared genes. Do both use gene symbols from the same annotation?", sample_id, shared)

  xy <- spot_xy(obj)[colnames(counts), , drop = FALSE]
  w  <- run_rctd(counts, xy, ref, params)
  props <- sweep(w, 1, pmax(rowSums(w), 1e-12), "/")          # each spot's proportions sum to 1
  props <- as.data.frame(props)
  colnames(props) <- paste0("RCTD_", make.names(colnames(props)))
  type_cols <- colnames(props)
  props$RCTD_dominant <- sub("^RCTD_", "", type_cols[max.col(as.matrix(props[, type_cols]), ties.method = "first")])
  props$RCTD_max_proportion <- apply(as.matrix(props[, type_cols]), 1, max)

  old <- grep("^RCTD_", colnames(obj@meta.data), value = TRUE)    # from an earlier run
  if (length(old)) obj@meta.data <- obj@meta.data[setdiff(colnames(obj@meta.data), old)]
  dropped <- setdiff(Cells(obj), rownames(props))
  if (length(dropped)) log_msg("[%s] RCTD skipped %d spots (too few UMIs for its own cutoff); they get NA", sample_id, length(dropped))
  obj <- AddMetaData(obj, props)
  list(obj = obj, props = props, type_cols = type_cols)
}

## =============================================================================
## Step 4: outputs
## =============================================================================
plot_sample <- function(obj, type_cols, sample_id, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  xy <- spot_xy(obj)
  md <- cbind(xy, obj[[]][rownames(xy), c(type_cols, "RCTD_dominant")])
  md <- md[!is.na(md$RCTD_dominant), ]
  base <- function(g) g + scale_y_reverse() + coord_fixed() + theme_void(base_size = 10)

  p1 <- base(ggplot(md, aes(x, y, colour = RCTD_dominant)) + geom_point(size = 0.8)) +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    labs(title = paste0(sample_id, ": dominant cell type per spot"), colour = NULL)
  ggsave(file.path(out_dir, "dominant_cell_type.png"), p1, width = 7, height = 5.5, dpi = 150)

  long <- do.call(rbind, lapply(type_cols, function(cl)
    data.frame(x = md$x, y = md$y, cell_type = sub("^RCTD_", "", cl), proportion = md[[cl]])))
  nc <- min(3, length(type_cols)); nr <- ceiling(length(type_cols) / nc)
  p2 <- base(ggplot(long, aes(x, y, colour = proportion)) + geom_point(size = 0.5)) +
    facet_wrap(~ cell_type, ncol = nc) + scale_colour_viridis_c(limits = c(0, 1)) +
    labs(title = paste0(sample_id, ": RCTD cell-type proportions"), colour = "proportion") +
    theme(strip.text = element_text(size = 10))
  ggsave(file.path(out_dir, "proportions.png"), p2, width = 4 * nc + 1, height = 3.6 * nr + 0.6, dpi = 150)
}

# Write to a temporary file first, so an interruption can't leave a half-written <sample>.rds
save_rds_safely <- function(obj, f) {
  tmp <- paste0(f, ".part")
  saveRDS(obj, tmp)
  if (file.exists(f)) file.remove(f)
  if (!file.rename(tmp, f)) stop("could not replace ", f, " (the new version is in ", tmp, ")")
  invisible(f)
}

# Replace one sample's row in a CSV, or change some of its columns (same as in the QC scripts)
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
update_csv_row <- function(file, id, values, key = "sample_id") {
  old <- if (file.exists(file)) read.csv(file, colClasses = "character", na.strings = "", check.names = FALSE)
  row <- if (!is.null(old) && id %in% old[[key]]) as.list(old[match(id, old[[key]]), , drop = FALSE]) else setNames(list(id), key)
  for (nm in names(values)) row[[nm]] <- values[[nm]]
  upsert_csv(file, row, key)
}

## =============================================================================
## main()
## =============================================================================
main <- function(qc_dir = QC_DIR, reference_rds = REFERENCE_RDS, celltype_column = CELLTYPE_COLUMN, params = PARAMS) {
  if (!dir.exists(qc_dir)) stop("QC folder not found: ", qc_dir, "\n(current working directory: ", getwd(), ")")
  if (!file.exists(reference_rds)) stop("Reference not found: ", reference_rds, "\n(current working directory: ", getwd(), ")")
  if (!requireNamespace("spacexr", quietly = TRUE))
    stop("Package 'spacexr' (RCTD) is not installed. See the README for installation.")
  if (spacexr_api() == "bioconductor" && !requireNamespace("SpatialExperiment", quietly = TRUE))
    stop("The Bioconductor spacexr needs 'SpatialExperiment': BiocManager::install(\"SpatialExperiment\")")
  qc_dir <- normalizePath(qc_dir, winslash = "/")
  start_log(qc_dir)

  files <- Sys.glob(file.path(qc_dir, "*", "*", "*.rds"))
  files <- files[basename(files) == paste0(basename(dirname(files)), ".rds")]   # <dataset>/<sample>/<sample>.rds
  if (!length(files)) stop("No QC'd samples (<dataset>/<sample>/<sample>.rds) in ", qc_dir, ". Run visium_qc.R first.")
  log_msg("Found %d QC'd sample(s) in %s", length(files), qc_dir)

  log_msg("Using spacexr %s (%s interface)", as.character(utils::packageVersion("spacexr")), spacexr_api())
  reference <- build_reference(reference_rds, celltype_column, params)
  ref_genes <- rownames(reference$counts)

  ok <- character(0); failures <- list()
  for (f in files) {
    sid <- basename(dirname(f)); sdir <- dirname(f)
    res <- tryCatch({
      obj <- readRDS(f)
      set.seed(params$seed)
      d <- deconvolve_sample(obj, reference, ref_genes, params, sid)
      save_rds_safely(d$obj, f)
      tab_file <- file.path(sdir, paste0(sid, "_spots.csv"))
      tab <- if (file.exists(tab_file)) read.csv(tab_file, check.names = FALSE, stringsAsFactors = FALSE) else data.frame(cell = colnames(obj))
      tab <- tab[!grepl("^RCTD_", names(tab))]
      i <- match(tab$cell, rownames(d$props))
      for (cl in names(d$props)) tab[[cl]] <- d$props[[cl]][i]
      write.csv(tab, tab_file, row.names = FALSE)
      unlink(file.path(sdir, "deconvolution"), recursive = TRUE)
      tryCatch(plot_sample(d$obj, d$type_cols, sid, file.path(sdir, "deconvolution")),
               error = function(e) log_msg("[%s] WARNING: plotting failed (%s)", sid, conditionMessage(e)))
      means <- sort(colMeans(d$props[, d$type_cols, drop = FALSE]), decreasing = TRUE)
      txt <- paste(sprintf("%s %.2f", sub("^RCTD_", "", names(means)), means), collapse = ", ")
      update_csv_row(file.path(qc_dir, "summary.csv"), sid,
                     list(deconvolved_at = format(Sys.time(), "%Y-%m-%d %H:%M"), rctd_mean_proportions = txt))
      log_msg("[%s] done (%d of %d spots). Mean proportions: %s", sid, nrow(d$props), ncol(obj), txt)
      TRUE
    }, error = function(e) { log_msg("[%s] FAILED: %s", sid, conditionMessage(e)); conditionMessage(e) })
    if (isTRUE(res)) ok <- c(ok, sid) else failures[[sid]] <- res
  }
  log_msg("Done: %d sample(s) OK, %d failed. Outputs in %s", length(ok), length(failures), qc_dir)
  invisible(list(ok = ok, failures = failures, qc_dir = qc_dir))
}

## =============================================================================
## Run
## =============================================================================
args <- commandArgs(trailingOnly = TRUE)
run <- main(qc_dir          = if (length(args) >= 1) args[1] else QC_DIR,
            reference_rds   = if (length(args) >= 2) args[2] else REFERENCE_RDS,
            celltype_column = if (length(args) >= 3) args[3] else CELLTYPE_COLUMN)
if (length(run$failures) && !interactive()) quit(status = 1)
