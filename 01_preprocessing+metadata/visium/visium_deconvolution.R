#!/usr/bin/env Rscript
###############################################################################
# visium_deconvolution.R — cell-type deconvolution of QC'd Visium spots with RCTD
#
# Usage
#   Rscript visium_deconvolution.R <qc_dir> <reference.rds> <celltype_column> [out_dir]
#     qc_dir           folder with the <sample>_qc.rds files written by visium_qc.R
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
#   4. saves <sample>_deconv.rds, <sample>_proportions.csv and spatial maps
# A summary table compares mean cell-type proportions across samples.
###############################################################################

## =============================================================================
## Settings (used when no command-line arguments are given)
## =============================================================================
# Relative to the repo root (where "bladder cancer.Rproj" is).
QC_DIR          <- "datasets/VISIUM_DATASET/qc_output/VISIUM_DATASET"   # <- folder of <sample>_qc.rds
REFERENCE_RDS   <- "datasets/REFERENCE/reference_annotated.rds"         # <- annotated scRNA-seq Seurat object
CELLTYPE_COLUMN <- "celltype"                                           # <- its cell-type column
OUT_DIR         <- NULL   # NULL = "<QC_DIR>/deconvolution"

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
## Logging: console + <out_dir>/deconvolution_log.txt
## =============================================================================
.dc_log <- new.env()
start_log <- function(out_dir) {
  .dc_log$file <- file.path(out_dir, "deconvolution_log.txt")
  cat("", file = .dc_log$file)
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

build_reference <- function(reference_rds, celltype_column, params, out_dir) {
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
  write.csv(type_table, file.path(out_dir, "reference_cell_types.csv"), row.names = FALSE)
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

  dropped <- setdiff(Cells(obj), rownames(props))
  if (length(dropped)) log_msg("[%s] RCTD skipped %d spots (too few UMIs for its own cutoff); they get NA", sample_id, length(dropped))
  obj <- AddMetaData(obj, props)
  list(obj = obj, props = props, type_cols = type_cols)
}

## =============================================================================
## Step 4: outputs
## =============================================================================
plot_sample <- function(obj, type_cols, sample_id, out_dir) {
  xy <- spot_xy(obj)
  md <- cbind(xy, obj[[]][rownames(xy), c(type_cols, "RCTD_dominant")])
  md <- md[!is.na(md$RCTD_dominant), ]
  base <- function(g) g + scale_y_reverse() + coord_fixed() + theme_void(base_size = 10)

  p1 <- base(ggplot(md, aes(x, y, colour = RCTD_dominant)) + geom_point(size = 0.8)) +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    labs(title = paste0(sample_id, ": dominant cell type per spot"), colour = NULL)
  ggsave(file.path(out_dir, paste0(sample_id, "_dominant_cell_type.png")), p1, width = 7, height = 5.5, dpi = 150)

  long <- do.call(rbind, lapply(type_cols, function(cl)
    data.frame(x = md$x, y = md$y, cell_type = sub("^RCTD_", "", cl), proportion = md[[cl]])))
  nc <- min(3, length(type_cols)); nr <- ceiling(length(type_cols) / nc)
  p2 <- base(ggplot(long, aes(x, y, colour = proportion)) + geom_point(size = 0.5)) +
    facet_wrap(~ cell_type, ncol = nc) + scale_colour_viridis_c(limits = c(0, 1)) +
    labs(title = paste0(sample_id, ": RCTD cell-type proportions"), colour = "proportion") +
    theme(strip.text = element_text(size = 10))
  ggsave(file.path(out_dir, paste0(sample_id, "_proportions.png")), p2, width = 4 * nc + 1, height = 3.6 * nr + 0.6, dpi = 150)
}

## =============================================================================
## main()
## =============================================================================
main <- function(qc_dir = QC_DIR, reference_rds = REFERENCE_RDS, celltype_column = CELLTYPE_COLUMN,
                 out_dir = OUT_DIR, params = PARAMS) {
  if (!dir.exists(qc_dir)) stop("QC folder not found: ", qc_dir, "\n(current working directory: ", getwd(), ")")
  if (!file.exists(reference_rds)) stop("Reference not found: ", reference_rds, "\n(current working directory: ", getwd(), ")")
  if (!requireNamespace("spacexr", quietly = TRUE))
    stop("Package 'spacexr' (RCTD) is not installed. See the README for installation.")
  if (spacexr_api() == "bioconductor" && !requireNamespace("SpatialExperiment", quietly = TRUE))
    stop("The Bioconductor spacexr needs 'SpatialExperiment': BiocManager::install(\"SpatialExperiment\")")
  if (is.null(out_dir)) out_dir <- file.path(qc_dir, "deconvolution")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  start_log(out_dir)

  files <- list.files(qc_dir, pattern = "_qc\\.rds$", full.names = TRUE)
  if (!length(files)) stop("No *_qc.rds files in ", qc_dir, ". Run visium_qc.R first.")
  log_msg("Found %d QC'd sample(s) in %s", length(files), qc_dir)

  log_msg("Using spacexr %s (%s interface)", as.character(utils::packageVersion("spacexr")), spacexr_api())
  reference <- build_reference(reference_rds, celltype_column, params, out_dir)
  ref_genes <- rownames(reference$counts)

  summaries <- list(); failures <- list()
  for (f in files) {
    sid <- sub("_qc\\.rds$", "", basename(f))
    res <- tryCatch({
      obj <- readRDS(f)
      set.seed(params$seed)
      d <- deconvolve_sample(obj, reference, ref_genes, params, sid)
      saveRDS(d$obj, file.path(out_dir, paste0(sid, "_deconv.rds")))
      write.csv(cbind(spot = rownames(d$props), d$props), file.path(out_dir, paste0(sid, "_proportions.csv")), row.names = FALSE)
      tryCatch(plot_sample(d$obj, d$type_cols, sid, out_dir),
               error = function(e) log_msg("[%s] WARNING: plotting failed (%s)", sid, conditionMessage(e)))
      means <- colMeans(d$props[, d$type_cols, drop = FALSE])
      log_msg("[%s] done. Mean proportions: %s", sid,
              paste(sprintf("%s %.2f", sub("^RCTD_", "", names(means)), means), collapse = ", "))
      data.frame(sample = sid, status = "ok", spots = ncol(obj), spots_deconvolved = nrow(d$props), t(means), check.names = FALSE)
    }, error = function(e) { log_msg("[%s] FAILED: %s", sid, conditionMessage(e)); conditionMessage(e) })
    if (is.character(res)) failures[[sid]] <- res else summaries[[sid]] <- res
  }

  for (sid in names(failures)) summaries[[sid]] <- data.frame(sample = sid, status = paste("failed:", failures[[sid]]))
  cols <- unique(unlist(lapply(summaries, names)))
  fill <- function(d) { for (cl in setdiff(cols, names(d))) d[[cl]] <- NA; d[cols] }
  overview <- do.call(rbind, lapply(summaries, fill)); rownames(overview) <- NULL
  write.csv(overview, file.path(out_dir, "deconvolution_summary.csv"), row.names = FALSE)
  log_msg("Done: %d sample(s) OK, %d failed. Outputs in %s", length(files) - length(failures), length(failures), out_dir)
  invisible(list(summary = overview, failures = failures, out_dir = out_dir))
}

## =============================================================================
## Run
## =============================================================================
args <- commandArgs(trailingOnly = TRUE)
run <- main(qc_dir          = if (length(args) >= 1) args[1] else QC_DIR,
            reference_rds   = if (length(args) >= 2) args[2] else REFERENCE_RDS,
            celltype_column = if (length(args) >= 3) args[3] else CELLTYPE_COLUMN,
            out_dir         = if (length(args) >= 4) args[4] else OUT_DIR)
if (length(run$failures) && !interactive()) quit(status = 1)
