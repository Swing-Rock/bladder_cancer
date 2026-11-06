# engine.R — loads your QC scripts and wraps them per modality.
#
# The app does NOT re-implement QC. It sources each script into its own
# environment (skipping only the final "Run" block that calls main()), then
# calls the script's own functions: load_sample(), compute_qc_metrics(),
# compute_thresholds(), apply_*_filters(), detect_doublets(), save_sample(),
# make_plots(), ... If you edit a script, the app picks up the change the next
# time you click Load (it re-sources when the file's modified time changes).

APP_DIR <- normalizePath(getwd(), winslash = "/")

# Where the scripts live. Defaults assume the repo layout:
#   qc_explorer/
#   ├── R_scripts/scRNA/scrna_qc.R
#   ├── R_scripts/visium/visium_qc.R
#   ├── R_scripts/xenium/xenium_qc.R
#   ├── R_scripts/annotation/annotate_cells.R
#   └── app_scripts/                  <- this app
# Override any of them with environment variables (e.g. in .Renviron):
#   QCAPP_SCRIPTS_DIR (the R_scripts folder), QCAPP_SCRNA_SCRIPT, QCAPP_VISIUM_SCRIPT,
#   QCAPP_XENIUM_SCRIPT, QCAPP_ANNOTATION_SCRIPT
qc_script_path <- function(modality) {
  base <- Sys.getenv("QCAPP_SCRIPTS_DIR", file.path(dirname(APP_DIR), "R_scripts"))
  p <- switch(modality,
    scrna  = Sys.getenv("QCAPP_SCRNA_SCRIPT",  file.path(base, "scRNA", "scrna_qc.R")),
    visium = Sys.getenv("QCAPP_VISIUM_SCRIPT", file.path(base, "visium", "visium_qc.R")),
    xenium = Sys.getenv("QCAPP_XENIUM_SCRIPT", file.path(base, "xenium", "xenium_qc.R")),
    annotation = Sys.getenv("QCAPP_ANNOTATION_SCRIPT", file.path(base, "annotation", "annotate_cells.R")))
  normalizePath(p, winslash = "/", mustWork = FALSE)
}

# Source a script's definitions without running it: every top-level expression
# is evaluated except the ones that read command-line arguments, call main(),
# or quit R (the "Run" block at the bottom).
QC_FUNCTIONS <- c("DEFAULTS", "resolve_sample", "load_sample", "compute_qc_metrics", "compute_thresholds",
                  "save_sample", "summarise_sample", "make_plots")
ANNOTATION_FUNCTIONS <- c("PARAMS", "REFERENCES", "check_packages", "prepare_references", "annotate_one_sample",
                          "make_sample_plots", "label_columns", "parameters_used")

source_definitions <- function(path, required = QC_FUNCTIONS, defaults_block = "DEFAULTS") {
  if (!file.exists(path)) stop("QC script not found: ", path)
  env <- new.env(parent = globalenv())
  for (e in parse(path, keep.source = FALSE)) {
    nm <- all.names(e)
    if (any(c("commandArgs", "quit") %in% nm)) next
    if (is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name("run"))) next
    eval(e, env)
  }
  for (fn in required)
    if (!exists(fn, envir = env, inherits = FALSE)) stop(basename(path), " has no ", fn, "; the app needs it")
  env$.script_path  <- path
  env$.script_mtime <- file.mtime(path)
  env$.param_help   <- default_comments(path, defaults_block)
  env
}

# The "# comment" after each entry of DEFAULTS, used as help text in the app
default_comments <- function(path, block = "DEFAULTS") {
  lines <- readLines(path, warn = FALSE)
  start <- grep(paste0("^", block, "\\s*<-\\s*list\\("), lines)[1]
  if (is.na(start)) return(list())
  end <- start + which(grepl("^\\)", lines[(start + 1):length(lines)]))[1]
  block <- lines[start:end]
  hit <- regmatches(block, regexec("^\\s*([A-Za-z_][A-Za-z0-9_.]*)\\s*=.*?#\\s*(.*)$", block))
  hit <- Filter(length, hit)
  setNames(lapply(hit, `[[`, 3), vapply(hit, `[[`, "", 2))
}

.engine_cache <- new.env()
get_engine <- function(modality) {
  path <- qc_script_path(modality)
  env <- .engine_cache[[modality]]
  if (is.null(env) || !identical(env$.script_path, path) ||
      !identical(env$.script_mtime, file.mtime(path))) {
    env <- if (modality == "annotation") source_definitions(path, ANNOTATION_FUNCTIONS, "PARAMS") else source_definitions(path)
    .engine_cache[[modality]] <- env
  }
  env
}

## =============================================================================
## Per-modality wrappers. Each one knows which script functions to call and
## which metadata columns the sliders act on. Everything else is the script.
## =============================================================================
COMMON_LIVE <- c("nmads_low", "nmads_high", "min_counts", "min_features", "max_counts", "max_features")

# Groups shown under "QC cutoffs -> Advanced settings"; every other group (and anything the
# script adds later) goes in the "Advanced settings" panel.
LOAD_QC_GROUP <- "Gene names, mito & cell calling (applied on Load)"   # QC settings that only take effect when the sample is loaded again
QC_GROUPS <- c("Automatic (MAD) cutoff rules", "Doublets & ambient RNA", "Gene filter", LOAD_QC_GROUP)

MODALITIES <- list(
  scrna = list(
    label = "Single-cell RNA-seq", unit = "cells",
    path_hint = "10x folder (matrix.mtx + barcodes + features), 10x .h5, dense .csv/.tsv(.gz), .rds (Seurat, SingleCellExperiment or a matrix), or .h5ad (AnnData)",
    meta_file = "cell_qc_metadata.csv", spatial = FALSE,
    sliders = list(
      list(id = "counts",   label = "UMIs per cell",            col = "nCount_RNA",   scale = "log", lo = "min_counts",   hi = "max_counts"),
      list(id = "features", label = "Genes per cell",           col = "nFeature_RNA", scale = "log", lo = "min_features", hi = "max_features"),
      list(id = "mito",     label = "% mitochondrial (max)",    col = "percent.mt",   scale = "linear", hi = "max_mito")),
    groups = list(
      "Reading the data"            = c("format", "is_filtered", "raw_path"),
      "Gene names, mito & cell calling (applied on Load)" = c("gene_column", "mito_pattern", "emptydrops_lower", "emptydrops_fdr"),
      "Automatic (MAD) cutoff rules" = c(COMMON_LIVE, "nmads_mito", "mito_floor", "mito_ceiling"),
      "Doublets & ambient RNA"      = c("run_doublets", "doublet_rate", "run_decontx", "max_contamination"),
      "Gene filter"                 = c("min_cells_gene"),
      "Reproducibility"             = c("seed")),
    extra_step = TRUE),
  visium = list(
    label = "Visium", unit = "spots",
    path_hint = "Space Ranger outs/ folder (or a GEO folder with the same files), a Seurat .rds, or an .h5ad with obsm['spatial'] coordinates",
    meta_file = "spot_qc_metadata.csv", spatial = TRUE,
    sliders = list(
      list(id = "counts",   label = "UMIs per spot",         col = "nCount_Spatial",   scale = "log", lo = "min_counts",   hi = "max_counts"),
      list(id = "features", label = "Genes per spot",        col = "nFeature_Spatial", scale = "log", lo = "min_features", hi = "max_features"),
      list(id = "mito",     label = "% mitochondrial (max)", col = "percent.mt",       scale = "linear", hi = "max_mito")),
    groups = list(
      "Reading the data"            = c("format"),
      "Gene names, mito & cell calling (applied on Load)" = c("gene_column", "mito_pattern"),
      "Automatic (MAD) cutoff rules" = c(COMMON_LIVE, "nmads_mito", "mito_floor", "mito_ceiling"),
      "Gene filter"                 = c("min_spots_gene"),
      "Reproducibility"             = c("seed")),
    extra_step = FALSE),
  xenium = list(
    label = "Xenium", unit = "cells",
    path_hint = "Xenium output folder (cell_feature_matrix.h5 or cell_feature_matrix/ + cells.csv.gz), a LoadXenium() .rds, or an .h5ad (e.g. from spatialdata-io)",
    meta_file = "cell_qc_metadata.csv", spatial = TRUE,
    sliders = list(
      list(id = "counts",   label = "Transcripts per cell",      col = "nCount_Xenium",   scale = "log", lo = "min_counts",    hi = "max_counts"),
      list(id = "features", label = "Genes per cell",            col = "nFeature_Xenium", scale = "log", lo = "min_features",  hi = "max_features"),
      list(id = "area",     label = "Cell area (um^2)",          col = "cell_area",       scale = "log", lo = "min_cell_area", hi = "max_cell_area"),
      list(id = "control",  label = "% negative-control counts (max)", col = "pct_control", scale = "linear", hi = "max_pct_control")),
    groups = list(
      "Reading the data"            = c("format"),
      "Automatic (MAD) cutoff rules" = c(COMMON_LIVE, "nmads_area", "min_cell_area", "max_cell_area", "max_pct_control"),
      "Gene filter"                 = c("min_cells_gene"),
      "Reproducibility"             = c("seed")),
    extra_step = FALSE)
)

live_params <- function(modality) MODALITIES[[modality]]$groups[["Automatic (MAD) cutoff rules"]]
load_params <- function(modality) unlist(MODALITIES[[modality]]$groups[c("Reading the data", LOAD_QC_GROUP)], use.names = FALSE)

# .h5ad files: the scripts' resolve_sample() doesn't know the format, so it checks the row as
# if it were an .rds (same path checks) and the app then reads it itself (R/h5ad.R)
is_h5ad <- function(path) grepl("\\.h5ad$", path %||% "", ignore.case = TRUE)
resolve_row <- function(env, row, i, dir) {
  h5ad <- isTRUE(is_h5ad(row$path)) && (is.null(row$format) || is.na(row$format) || row$format %in% c("auto", "h5ad", ""))
  if (h5ad) row$format <- "rds"
  res <- env$resolve_sample(row, i, dir)
  if (h5ad) res$params$format <- "h5ad"
  res
}

# Step 1 of each script: read the data and measure every cell/spot. Nothing is filtered, and
# for scRNA no cells are called yet (prepare_qc() does that when you click Apply QC).
# Also notes what the input looks like (raw droplets, cells only, already QC'd).
load_data <- function(modality, env, p) {
  id <- p$sample_id
  if (modality == "scrna") {
    if (exists("check_packages", envir = env)) env$check_packages(list(p))
    prior <- NULL
    if (identical(p$format, "h5ad")) {                     # read by the app (R/h5ad.R), not the script
      h <- h5ad_scrna(p); m <- h$counts; prior_cols <- h$obs_cols
    } else {
      m <- with_rds_cache(env, p, function(obj) prior <<- prior_metadata(obj), env$load_sample(p))
      prior_cols <- names(prior)
    }
    colnames(m) <- dedup_prefix(colnames(m), id)
    meta <- env$compute_qc_metrics(m, p)
    if (!identical(p$format, "h5ad")) meta <- add_prior_columns(meta, prior, id)
    raw <- env$detect_raw(m, p$emptydrops_lower); raw$n_total <- ncol(m)
    looks <- describe_input(modality, p$path, prior_cols, raw, attr(prior, "assays"))
    return(list(meta = meta, counts = m, background = NULL, n_barcodes = ncol(m),
                matrix_type = "not called yet", looks = looks))
  }
  if (modality == "visium") {
    if (identical(p$format, "h5ad")) { h <- h5ad_visium(p); obj <- h$obj; prior_cols <- h$obs_cols }
    else { obj <- env$load_sample(p); prior_cols <- colnames(obj[[]]) }
    new <- dedup_prefix(colnames(obj), id)
    if (!identical(new, colnames(obj))) obj <- SeuratObject::RenameCells(obj, new.names = new)
    obj <- env$compute_qc_metrics(obj, p)
    n_mito <- attr(obj, "n_mito_genes")
    md <- obj[[]]
    xy <- env$spot_xy(obj)
    md$x <- xy$x; md$y <- xy$y; md$cell <- rownames(md)
    looks <- describe_input(modality, p$path, prior_cols, assays = names(obj@assays))
    return(list(meta = md, obj = obj, n_mito = n_mito, image = visium_image(obj, md), looks = looks))
  }
  if (modality == "xenium") {
    if (exists("check_packages", envir = env) && !identical(p$format, "h5ad")) env$check_packages(list(p))
    x <- if (identical(p$format, "h5ad")) h5ad_xenium(p) else env$load_sample(p)
    new <- dedup_prefix(colnames(x$counts), id)
    if (!identical(new, colnames(x$counts))) {
      colnames(x$counts) <- new; rownames(x$info) <- new
      if (!is.null(x$obj)) x$obj <- SeuratObject::RenameCells(x$obj, new.names = new)
    }
    meta <- env$compute_qc_metrics(x, p)
    has_obj <- !is.null(x$obj)
    looks <- describe_input(modality, p$path, if (has_obj) colnames(x$obj[[]]) else x$obs_cols %||% character(0),
                            assays = if (has_obj) names(x$obj@assays) else character(0))
    return(list(meta = meta, data = x, looks = looks))
  }
}

# What "Apply QC" does before the cutoffs: for scRNA, scrna_qc.R's step 2 (cell calling with
# emptyDrops on a raw matrix; DecontX background from raw_path) and step 3 on the called cells.
# Visium and Xenium need nothing more. Returns the data in the form the rest of the app uses.
prepare_qc <- function(modality, env, d, p) {
  if (modality != "scrna") return(d)
  if (exists("check_packages", envir = env)) env$check_packages(list(p))
  called <- env$call_cells(d$counts, p)
  meta <- env$compute_qc_metrics(called$counts, p)
  keep_cols <- setdiff(intersect(names(d$meta), c("scDblFinder.score", "scDblFinder.class", "decontX_contamination",
                                                   ANNOT_MARKER_COLS)), names(meta))
  for (cl in keep_cols) meta[[cl]] <- d$meta[rownames(meta), cl]
  list(meta = meta, counts = called$counts, background = called$background,
       n_barcodes = called$n_barcodes, matrix_type = called$matrix_type, looks = d$looks)
}

## ---- what the input looks like -----------------------------------------------------------
# Columns that mean the data were QC'd / annotated before (Seurat + our scripts, and Scanpy / cellxgene)
QC_MARKER_COLS    <- c("qc_reason", "scDblFinder.score", "scDblFinder.class", "decontX_contamination",
                       "doublet_score", "predicted_doublet")
ANNOT_MARKER_COLS <- c("cell_label", "singler_label", "cnv_status", "cell_type", "celltypist_majority_voting", "majority_voting")

describe_input <- function(modality, path, prior_cols, raw = NULL, assays = character(0)) {
  qc_file <- grepl("_(qc|annotated)\\.(rds|h5ad)$", basename(path %||% ""), ignore.case = TRUE)
  qc_cols <- intersect(QC_MARKER_COLS, prior_cols)
  an_cols <- intersect(ANNOT_MARKER_COLS, prior_cols)
  if (length(an_cols))
    return(list(kind = "annotated", short = "already annotated",
                text = sprintf("already QC'd and annotated (has %s)", paste(an_cols, collapse = ", "))))
  if (qc_file || length(qc_cols) || "decontX" %in% assays) {
    why <- c(if (qc_file) sprintf("file name %s", basename(path)), if (length(qc_cols)) sprintf("has %s", paste(qc_cols, collapse = ", ")),
             if ("decontX" %in% assays) "has a decontX assay")
    return(list(kind = "qc", short = "already QC'd", text = sprintf("already QC'd (%s)", paste(why, collapse = "; "))))
  }
  if (!is.null(raw) && isTRUE(raw$raw))
    return(list(kind = "raw", short = "raw droplets",
                text = sprintf("a raw droplet matrix: %s of %s barcodes have only a few UMIs (%.0f%%)",
                               format(raw$n_empty, big.mark = ","), format(raw$n_total, big.mark = ","),
                               100 * raw$frac_empty)))
  unit <- if (modality == "visium") "spots" else "cells"
  list(kind = "cells", short = paste(unit, "only"),
       text = sprintf("%s only, not QC'd yet%s", unit,
                      if (any(c("percent.mt", "pct_counts_mt") %in% prior_cols)) " (it has a % mito column, so it may already have been filtered)" else ""))
}

# Cell names: the QC scripts prefix them with the sample ID on loading. An object saved by
# the scripts is already prefixed, so drop the second copy.
dedup_prefix <- function(x, id) {
  pre <- paste0(id, "_")
  dbl <- startsWith(x, paste0(pre, pre))
  x[dbl] <- substring(x[dbl], nchar(pre) + 1)
  x
}

# scRNA .rds input: read the object once, look at its metadata, and hand the same object to the
# script's load_sample() (which would otherwise read the file a second time)
with_rds_cache <- function(env, p, look, expr) {
  if (!identical(p$format, "rds")) return(expr)
  obj <- readRDS(p$path)
  look(obj)
  path <- normalizePath(p$path, winslash = "/", mustWork = FALSE)
  assign("readRDS", function(file, ...) {
    if (identical(normalizePath(file, winslash = "/", mustWork = FALSE), path)) obj else base::readRDS(file, ...)
  }, envir = env)
  on.exit(rm("readRDS", envir = env), add = TRUE)
  force(expr)
}

prior_metadata <- function(obj) {
  md <- if (inherits(obj, "Seurat")) obj[[]] else if (inherits(obj, "SingleCellExperiment"))
    as.data.frame(SummarizedExperiment::colData(obj)) else NULL
  if (is.null(md)) return(NULL)
  attr(md, "assays") <- if (inherits(obj, "Seurat")) names(obj@assays) else character(0)
  md
}

# Keep earlier doublet / ambient / annotation results from an already-processed object
add_prior_columns <- function(meta, prior, id) {
  if (is.null(prior)) return(meta)
  cols <- intersect(c("scDblFinder.score", "scDblFinder.class", "decontX_contamination", ANNOT_MARKER_COLS), names(prior))
  if (!length(cols)) return(meta)
  rn <- dedup_prefix(paste0(id, "_", rownames(prior)), id)
  i <- match(rownames(meta), rn)
  if (all(is.na(i))) i <- match(rownames(meta), rownames(prior))
  for (cl in cols) meta[[cl]] <- prior[[cl]][i]
  meta
}

# The script's own MAD cutoffs for parameters p
auto_thresholds <- function(modality, env, d, p) {
  switch(modality,
    scrna  = { th <- env$compute_thresholds(d$meta, p, d$n_barcodes); th$input_matrix <- d$matrix_type; th },
    visium = env$compute_thresholds(d$meta, p, d$n_mito),
    xenium = env$compute_thresholds(d$meta, p))
}

apply_filters <- function(modality, env, d, th, p) {
  switch(modality,
    scrna  = env$apply_cell_filters(d$meta, th, p),
    visium = env$apply_spot_filters(d$meta, th, p),
    xenium = env$apply_cell_filters(d$meta, th, p))
}

# scRNA steps 5-6 (doublets, ambient RNA) on the cells that pass the cell filters.
# max_contamination is applied afterwards (apply_contamination) so the slow part
# can be cached while you change it.
run_extra_steps <- function(env, d, meta, p) {
  dbl <- env$detect_doublets(meta, d$counts, p)
  p2 <- p; p2$max_contamination <- NA
  amb <- env$correct_ambient(dbl$meta, d$counts, d$background, p2)
  list(meta = amb$meta, decontx = amb$decontx, doublet_threshold = dbl$threshold,
       key = extra_key(meta, p))
}

extra_key <- function(meta, p) list(keep = meta$keep, p = p[c("run_doublets", "doublet_rate", "run_decontx", "seed")])

# Same rule as correct_ambient() in scrna_qc.R
apply_contamination <- function(meta, p) {
  if (is.na(p$max_contamination)) return(meta)
  bad <- meta$keep & !is.na(meta$decontX_contamination) & meta$decontX_contamination > p$max_contamination
  meta$qc_reason[bad] <- "contamination"
  meta$keep[bad] <- FALSE
  meta
}

# Doublet columns of qc_thresholds.csv, computed the same way as qc_one_sample() in scrna_qc.R
add_doublet_columns <- function(th, meta, threshold) {
  th$doublet_threshold <- round(threshold, 4)
  scored <- !is.na(meta$scDblFinder.class)
  th$doublet_pct <- if (any(scored)) round(100 * mean(meta$qc_reason[scored] == "doublet"), 2) else NA_real_
  th
}

# The result object each script's save_sample()/summarise_sample() expects
result_object <- function(modality, d, meta, th, p, extra = NULL) {
  switch(modality,
    scrna = {
      kept <- meta$cell[meta$keep]
      list(sample_id = p$sample_id, dataset = p$dataset, meta = meta, thresholds = th,
           counts  = d$counts[, kept, drop = FALSE],
           decontx = if (is.null(extra$decontx)) NULL else extra$decontx[, intersect(kept, colnames(extra$decontx)), drop = FALSE],
           min_cells_gene = p$min_cells_gene)
    },
    visium = {
      obj <- d$obj
      obj$qc_reason <- meta$qc_reason
      obj$dataset <- p$dataset
      list(sample_id = p$sample_id, dataset = p$dataset, obj = obj, meta = meta, thresholds = th,
           min_spots_gene = p$min_spots_gene)
    },
    xenium = list(sample_id = p$sample_id, dataset = p$dataset, data = d$data, meta = meta, thresholds = th,
                  min_cells_gene = p$min_cells_gene))
}

## =============================================================================
## Visium tissue image (for drawing spots on the H&E in the app)
## =============================================================================
visium_image <- function(obj, md) {
  tryCatch({
    im <- obj@images[[1]]
    img <- im@image
    sf <- im@scale.factors
    if (is.null(img) || !all(c("pxl_col", "pxl_row") %in% names(md))) return(NULL)
    # The script reads the lowres image if there is one, else the hires one; pick the scale that fits
    cand <- unlist(sf[intersect(c("lowres", "hires"), names(sf))])
    fits <- vapply(cand, function(s) max(md$pxl_col, na.rm = TRUE) * s <= ncol(img) * 1.02 &&
                                      max(md$pxl_row, na.rm = TRUE) * s <= nrow(img) * 1.02, TRUE)
    if (!any(fits)) return(NULL)
    s <- max(cand[fits])
    list(raster = as.raster(img), width = ncol(img) / s, height = nrow(img) / s)
  }, error = function(e) NULL)
}
