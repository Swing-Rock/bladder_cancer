#!/usr/bin/env Rscript
###############################################################################
# annotate_cells.R — cell-type annotation (SingleR) and malignant-cell calling
# (copyKAT / inferCNV) for QC'd samples
#
# Usage
#   Rscript annotate_cells.R <qc output folder>      (or that folder's samplesheet.csv)
#     The output folder of scrna_qc.R, visium_qc.R, xenium_qc.R or the QC app. Every sample in
#     its samplesheet.csv is annotated (without a samplesheet.csv: every <dataset>/<sample>/<sample>.rds).
#     Settings: PARAMS below, for every sample; a sample-sheet column annot_<name> changes one
#     setting for that row's sample (e.g. annot_malignant_method = none; several values: a;b).
#   ...or open bladder cancer.Rproj in RStudio, set INPUT below, and click Source.
#   The QC app (app_scripts/) calls the functions in this file directly.
#
# Per sample, the script:
#   1. reads the QC'd Seurat object (raw counts) and log-normalizes it
#   2. labels every cell with the methods in `label_methods`: SingleR against each reference
#      you list (built-in celldex atlases and/or your own annotated single-cell object), and/or
#      marker-gene scores (built-in bladder TME markers and/or your own list). The final label,
#      `cell_label`, comes from `label_priority`'s method, the other filling its gaps.
#   3. optionally clusters the sample and labels whole clusters too, plus the
#      majority per-cell label in each cluster (a check on the per-cell calls)
#   4. scRNA-seq only: calls each cell malignant / normal (`cnv_status`) with the
#      method(s) chosen in `malignant_method` (copyKAT and/or inferCNV for now; see
#      MALIGNANT_METHODS to add one), using the sample's own immune + stromal cells
#      as the normal reference
#   5. writes the results into the QC output folder, per sample:
#        <dataset>/<sample>/<sample>.rds        the labels are added to the QC'd object (raw counts untouched)
#        <dataset>/<sample>/<sample>_cells.csv  the label columns are added to the per-cell table
#        <dataset>/<sample>/annotation/         diagnostic plots
#      and the settings used (annot_* columns of samplesheet.csv), the labels in summary.csv, and log.txt.
#      copyKAT / inferCNV working files go in working_files/<dataset>/<sample>/ (safe to delete).
#
# Interrupted runs: the labels are saved into <sample>.rds right after step 3, before the
# slow malignant-cell step, and each method's result is cached in working_files/. With
# reuse_previous = TRUE, a rerun with the same cells and settings picks up those labels
# and cached results instead of computing them again. Running the QC script (or saving from
# the QC app) on a sample again removes its annotation: annotate it again afterwards.
#
# By data type
#   scRNA-seq  everything above.
#   Xenium     SingleR works (single cells), but only on the panel's genes: a few
#              hundred genes give coarse labels; a custom reference with the same
#              genes works best. No CNV step (a targeted panel can't show CNVs).
#   Visium     a spot holds several cells, so SingleR can only name each spot's
#              dominant type. For proportions use visium/visium_deconvolution.R (RCTD).
#              No CNV step.
#
# Nothing is merged or batch-corrected here, on purpose. SingleR labels each cell
# independently of every other cell, and copyKAT/inferCNV need each tumour compared
# with normal cells from the same sample, so both run per sample. Merge + Harmony after.
#
# What SingleR can't do: name states that aren't in the reference (myCAF / iCAF are
# not in the celldex atlases; use a bladder/tumour single-cell reference for those).
###############################################################################

## =============================================================================
## Settings (used when no command-line arguments are given)
## =============================================================================
# Relative to the repo root (where "bladder cancer.Rproj" is).
INPUT <- "datasets/GSE130001/qc_output"   # <- the QC output folder (or its samplesheet.csv)

PARAMS <- list(
  # labelling methods
  label_methods       = c("SingleR"), # SingleR (reference-based) and/or markers (marker-gene scores); each writes its own columns
  label_priority      = "SingleR",  # when both run: the method whose label fills `cell_label`; the other fills the cells it leaves unlabelled
  # references (SingleR)
  references          = c("HPCA"),  # built-in references (names in REFERENCES below); empty = none. Project default: HPCA only
  label_level         = "main",     # main = broad types (T cells, Fibroblasts, ...); fine = subtypes
  custom_reference    = NA,         # optional: path to an annotated single-cell Seurat .rds (raw counts)
  custom_label_col    = NA,         # its cell-type metadata column
  custom_max_per_type = 500,        # downsample each custom-reference cell type to this many cells (speed)
  custom_min_per_type = 10,         # drop custom-reference cell types with fewer cells
  primary             = NA,         # reference that fills `singler_label`; blank = custom if given, else the first built-in
  # marker genes (markers)
  custom_markers      = NA,         # optional .csv with columns cell_type, gene: your own marker list (relative paths: see repo_path)
  marker_list_mode    = "add",      # add = your types are added to the built-in list (a type with the same name replaces it); replace = use only yours
  marker_min_score    = 0.1,        # a cell's top marker score must reach this (0-1) to get a label ...
  marker_min_delta    = 0.05,       # ... and beat the runner-up type by this much; otherwise `marker_pruned` is NA
  # test data
  assay               = "RNA",      # RNA = raw counts; decontX = ambient-corrected counts from scrna_qc.R (others fall back to the object's default assay)
  min_shared_genes    = 5000,       # warn if a sample and a reference share fewer genes than this (scaled down for small panels)
  # cluster-level labels
  cluster_level       = TRUE,       # also cluster each sample and label whole clusters
  cluster_resolution  = 0.8,        # Seurat FindClusters resolution
  n_pcs               = 30,         # principal components for clustering / UMAP
  n_hvg               = 2000,       # variable genes for clustering
  # malignant cells (scRNA-seq only)
  malignant_method    = "infercnv", # none | copykat | infercnv | copykat+infercnv (the choices come from MALIGNANT_METHODS below); copyKAT is much slower, see COPYKAT_WARNING
  cnv_reference_types = "T.?cell|T_cell|B.?cell|B_cell|Plasma|NK|Macrophage|Monocyte|Myeloid|DC|Dendritic|Neutrophil|Mast|Fibroblast|Endothelial|Smooth.?muscle|Pericyte|CD4|CD8|Treg|Erythro", # regex on cell_label (the final label): normal reference cells
  cnv_min_reference   = 50,         # fewer reference cells -> copyKAT finds its own; inferCNV is skipped
  cnv_score_quantile  = 0.95,       # inferCNV: malignant if CNV score > this quantile of the reference cells' scores ...
  cnv_cor_cutoff      = 0.4,        # ... and CNV profile correlates with the sample's tumour profile above this
  infercnv_gene_order =  "references/hg38_gencode_v27.txt",         # inferCNV gene positions file (gene, chr, start, end); blank = copyKAT's built-in hg38 table
  cnv_max_cells       = 15000,      # larger samples: CNV runs on a random subset of this many cells (copyKAT is slow)
  # compute
  n_cores             = 4,
  seed                = 1234,
  reuse_previous      = TRUE        # reuse saved labels and malignant-cell results when the cells and settings are unchanged; FALSE = compute everything again
)

# Built-in references from the celldex package (downloaded once, then cached by ExperimentHub).
# All human, bulk-sorted cell populations.
REFERENCES <- list(
  HPCA        = function() celldex::HumanPrimaryCellAtlasData(),         # broad: all major tissue + immune types
  Blueprint   = function() celldex::BlueprintEncodeData(),               # broad: stroma, endothelium, immune
  Monaco      = function() celldex::MonacoImmuneData(),                  # immune only, fine T/B/myeloid subsets
  DICE        = function() celldex::DatabaseImmuneCellExpressionData(),  # immune only
  Novershtern = function() celldex::NovershternHematopoieticData()       # haematopoietic only
)
LABEL_METHODS <- c("SingleR", "markers")

# Built-in marker genes for the major cell types of bladder tumours (human gene symbols).
# Canonical lineage markers; a type needs at least 2 of its markers in the data to be scored.
# Add or replace types with `custom_markers` (a .csv with columns cell_type, gene).
MARKERS <- list(
  Epithelial              = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "KRT5", "KRT13", "UPK1B", "UPK2", "UPK3A", "KRT20"),
  T_cells                 = c("CD3D", "CD3E", "CD3G", "CD2", "TRAC", "CD247"),
  NK_cells                = c("KLRD1", "KLRF1", "NCAM1", "GNLY", "NKG7"),
  B_cells                 = c("MS4A1", "CD79A", "CD79B", "CD19", "BANK1"),
  Plasma_cells            = c("JCHAIN", "MZB1", "IGHG1", "IGKC"),
  Myeloid                 = c("CD14", "LYZ", "CD68", "CSF1R", "C1QA", "AIF1", "TYROBP"),
  Mast_cells              = c("TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2"),
  Fibroblasts             = c("COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA", "FAP"),
  Endothelial             = c("PECAM1", "VWF", "CDH5", "CLDN5", "KDR", "PLVAP"),
  Smooth_muscle_pericytes = c("ACTA2", "MYH11", "DES", "RGS5", "CNN1", "TAGLN")
)

# Broad lineages, to compare labels from different methods and references (their names differ:
# "T_cells", "CD8+ T-cells", "T cells", ...). First match wins; anything else is compared as is.
LINEAGES <- c(
  "Epithelial"    = "epithel|urothel|keratinocyte|carcinoma|tumou?r|malignant",
  "NK cells"      = "^NK|NK.?cell|natural.?killer",
  "T cells"       = "T.?cell|T_cell|CD4|CD8|Treg|thymocyte",
  "Plasma cells"  = "plasma",
  "B cells"       = "B.?cell|B_cell",
  "Mast cells"    = "mast",
  "Myeloid"       = "macrophage|monocyte|myeloid|DC|dendritic|neutrophil",
  "Endothelial"   = "endothel",
  "Smooth muscle" = "smooth.?muscle|pericyte|myocyte",
  "Fibroblasts"   = "fibroblast|CAF"
)
lineage_of <- function(labels) {
  out <- as.character(labels)
  done <- is.na(out)
  for (nm in names(LINEAGES)) {
    hit <- !done & grepl(LINEAGES[[nm]], out, ignore.case = TRUE)
    out[hit] <- nm; done <- done | hit
  }
  out
}
# The labelling methods a setting asks for, in priority order
label_methods <- function(params) {
  m <- intersect(params$label_methods %||% "SingleR", LABEL_METHODS)
  p <- params$label_priority %||% "SingleR"
  c(intersect(p, m), setdiff(m, p))
}

# Malignant-cell methods. Each entry becomes a choice in `malignant_method` (and in the
# app's dropdown). To add a method:
#   1. write run_<method>(counts, ref_cells, sample_id, work_dir, params) returning NULL
#      (skipped) or list(call = named vector over the tested cells, "malignant" / "normal"
#      / NA; columns = data.frame with one row per tested cell (row names = cells) holding
#      the method's own output; heatmap = path to a figure or NA)
#   2. add an entry here. `columns` lists the metadata columns it adds, with their blank
#      value; `params` lists the PARAMS its result depends on (a change re-runs it).
#   3. optional: add a combination to MALIGNANT_COMBOS (cells where the methods disagree
#      become "uncertain").
MALIGNANT_METHODS <- list(
  copykat = list(
    label    = "copyKAT",
    packages = c(copykat = 'remotes::install_github("navinlabcode/copykat")'),
    params   = c("cnv_min_reference"),
    columns  = list(copykat_call = NA_character_),
    run      = function(...) run_copykat(...)),
  infercnv = list(
    label    = "inferCNV",
    packages = c(infercnv = 'BiocManager::install("infercnv")  # also needs JAGS: https://mcmc-jags.sourceforge.io'),
    params   = c("cnv_min_reference", "cnv_score_quantile", "cnv_cor_cutoff", "infercnv_gene_order"),
    columns  = list(cnv_score = NA_real_, cnv_cor = NA_real_, infercnv_call = NA_character_),
    run      = function(...) run_infercnv(...))
)
MALIGNANT_COMBOS <- c("copykat+infercnv")

# Shown (log + the QC app) whenever copyKAT is chosen
COPYKAT_WARNING <- "copyKAT is resource-intensive: it can take an hour or more per sample (~10k cells) and use most of a laptop's memory, and on Windows it runs on a single core. For large datasets, run it on Expanse with R_scripts/annotation/run_copykat.R (see the annotation README)."

# Choices for `malignant_method`, named for display ("copyKAT", "copyKAT + inferCNV", ...)
malignant_choices <- function() {
  lab <- function(x) paste(vapply(strsplit(x, "+", fixed = TRUE)[[1]], function(m) MALIGNANT_METHODS[[m]]$label, ""), collapse = " + ")
  ch <- c("none", names(MALIGNANT_METHODS), MALIGNANT_COMBOS)
  setNames(ch, c("none", vapply(ch[-1], lab, "")))
}
# The methods a setting asks for ("copykat+infercnv" -> c("copykat", "infercnv")).
# Also reads the old setting name `cnv_method` and its value "both".
malignant_methods <- function(params) {
  v <- params$malignant_method %||% params$cnv_method %||% "none"
  if (length(v) != 1 || is.na(v) || !nzchar(v) || v == "none") return(character(0))
  if (v == "both") v <- "copykat+infercnv"
  unique(trimws(strsplit(v, "+", fixed = TRUE)[[1]]))
}
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

# File paths in PARAMS (custom_reference, infercnv_gene_order) can be absolute or relative. A
# relative path is looked up in the folder R is running in, then in qc_explorer/, then in the
# folder that holds it ("bladder cancer", where bladder cancer.Rproj is). So one setting works
# both for this script, run from the project folder, and for the QC app, which runs from
# qc_explorer/app_scripts.
path_roots <- function() {
  sp <- get0(".script_path", inherits = TRUE)   # set by the QC app when it loads this file
  qc <- if (!is.null(sp)) normalizePath(file.path(dirname(sp), "..", ".."), winslash = "/", mustWork = FALSE)
  unique(c(normalizePath(getwd(), winslash = "/"), qc, if (!is.null(qc)) dirname(qc)))
}
repo_path <- function(p) {
  if (length(p) != 1 || is.na(p) || !nzchar(p) || grepl("^([A-Za-z]:|/|\\\\|~)", p) || file.exists(p)) return(p)
  for (root in path_roots()) { cand <- file.path(root, p); if (file.exists(cand)) return(cand) }
  p
}
path_not_found <- function(setting, p)
  stop(setting, " not found: ", p, "\n(looked in: ", paste(path_roots(), collapse = "; "),
       ". A full path such as C:/.../hg38_gencode_v27.txt always works)")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

## =============================================================================
## Logging: console + <QC output folder>/log.txt (added to, never wiped)
## =============================================================================
.an_log <- new.env()
start_log <- function(out_dir, what = "") {
  .an_log$file <- file.path(normalizePath(out_dir, winslash = "/"), "log.txt")   # absolute: copyKAT changes the working directory
  cat(sprintf("\n==== %s  annotate_cells.R %s ====\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), what),
      file = .an_log$file, append = TRUE)
}
log_msg <- function(...) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(...))
  message(line)
  if (!is.null(.an_log$file)) cat(line, "\n", file = .an_log$file, append = TRUE, sep = "")
}

## =============================================================================
## Checks and inputs
## =============================================================================
check_packages <- function(params) {
  methods <- malignant_methods(params)
  unknown_m <- setdiff(methods, names(MALIGNANT_METHODS))
  if (length(unknown_m)) stop("malignant_method: unknown method(s) ", paste(unknown_m, collapse = ", "),
                              ". Choose from: ", paste(malignant_choices(), collapse = ", "))
  bad_lm <- setdiff(params$label_methods, LABEL_METHODS)
  if (length(bad_lm)) stop("label_methods: unknown method(s) ", paste(bad_lm, collapse = ", "), ". Choose from: ", paste(LABEL_METHODS, collapse = ", "))
  lm <- label_methods(params)
  if (!length(lm)) stop("No labelling method: set label_methods to SingleR, markers, or both")
  use_singler <- "SingleR" %in% lm
  need <- c("pheatmap", "viridis", "scales")
  if (use_singler) need <- c(need, "SingleR", "BiocParallel", "SummarizedExperiment")
  if (use_singler && length(params$references)) need <- c(need, "celldex")
  if (use_singler && !is.na(params$custom_reference)) need <- c(need, "scran")   # Wilcoxon marker detection for single-cell references
  how <- unlist(lapply(unname(MALIGNANT_METHODS[methods]), `[[`, "packages"))
  need <- unique(c(need, names(how)))
  missing_pkgs <- need[!vapply(need, requireNamespace, TRUE, quietly = TRUE)]
  if (length(missing_pkgs)) {
    other <- setdiff(missing_pkgs, names(how))
    stop("Missing packages: ", paste(missing_pkgs, collapse = ", "), "\nInstall with:",
         if (length(other)) paste0("\n  BiocManager::install(c(", paste0('"', other, '"', collapse = ", "), "))"),
         paste0("\n  ", how[intersect(missing_pkgs, names(how))], collapse = ""),
         "\nor change `label_methods` / `references` / `malignant_method`.")
  }
  if (use_singler) {
    unknown <- setdiff(params$references, names(REFERENCES))
    if (length(unknown)) stop("Unknown reference(s): ", paste(unknown, collapse = ", "),
                              ". Choose from: ", paste(names(REFERENCES), collapse = ", "))
    if (!params$label_level %in% c("main", "fine")) stop("label_level must be 'main' or 'fine'")
    if (!is.na(params$custom_reference)) {
      if (!file.exists(repo_path(params$custom_reference))) path_not_found("custom_reference", params$custom_reference)
      if (is.na(params$custom_label_col)) stop("Set custom_label_col: the cell-type column of the custom reference")
    }
    if (!length(params$references) && is.na(params$custom_reference))
      stop("No references for SingleR: list at least one in `references`, set `custom_reference`, or use label_methods = markers only")
  }
  if ("markers" %in% lm) {
    if (!params$marker_list_mode %in% c("add", "replace")) stop("marker_list_mode must be 'add' or 'replace'")
    if (!is.na(params$custom_markers) && !file.exists(repo_path(params$custom_markers))) path_not_found("custom_markers", params$custom_markers)
    if (params$marker_list_mode == "replace" && is.na(params$custom_markers))
      stop("marker_list_mode = 'replace' needs a custom_markers file")
  }
  if ("infercnv" %in% methods && !is.na(params$infercnv_gene_order) && !file.exists(repo_path(params$infercnv_gene_order)))
    path_not_found("infercnv_gene_order", params$infercnv_gene_order)
}

# The samples of a QC output folder (the folder or its samplesheet.csv): one row per sample with
# sample_id, dataset, path (its <sample>.rds), patient and the sheet's annot_* columns.
# Without a samplesheet.csv, every <dataset>/<sample>/<sample>.rds in the folder.
list_inputs <- function(input) {
  if (!file.exists(input)) stop("Input not found: ", input, "\n(current working directory: ", getwd(), ")")
  out_dir <- normalizePath(if (dir.exists(input)) input else dirname(input), winslash = "/")
  sheet_file <- if (dir.exists(input)) file.path(out_dir, "samplesheet.csv") else input
  if (file.exists(sheet_file)) {
    sh <- read.csv(sheet_file, colClasses = "character", na.strings = c("", "NA"), check.names = FALSE, strip.white = TRUE)
    names(sh) <- tolower(trimws(names(sh)))
    if (!"sample_id" %in% names(sh)) stop(sheet_file, " has no sample_id column")
    sh <- sh[!is.na(sh$sample_id), , drop = FALSE]
    if (!"dataset" %in% names(sh)) sh$dataset <- NA_character_
    sh$dataset[is.na(sh$dataset)] <- "dataset1"
  } else {
    f <- Sys.glob(file.path(out_dir, "*", "*", "*.rds"))
    f <- f[basename(f) == paste0(basename(dirname(f)), ".rds")]
    sh <- data.frame(sample_id = basename(dirname(f)), dataset = basename(dirname(dirname(f))), stringsAsFactors = FALSE)
  }
  pc <- intersect(c("patient", "patient_id"), names(sh))
  out <- data.frame(sample_id = sh$sample_id, dataset = sh$dataset,
                    path = file.path(out_dir, sh$dataset, sh$sample_id, paste0(sh$sample_id, ".rds")),
                    patient = if (length(pc)) sh[[pc[1]]] else NA_character_, stringsAsFactors = FALSE)
  for (cl in grep("^annot_", names(sh), value = TRUE)) out[[cl]] <- sh[[cl]]
  missing <- !file.exists(out$path)
  if (all(missing)) stop("No QC'd samples in ", out_dir, " (looked for <dataset>/<sample>/<sample>.rds",
                         if (file.exists(sheet_file)) " for the samples in samplesheet.csv", "). Run the QC script or the QC app first.")
  if (any(missing)) log_msg("Not QC'd yet, skipped: %s", paste(out$sample_id[missing], collapse = ", "))
  structure(out[!missing, , drop = FALSE], out_dir = out_dir, sheet = if (file.exists(sheet_file)) sheet_file)
}

# One sample's settings: `params` with the sample sheet's annot_<name> columns on top.
# Several values (label_methods, references) are separated by ";"; "none" = an empty list.
MULTI_VALUE_PARAMS <- c("label_methods", "references")
row_params <- function(row, params = PARAMS) {
  problems <- character(0)
  for (nm in names(PARAMS)) {
    col <- paste0("annot_", nm)
    if (!col %in% names(row)) next
    v <- as.character(row[[col]])
    if (!length(v) || is.na(v) || !nzchar(trimws(v))) next
    v <- trimws(v); d <- PARAMS[[nm]]
    val <- if (nm %in% MULTI_VALUE_PARAMS || length(d) > 1) {
      if (tolower(v) == "none") character(0) else trimws(strsplit(v, ";", fixed = TRUE)[[1]])
    } else if (is.logical(d) && !is.na(d)) {
      if (toupper(v) %in% c("TRUE", "T", "YES", "Y", "1")) TRUE else if (toupper(v) %in% c("FALSE", "F", "NO", "N", "0")) FALSE else NA
    } else if (is.numeric(d)) suppressWarnings(as.numeric(v)) else v
    if (length(val) == 1 && is.na(val)) { problems <- c(problems, sprintf("[%s] %s: can't read '%s'", row$sample_id, col, v)); next }
    params[[nm]] <- val
  }
  list(params = params, problems = problems)
}

# The settings as samplesheet.csv columns (annot_<name>), the way row_params() reads them back
annot_sheet_values <- function(params) {
  vals <- lapply(names(PARAMS), function(nm) {
    v <- params[[nm]]
    if (nm %in% MULTI_VALUE_PARAMS && !length(v)) return("none")
    if (is.null(v) || !length(v) || all(is.na(v))) return(NA_character_)
    paste(v, collapse = ";")
  })
  setNames(vals, paste0("annot_", names(PARAMS)))
}

bp_param <- function(n) {
  if (n <= 1) return(BiocParallel::SerialParam())
  if (.Platform$OS.type == "windows") BiocParallel::SnowParam(n) else BiocParallel::MulticoreParam(n)
}

# What kind of data an object holds, from its assays
data_type <- function(obj) {
  a <- names(obj@assays)
  if ("Spatial" %in% a) return("visium")
  if ("Xenium" %in% a) return("xenium")
  "scrna"
}

## =============================================================================
## Step 0: references (prepared once, used for every sample)
## =============================================================================
# Returns list(name = list(expr = log-expression matrix, labels = character, de_method = ...))
load_builtin_reference <- function(name, level) {
  log_msg("Loading built-in reference %s (first use downloads it; later runs use the cache)", name)
  se <- REFERENCES[[name]]()
  labels <- SummarizedExperiment::colData(se)[[paste0("label.", level)]]
  list(expr = SummarizedExperiment::assay(se, "logcounts"), labels = as.character(labels),
       de_method = "classic", source = sprintf("celldex %s (label.%s)", name, level))
}

load_custom_reference <- function(path, label_col, params) {
  log_msg("Loading custom reference %s (labels: %s)", path, label_col)
  ref <- readRDS(path)
  if (!inherits(ref, "Seurat")) stop("custom_reference must be a Seurat object, not ", class(ref)[1])
  if (!label_col %in% colnames(ref[[]]))
    stop("Column '", label_col, "' not in the custom reference. Available: ", paste(colnames(ref[[]]), collapse = ", "))
  assay <- if ("RNA" %in% names(ref@assays)) "RNA" else DefaultAssay(ref)
  ref[[assay]] <- JoinLayers(ref[[assay]])
  counts <- LayerData(ref, assay = assay, layer = "counts")
  labels <- as.character(ref[[label_col, drop = TRUE]])
  ok <- !is.na(labels) & nzchar(labels) & labels != "NA"
  counts <- counts[, ok, drop = FALSE]; labels <- labels[ok]
  tab <- table(labels)
  small <- names(tab)[tab < params$custom_min_per_type]
  if (length(small)) log_msg("Custom reference: dropping types with < %d cells: %s", params$custom_min_per_type,
                             paste(sprintf("%s (%d)", small, tab[small]), collapse = ", "))
  keep <- !labels %in% small
  set.seed(params$seed)
  idx <- unlist(lapply(split(which(keep), labels[keep]), function(i)
    if (length(i) > params$custom_max_per_type) sample(i, params$custom_max_per_type) else i), use.names = FALSE)
  counts <- counts[, idx, drop = FALSE]; labels <- labels[idx]
  log_msg("Custom reference: %d cells, %d types", ncol(counts), length(unique(labels)))
  # Single-cell references: pairwise Wilcoxon marker detection is the SingleR-recommended setting
  list(expr = NormalizeData(counts, verbose = FALSE), labels = labels, de_method = "wilcox",
       source = sprintf("custom %s (%s)", basename(path), label_col))
}

prepare_references <- function(params) {
  refs <- list()
  if (!"SingleR" %in% label_methods(params)) return(structure(refs, primary = NA_character_))   # markers only
  for (nm in params$references) refs[[nm]] <- load_builtin_reference(nm, params$label_level)
  if (!is.na(params$custom_reference))
    refs[["custom"]] <- load_custom_reference(repo_path(params$custom_reference), params$custom_label_col, params)
  primary <- if (!is.na(params$primary)) params$primary else if ("custom" %in% names(refs)) "custom" else names(refs)[1]
  if (!primary %in% names(refs)) stop("primary = '", primary, "' is not one of the references: ", paste(names(refs), collapse = ", "))
  attr(refs, "primary") <- primary
  refs
}

## =============================================================================
## Steps 1-3: SingleR
## =============================================================================
# Step 1: read the QC'd object and log-normalize the counts SingleR will use
load_sample <- function(row, params) {
  obj <- if (inherits(row$obj[[1]], "Seurat")) row$obj[[1]] else readRDS(row$path)
  if (!inherits(obj, "Seurat")) stop(basename(row$path), " holds a ", class(obj)[1], ", not a Seurat object")
  assay <- params$assay
  if (!assay %in% names(obj@assays)) {
    if (assay != "RNA") log_msg("[%s] no '%s' assay; using %s", row$sample_id, assay, DefaultAssay(obj))
    assay <- DefaultAssay(obj)
  }
  if (inherits(obj[[assay]], "Assay5")) obj[[assay]] <- JoinLayers(obj[[assay]])
  counts <- LayerData(obj, assay = assay, layer = "counts")
  type <- data_type(obj)
  log_msg("[%s] %d genes x %d %s (%s, assay %s)", row$sample_id, nrow(counts), ncol(counts),
          if (type == "visium") "spots" else "cells", type, assay)
  if (type == "visium") log_msg("[%s] Visium: SingleR gives each spot's dominant type only; use RCTD for proportions", row$sample_id)
  list(obj = obj, assay = assay, counts = counts, type = type, logexpr = NormalizeData(counts, verbose = FALSE))
}

# Step 2: SingleR per cell against one reference; returns a data.frame (one row per cell) + the raw result
label_cells <- function(logexpr, ref, ref_name, sample_id, params) {
  shared <- length(intersect(rownames(logexpr), rownames(ref$expr)))
  if (shared == 0) stop(sprintf("no genes shared with reference %s (gene symbols vs Ensembl IDs? check gene_column in QC)", ref_name))
  floor_genes <- min(params$min_shared_genes, 0.5 * nrow(logexpr))   # small panels: judge against the panel size
  if (shared < floor_genes)
    log_msg("[%s] WARNING: only %d genes shared with %s (gene symbols vs Ensembl IDs?)", sample_id, shared, ref_name)
  res <- SingleR::SingleR(test = logexpr, ref = ref$expr, labels = ref$labels, de.method = ref$de_method,
                          BPPARAM = bp_param(params$n_cores))
  out <- data.frame(label  = res$labels, pruned = res$pruned.labels,
                    score  = apply(res$scores, 1, max),
                    delta  = if ("delta.next" %in% colnames(res)) res$delta.next else NA_real_,
                    row.names = rownames(res), stringsAsFactors = FALSE)
  names(out) <- paste0("singler_", ref_name, "_", names(out))
  list(table = out, result = res)
}

# Step 3: cluster the sample (standard Seurat steps on this sample only)
cluster_sample <- function(obj, assay, params) {
  DefaultAssay(obj) <- assay
  n_pcs <- min(params$n_pcs, ncol(obj) - 1, nrow(obj) - 1)
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = min(params$n_hvg, nrow(obj)), verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = n_pcs, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = seq_len(n_pcs), verbose = FALSE)
  obj <- FindClusters(obj, resolution = params$cluster_resolution, random.seed = params$seed, verbose = FALSE)
  obj <- RunUMAP(obj, dims = seq_len(n_pcs), seed.use = params$seed, verbose = FALSE)
  obj
}

# Step 3b: label whole clusters (SingleR on cluster-aggregated expression) and the majority per-cell label
label_clusters <- function(logexpr, clusters, ref, ref_name, cell_pruned, params) {
  res <- SingleR::SingleR(test = logexpr, ref = ref$expr, labels = ref$labels, de.method = ref$de_method,
                          clusters = clusters, BPPARAM = bp_param(params$n_cores))
  cl_label <- setNames(res$labels, rownames(res))
  maj <- tapply(cell_pruned, clusters, function(x) {
    x <- x[!is.na(x)]; if (!length(x)) return(c(NA, NA))
    t <- sort(table(x), decreasing = TRUE); c(names(t)[1], round(t[[1]] / length(x), 3))
  })
  maj_label <- vapply(maj, `[`, "", 1); maj_frac <- as.numeric(vapply(maj, `[`, "", 2))
  out <- data.frame(cluster_label = unname(cl_label[as.character(clusters)]),
                    cluster_majority = unname(maj_label[as.character(clusters)]),
                    cluster_majority_frac = unname(setNames(maj_frac, names(maj))[as.character(clusters)]),
                    stringsAsFactors = FALSE)
  names(out) <- paste0("singler_", ref_name, "_", names(out))
  out
}

## =============================================================================
## Steps 2-3, marker genes: score each cell (and cluster) for each type's markers
## =============================================================================
# The marker list used: MARKERS above, plus or replaced by custom_markers (.csv: cell_type, gene)
load_markers <- function(params) {
  mk <- if (params$marker_list_mode == "replace") list() else MARKERS
  src <- if (params$marker_list_mode == "replace") character(0) else "built-in"
  if (!is.na(params$custom_markers)) {
    f <- repo_path(params$custom_markers)
    tab <- read.csv(f, colClasses = "character", check.names = FALSE, strip.white = TRUE)
    names(tab) <- tolower(trimws(names(tab)))
    if (!all(c("cell_type", "gene") %in% names(tab)))
      stop("custom_markers needs columns cell_type and gene; ", basename(f), " has: ", paste(names(tab), collapse = ", "))
    tab <- tab[!is.na(tab$cell_type) & nzchar(tab$cell_type) & !is.na(tab$gene) & nzchar(tab$gene), ]
    custom <- lapply(split(tab$gene, tab$cell_type), unique)
    mk[names(custom)] <- custom
    src <- c(src, sprintf("%s (%d types)", basename(f), length(custom)))
  }
  structure(mk, source = paste(src, collapse = " + "))
}

# Marker scores, 0-1: each marker gene's log-expression is divided by its 99th percentile in this
# sample (capped at 1), so every gene counts equally whatever its expression level; a type's score
# is the mean over its markers. Returns list(scores = types x cells, used = markers found per type).
marker_scores <- function(logexpr, markers, sample_id) {
  used <- lapply(markers, function(g) intersect(unique(g), rownames(logexpr)))
  few <- vapply(used, length, 0L) < 2
  if (any(few)) log_msg("[%s] markers: skipping %s (fewer than 2 of their markers in the data)", sample_id,
                        paste(sprintf("%s (%d)", names(used)[few], vapply(used[few], length, 0L)), collapse = ", "))
  used <- used[!few]
  if (!length(used)) stop("none of the marker types has 2 or more of its genes in the data (gene symbols vs Ensembl IDs?)")
  genes <- unique(unlist(used))
  m <- as.matrix(logexpr[genes, , drop = FALSE])
  top <- apply(m, 1, function(x) { q <- stats::quantile(x, 0.99, names = FALSE); if (q > 0) q else max(x) })
  top[top <= 0] <- 1                                    # gene never detected: all zeros stay zero
  m <- pmin(m / top, 1)
  scores <- t(vapply(used, function(g) colMeans(m[g, , drop = FALSE]), numeric(ncol(m))))
  if (ncol(m) == 1) scores <- matrix(scores, ncol = 1, dimnames = list(names(used), colnames(m)))
  colnames(scores) <- colnames(logexpr)
  list(scores = scores, used = used, scaled = m)
}

# Per cell: the best-scoring type; pruned to NA when its score is below marker_min_score or it
# beats the runner-up by less than marker_min_delta
label_cells_markers <- function(ms, params) {
  sc <- ms$scores
  ord <- apply(sc, 2, order, decreasing = TRUE)
  ord <- matrix(ord, nrow = nrow(sc))
  best <- rownames(sc)[ord[1, ]]
  top <- sc[cbind(ord[1, ], seq_len(ncol(sc)))]
  second <- if (nrow(sc) > 1) sc[cbind(ord[2, ], seq_len(ncol(sc)))] else 0
  delta <- top - second
  pruned <- ifelse(top >= params$marker_min_score & delta >= params$marker_min_delta, best, NA_character_)
  data.frame(marker_label = best, marker_pruned = pruned, marker_score = round(top, 4), marker_delta = round(delta, 4),
             row.names = colnames(sc), stringsAsFactors = FALSE)
}

# Per cluster: the type whose markers the cluster expresses most on average (NA if even that is
# below marker_min_score), plus the majority per-cell marker label
label_clusters_markers <- function(ms, clusters, cell_pruned, params) {
  cl <- as.character(clusters)
  avg <- vapply(split(seq_along(cl), cl), function(i) rowMeans(ms$scaled[, i, drop = FALSE]), numeric(nrow(ms$scaled)))
  avg <- matrix(avg, nrow = nrow(ms$scaled), dimnames = list(rownames(ms$scaled), sort(unique(cl))))
  type_sc <- vapply(ms$used, function(g) colMeans(avg[g, , drop = FALSE]), numeric(ncol(avg)))
  type_sc <- matrix(type_sc, nrow = ncol(avg), dimnames = list(colnames(avg), names(ms$used)))
  cl_label <- apply(type_sc, 1, function(x) if (max(x) >= params$marker_min_score) names(x)[which.max(x)] else NA_character_)
  maj <- tapply(cell_pruned, cl, function(x) {
    x <- x[!is.na(x)]; if (!length(x)) return(c(NA, NA))
    t <- sort(table(x), decreasing = TRUE); c(names(t)[1], round(t[[1]] / length(x), 3))
  })
  maj_label <- vapply(maj, `[`, "", 1); maj_frac <- as.numeric(vapply(maj, `[`, "", 2))
  data.frame(marker_cluster_label = unname(cl_label[cl]),
             marker_cluster_majority = unname(maj_label[cl]),
             marker_cluster_majority_frac = unname(setNames(maj_frac, names(maj))[cl]),
             stringsAsFactors = FALSE)
}

# The final label from the method(s) run, in priority order:
#   cell_label        the priority method's confident label; where it has none, the other method's
#   cell_label_source which method the label came from
#   cell_lineage      cell_label as a broad lineage (LINEAGES), so labels from different methods
#                     and references share one set of names ("T_cells", "CD8+ T-cells" -> "T cells")
#   label_agreement   both methods ran: "agree" / "disagree" (same broad lineage or not, see
#                     LINEAGES), "one method only" (just one made a confident call) or NA (neither)
combine_labels <- function(obj, methods) {
  md <- obj@meta.data
  col_of <- c(SingleR = "singler_label", markers = "marker_pruned")
  lab <- source <- rep(NA_character_, nrow(md))
  for (m in methods) {
    v <- md[[col_of[[m]]]]
    fill <- is.na(lab) & !is.na(v)
    lab[fill] <- v[fill]; source[fill] <- m
  }
  obj$cell_label <- lab
  obj$cell_label_source <- source
  obj$cell_lineage <- lineage_of(lab)   # the same broad names whichever method/reference the label came from
  if (all(c("SingleR", "markers") %in% methods)) {
    a <- lineage_of(md$singler_label); b <- lineage_of(md$marker_pruned)
    obj$label_agreement <- ifelse(is.na(a) & is.na(b), NA_character_,
                           ifelse(is.na(a) | is.na(b), "one method only", ifelse(a == b, "agree", "disagree")))
  }
  obj
}

## =============================================================================
## Step 4: malignant cells from inferred copy-number changes (scRNA-seq only)
## =============================================================================
# Normal reference = cells whose final label (cell_label) matches cnv_reference_types (immune + stroma)
cnv_reference_cells <- function(obj, params) {
  lab <- if ("cell_label" %in% colnames(obj@meta.data)) obj$cell_label else obj$singler_label
  cells <- colnames(obj)[!is.na(lab) & grepl(params$cnv_reference_types, lab, ignore.case = TRUE)]
  cells
}

# CNV is only meaningful on whole-transcriptome single cells
cnv_applicable <- function(s, params) {
  if (!length(malignant_methods(params))) return("malignant_method = none")
  if (s$type != "scrna") return(sprintf("%s data (spots or a targeted panel)", s$type))
  if (nrow(s$counts) < 8000) return(sprintf("only %d genes (needs a whole transcriptome)", nrow(s$counts)))
  if (ncol(s$counts) < 100) return(sprintf("only %d cells", ncol(s$counts)))
  NA_character_
}

# copyKAT: aneuploid (tumour) vs diploid (normal) per cell. Writes its own files to work_dir.
run_copykat <- function(counts, ref_cells, sample_id, work_dir, params) {
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  work_dir <- normalizePath(work_dir, winslash = "/")
  old <- setwd(work_dir); on.exit(setwd(old), add = TRUE)   # copyKAT writes its files to the working directory
  det <- Matrix::rowMeans(counts > 0)
  m <- as.matrix(counts[det >= 0.01, , drop = FALSE])   # copyKAT needs a dense matrix; rare genes add nothing
  norm <- if (length(ref_cells) >= params$cnv_min_reference) ref_cells else ""
  log_msg("[%s] copyKAT on %d cells x %d genes (%s)", sample_id, ncol(m), nrow(m),
          if (identical(norm, "")) "no reference given: copyKAT finds diploid cells itself" else sprintf("%d reference cells", length(norm)))
  # copyKAT reads its built-in gene tables (full.anno, ...) by name, which only works when the
  # package is attached, not just called with copykat::
  if (!"package:copykat" %in% search()) suppressPackageStartupMessages(library(copykat))
  # copyKAT parallelizes with mclapply, which can't use more than one core on Windows
  cores <- if (.Platform$OS.type == "windows") 1 else max(1, params$n_cores)
  if (cores < params$n_cores) log_msg("[%s] copyKAT runs on 1 core on Windows (n_cores is used elsewhere)", sample_id)
  set.seed(params$seed)
  res <- copykat::copykat(rawmat = m, id.type = "S", ngene.chr = 5, win.size = 25, KS.cut = 0.1,
                          sam.name = sample_id, distance = "euclidean", norm.cell.names = norm,
                          output.seg = "FALSE", plot.genes = "FALSE", genome = "hg20",
                          n.cores = cores)
  pred <- res$prediction
  call <- ifelse(grepl("aneuploid", pred$copykat.pred), "aneuploid",
                 ifelse(grepl("diploid", pred$copykat.pred), "diploid", NA_character_))
  out <- setNames(call, pred$cell.names)
  heat <- list.files(work_dir, pattern = paste0("^", sample_id, ".*heatmap.*\\.(jpe?g|png)$"), full.names = TRUE)
  log_msg("[%s] copyKAT: %d aneuploid, %d diploid, %d undefined", sample_id,
          sum(out == "aneuploid", na.rm = TRUE), sum(out == "diploid", na.rm = TRUE), sum(is.na(out)))
  list(call = setNames(ifelse(out == "aneuploid", "malignant", ifelse(out == "diploid", "normal", NA_character_)), names(out)),
       columns = data.frame(copykat_call = unname(out), row.names = names(out), stringsAsFactors = FALSE),
       heatmap = if (length(heat)) heat[1] else NA_character_)
}

# Gene positions for inferCNV: a user file (gene, chr, start, end) or copyKAT's built-in table
infercnv_gene_order <- function(params, genes) {
  if (!is.na(params$infercnv_gene_order)) {
    go <- read.delim(repo_path(params$infercnv_gene_order), header = FALSE, stringsAsFactors = FALSE)[, 1:4]
  } else {
    if (!requireNamespace("copykat", quietly = TRUE))
      stop("inferCNV needs gene positions: set infercnv_gene_order (e.g. hg38_gencode_v27.txt from ",
           "https://data.broadinstitute.org/Trinity/CTAT/cnv/) or install copykat to use its table")
    # copyKAT keeps full.anno inside the package (sysdata.rda), not as a dataset, so
    # utils::data() can't see it; read it from the package namespace instead
    a <- tryCatch(utils::getFromNamespace("full.anno", "copykat"), error = function(e) NULL)
    if (is.null(a) || !all(c("hgnc_symbol", "chromosome_name", "start_position", "end_position") %in% names(a)))
      stop("copykat's gene table has an unexpected format; set infercnv_gene_order to a gene positions file")
    chr <- as.character(a$chromosome_name); chr[chr == "23"] <- "X"; chr[chr == "24"] <- "Y"
    go <- data.frame(a$hgnc_symbol, paste0("chr", chr), a$start_position, a$end_position, stringsAsFactors = FALSE)
  }
  names(go) <- c("gene", "chr", "start", "stop")
  go <- go[!is.na(go$gene) & nzchar(go$gene) & !duplicated(go$gene) & go$gene %in% genes, ]
  chr_rank <- suppressWarnings(as.numeric(sub("^chr", "", go$chr)))
  chr_rank[is.na(chr_rank)] <- 100 + match(go$chr[is.na(chr_rank)], unique(go$chr[is.na(chr_rank)]))
  go <- go[order(chr_rank, go$start), ]
  data.frame(chr = go$chr, start = go$start, stop = go$stop, row.names = go$gene)
}

# Score each cell's inferred CNV profile (expr: genes x cells, centred on 1, from inferCNV):
#   score = mean squared deviation from 1 (how much CNV signal the cell has)
#   cor   = correlation with the average profile of the sample's top-scoring non-reference cells
# Malignant = score above the cnv_score_quantile of the reference cells AND cor > cnv_cor_cutoff
# (the approach of Tirosh et al. 2016 / Puram et al. 2017).
score_cnv <- function(expr, ref_cells, params) {
  score <- colMeans((expr - 1)^2)
  test <- setdiff(colnames(expr), ref_cells)
  top <- test[score[test] >= quantile(score[test], 0.9)]
  profile <- rowMeans(expr[, top, drop = FALSE])
  cors <- as.numeric(cor(expr, profile))
  ref_q <- if (length(ref_cells)) quantile(score[ref_cells], params$cnv_score_quantile) else quantile(score, 0.5)
  call <- ifelse(colnames(expr) %in% ref_cells, "reference",
                 ifelse(score > ref_q & cors > params$cnv_cor_cutoff, "malignant", "normal"))
  data.frame(cnv_score = score, cnv_cor = cors, infercnv_call = call, row.names = colnames(expr))
}

run_infercnv <- function(counts, ref_cells, sample_id, work_dir, params) {
  if (length(ref_cells) < params$cnv_min_reference) {
    log_msg("[%s] inferCNV skipped: only %d reference (immune/stromal) cells", sample_id, length(ref_cells))
    return(NULL)
  }
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  go <- infercnv_gene_order(params, rownames(counts))
  counts <- counts[rownames(go), , drop = FALSE]
  annot <- data.frame(group = ifelse(colnames(counts) %in% ref_cells, "reference", "test"), row.names = colnames(counts))
  log_msg("[%s] inferCNV on %d cells x %d genes (%d reference cells)", sample_id, ncol(counts), nrow(counts), length(ref_cells))
  obj <- infercnv::CreateInfercnvObject(raw_counts_matrix = counts, annotations_file = annot,
                                        gene_order_file = go, ref_group_names = "reference")
  obj <- infercnv::run(obj, cutoff = 0.1, out_dir = work_dir, cluster_by_groups = TRUE, denoise = TRUE,
                       HMM = FALSE, num_threads = max(1, params$n_cores), no_plot = FALSE,
                       output_format = "png", resume_mode = FALSE)
  sc <- score_cnv(obj@expr.data, intersect(ref_cells, colnames(obj@expr.data)), params)
  log_msg("[%s] inferCNV: %d malignant, %d normal (+%d reference)", sample_id,
          sum(sc$infercnv_call == "malignant"), sum(sc$infercnv_call == "normal"), sum(sc$infercnv_call == "reference"))
  heat <- list.files(work_dir, pattern = "^infercnv\\.png$", full.names = TRUE)
  list(call = setNames(ifelse(sc$infercnv_call == "malignant", "malignant", "normal"), rownames(sc)),   # reference cells count as normal
       columns = sc, heatmap = if (length(heat)) heat[1] else NA_character_)
}

# Every metadata column the malignant-cell methods can add (plus the combined cnv_status)
malignant_columns <- function() c(unlist(lapply(unname(MALIGNANT_METHODS), function(d) names(d$columns))), "cnv_status")

# Combine the methods' calls into one column per cell. res: list(<method> = run_<method>() result).
# Every method's columns are always written (blank for methods not run), so the output has
# the same columns whatever was chosen.
combine_calls <- function(cells, ref_cells, res) {
  out <- data.frame(row.names = cells)
  for (m in names(MALIGNANT_METHODS)) {
    tab <- res[[m]]$columns
    for (col in names(MALIGNANT_METHODS[[m]]$columns))
      out[[col]] <- if (!is.null(tab) && col %in% names(tab)) tab[cells, col] else rep(MALIGNANT_METHODS[[m]]$columns[[col]], length(cells))
  }
  calls <- vapply(names(res), function(m) unname(res[[m]]$call[cells]), character(length(cells)))
  calls <- matrix(calls, nrow = length(cells))
  # one method: its call; several: their shared call, "uncertain" if they disagree; no call = "unresolved"
  status <- apply(calls, 1, function(v) if (anyNA(v)) "unresolved" else if (all(v == v[1])) v[1] else "uncertain")
  status[cells %in% ref_cells & status != "malignant"] <- "normal (reference)"
  out$cnv_status <- status
  out
}

# Settings as plain text, so a cache key doesn't depend on how a value was typed (NA vs "", 30L vs 30)
param_text <- function(params, nms) vapply(nms, function(nm) {
  v <- params[[nm]]
  if (is.null(v) || !length(v) || all(is.na(v))) ""
  else paste(if (is.numeric(v)) format(v, digits = 15, trim = TRUE) else as.character(v), collapse = ";")
}, "")

# One method on one sample. Its result is cached in <sample_dir>/<method>/result_cache.rds together
# with what it depended on (cells, reference cells, counts, the method's settings); a later run
# with the same inputs reuses it (reuse_previous = TRUE) instead of running the method again.
run_malignant_method <- function(m, counts, ref_cells, sample_id, sample_dir, params) {
  def <- MALIGNANT_METHODS[[m]]
  dir <- file.path(sample_dir, m)
  cache <- file.path(dir, "result_cache.rds")
  key <- list(method = m, cells = colnames(counts), ref_cells = sort(ref_cells), genes = nrow(counts),
              total = sum(counts), params = param_text(params, c(def$params, "assay", "seed")))
  if (isTRUE(params$reuse_previous) && file.exists(cache)) {
    hit <- tryCatch(readRDS(cache), error = function(e) NULL)
    if (!is.null(hit) && identical(hit$key, key)) {
      log_msg("[%s] %s: reusing the earlier result (same cells and settings; %s)", sample_id, def$label, cache)
      return(hit$result)
    }
  }
  res <- tryCatch(def$run(counts, ref_cells, sample_id, dir, params),
                  error = function(e) { log_msg("[%s] WARNING: %s failed (%s)", sample_id, def$label, conditionMessage(e)); NULL })
  if (!is.null(res)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(list(key = key, result = res), cache)
  }
  res
}

# Step 4: malignant vs normal with the method(s) in `malignant_method`, using the SingleR labels
# already in `obj` to pick the normal reference cells. Also run on its own by the QC app,
# after its SingleR step. Returns list(obj = obj with the methods' columns + cnv_status,
# cnv = list(skipped, methods, <method> = each method's result)).
call_malignant <- function(obj, counts, type, sample_id, params, work_dir = tempdir()) {
  set.seed(params$seed)
  methods <- malignant_methods(params)
  cnv <- list(skipped = cnv_applicable(list(type = type, counts = counts), params), methods = methods)
  if (is.na(cnv$skipped)) {
    cells <- colnames(obj)
    if (length(cells) > params$cnv_max_cells) {
      set.seed(params$seed); cells <- sort(sample(cells, params$cnv_max_cells))
      log_msg("[%s] malignant-cell calling on a random %d of %d cells (cnv_max_cells)", sample_id, length(cells), ncol(obj))
    }
    ref_cells <- intersect(cnv_reference_cells(obj, params), cells)
    log_msg("[%s] CNV reference: %d immune/stromal cells (cnv_reference_types)", sample_id, length(ref_cells))
    counts <- counts[, cells, drop = FALSE]
    sample_dir <- file.path(work_dir, paste0("cnv_", sample_id))
    res <- list()
    for (m in methods) {
      r <- run_malignant_method(m, counts, ref_cells, sample_id, sample_dir, params)
      if (!is.null(r)) res[[m]] <- r
    }
    cnv <- c(cnv, res)
    if (length(res)) {
      tab <- combine_calls(colnames(obj), ref_cells, res)
      tab$cnv_status[!colnames(obj) %in% cells] <- "not tested"
      obj <- AddMetaData(obj, tab)
      t2 <- table(obj$cnv_status)
      log_msg("[%s] cnv_status: %s", sample_id, paste(sprintf("%s %d", names(t2), t2), collapse = ", "))
    }
  } else if (length(methods)) {
    log_msg("[%s] malignant-cell step skipped: %s", sample_id, cnv$skipped)
  }
  list(obj = obj, cnv = cnv)
}

## =============================================================================
## One sample, start to finish
## =============================================================================
# What the labels (steps 1-3) depend on: the cells and counts, the references and their settings
LABEL_KEY_PARAMS <- c("label_methods", "label_priority", "references", "label_level", "custom_reference", "custom_label_col",
                      "custom_max_per_type", "custom_min_per_type", "primary", "custom_markers", "marker_list_mode",
                      "marker_min_score", "marker_min_delta", "assay", "cluster_level", "cluster_resolution", "n_pcs", "n_hvg", "seed")
labels_key <- function(s, refs, params) list(
  cells = colnames(s$counts), genes = rownames(s$counts), total = sum(s$counts), type = s$type, assay = s$assay,
  params = param_text(params, LABEL_KEY_PARAMS),
  references = vapply(refs, `[[`, "", "source"), primary = attr(refs, "primary"),
  custom_mtime = if (!is.na(params$custom_reference)) as.numeric(file.mtime(repo_path(params$custom_reference))) else NA_real_,
  markers = if ("markers" %in% label_methods(params)) unclass(load_markers(params)))

# A labelled object saved earlier (by this script before its malignant-cell step, or when it
# finished), if its labels came from the same cells and settings; else NULL. obj: the object
# already read from `path`, if it was
read_checkpoint <- function(path, key, obj = NULL) {
  if (is.null(obj) && (is.null(path) || !file.exists(path))) return(NULL)
  if (is.null(obj)) obj <- tryCatch(readRDS(path), error = function(e) NULL)
  if (!inherits(obj, "Seurat") || !identical(obj@misc$annotation$labels_key, key)) return(NULL)
  obj
}

# Steps 1-3 for one sample, with each method in label_methods: per-cell labels (SingleR against
# each reference; marker scores), then clusters and cluster-level labels, then the final
# cell_label in priority order. Returns list(obj, results = the raw SingleR result per reference,
# markers = the marker list used).
label_sample <- function(s, refs, row, params) {
  obj <- s$obj
  methods <- label_methods(params)
  tab_msg <- function(x) { t <- sort(table(x, useNA = "ifany"), decreasing = TRUE); paste(sprintf("%s %d", names(t), t), collapse = ", ") }
  results <- list(); mk <- NULL; ms <- NULL
  if ("SingleR" %in% methods) {
    primary <- attr(refs, "primary")
    for (nm in names(refs)) {
      lc <- label_cells(s$logexpr, refs[[nm]], nm, row$sample_id, params)
      obj <- AddMetaData(obj, lc$table[colnames(obj), , drop = FALSE])
      results[[nm]] <- lc$result
      log_msg("[%s] SingleR %s: %s", row$sample_id, nm, tab_msg(lc$table[[paste0("singler_", nm, "_pruned")]]))
    }
    obj$singler_label <- obj[[paste0("singler_", primary, "_pruned"), drop = TRUE]]
    obj$singler_primary_reference <- primary
  }
  if ("markers" %in% methods) {
    mk <- load_markers(params)
    ms <- marker_scores(s$logexpr, mk, row$sample_id)
    lm <- label_cells_markers(ms, params)
    obj <- AddMetaData(obj, lm[colnames(obj), , drop = FALSE])
    log_msg("[%s] markers (%s; %d types): %s", row$sample_id, attr(mk, "source"), length(ms$used), tab_msg(lm$marker_pruned))
  }

  if (params$cluster_level) {
    if (ncol(obj) < 50) {
      log_msg("[%s] only %d cells; skipping cluster-level labels", row$sample_id, ncol(obj))
    } else {
      obj <- cluster_sample(obj, s$assay, params)
      obj$sample_cluster <- obj$seurat_clusters
      if ("SingleR" %in% methods) for (nm in names(refs)) {
        cl <- label_clusters(s$logexpr[, colnames(obj)], obj$seurat_clusters, refs[[nm]], nm,
                             obj[[paste0("singler_", nm, "_pruned"), drop = TRUE]], params)
        rownames(cl) <- colnames(obj)
        obj <- AddMetaData(obj, cl)
      }
      if ("markers" %in% methods) {
        cl <- label_clusters_markers(ms, obj$seurat_clusters, obj$marker_pruned, params)
        rownames(cl) <- colnames(obj)
        obj <- AddMetaData(obj, cl)
      }
      log_msg("[%s] %d clusters (resolution %.2f)", row$sample_id, length(unique(obj$seurat_clusters)), params$cluster_resolution)
    }
  }
  obj <- combine_labels(obj, methods)
  log_msg("[%s] cell_label (%s): %s", row$sample_id, paste(methods, collapse = " > "), tab_msg(obj$cell_label))
  if ("label_agreement" %in% colnames(obj@meta.data))
    log_msg("[%s] SingleR vs markers: %s", row$sample_id, tab_msg(obj$label_agreement))
  list(obj = obj, results = results, markers = if (!is.null(ms)) ms$used)
}

# One sample, start to finish.
# row: data.frame with sample_id, dataset, path, patient (or an `obj` list-column holding a Seurat object)
# checkpoint: optional .rds path (main() uses the sample's own <sample>.rds). The labelled object
#   (steps 1-3) is saved there before the slow malignant-cell step, so an interrupted run keeps its
#   labels; with reuse_previous = TRUE a labelled object already there (same cells and settings)
#   is used instead of labelling again.
annotate_one_sample <- function(row, refs, params, work_dir = tempdir(), checkpoint = NULL) {
  set.seed(params$seed)
  s <- load_sample(row, params)
  key <- labels_key(s, refs, params)
  same_file <- !is.null(checkpoint) && !is.null(row$path) &&
    identical(normalizePath(checkpoint, winslash = "/", mustWork = FALSE), normalizePath(row$path, winslash = "/", mustWork = FALSE))
  obj <- if (isTRUE(params$reuse_previous)) read_checkpoint(checkpoint, key, if (same_file) s$obj) else NULL
  fresh <- is.null(obj)
  if (!fresh) {
    log_msg("[%s] reusing the labels and clusters saved in %s (same cells and settings)", row$sample_id, checkpoint)
    results <- obj@misc$annotation$singler
    old <- intersect(malignant_columns(), colnames(obj@meta.data))   # recomputed below with the current settings
    if (length(old)) obj@meta.data[old] <- NULL
  } else {
    s$obj <- strip_annotation(s$obj)   # labels from an earlier run (other settings) go
    l <- label_sample(s, refs, row, params)
    obj <- l$obj; results <- l$results
    obj@misc$annotation <- list(labels_key = key, singler = results, markers = l$markers, stage = "labels")
  }
  if (!"sample" %in% colnames(obj[[]])) obj$sample <- row$sample_id
  if (!"dataset" %in% colnames(obj[[]])) obj$dataset <- row$dataset
  if ("patient" %in% names(row) && !is.na(row$patient) && nzchar(row$patient)) obj$patient <- row$patient

  # Save new labels before the slow step, so an interruption there doesn't lose them
  # (reused labels are already in the checkpoint file)
  methods <- malignant_methods(params)
  if (fresh && !is.null(checkpoint)) {
    dir.create(dirname(checkpoint), recursive = TRUE, showWarnings = FALSE)
    obj@misc$annotation$stage <- "labels saved"
    save_rds_safely(drop_scaled(obj), checkpoint)
    if (length(methods))
      log_msg("[%s] labels saved to %s; malignant-cell calling (%s) next", row$sample_id, checkpoint,
              paste(vapply(methods, function(m) MALIGNANT_METHODS[[m]]$label %||% m, ""), collapse = " + "))
  }

  # Step 4: malignant vs normal
  m <- call_malignant(obj, s$counts, s$type, row$sample_id, params, work_dir)
  obj <- m$obj; cnv <- m$cnv
  obj@misc$annotation$stage <- "complete"
  obj@misc$annotation$malignant_method <- params$malignant_method
  list(sample_id = row$sample_id, dataset = row$dataset, obj = obj, results = results, markers = obj@misc$annotation$markers, cnv = cnv, type = s$type)
}

## =============================================================================
## Outputs (into the QC output folder)
## =============================================================================
label_columns <- function(obj) intersect(c("sample", "dataset", "patient", "cell_label", "cell_lineage", "cell_label_source", "label_agreement",
  grep("^singler_|^marker_|^sample_cluster$", colnames(obj[[]]), value = TRUE), malignant_columns()), colnames(obj[[]]))

# The final label column of an annotated object (older objects only have singler_label)
final_label_col <- function(md) if ("cell_label" %in% colnames(md)) "cell_label" else "singler_label"

save_plot <- function(p, file, width, height) ggsave(file, p, width = width, height = height, dpi = 150, limitsize = FALSE)

make_sample_plots <- function(r, plot_dir) {
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  for (nm in names(r$results)) {
    res <- r$results[[nm]]
    # 01: score heatmap: each cell's correlation with every reference label (a clean call = one bright row)
    hm <- tryCatch(SingleR::plotScoreHeatmap(res, show_colnames = FALSE, silent = TRUE,
                                             max.labels = 30, fontsize_row = 8),
                   error = function(e) { log_msg("[%s] score heatmap skipped: %s", r$sample_id, conditionMessage(e)); NULL })
    if (!is.null(hm)) {
      n_lab <- min(30, length(unique(res$labels)))
      png(file.path(plot_dir, sprintf("01_score_heatmap_%s.png", nm)), width = 9, height = max(4, 0.22 * n_lab + 2), units = "in", res = 150)
      grid::grid.newpage(); grid::grid.draw(hm$gtable); dev.off()
    }
    # 02: delta distributions: how far the winning label is ahead of the rest; pruned cells marked
    dd <- tryCatch(SingleR::plotDeltaDistribution(res, ncol = 4), error = function(e) NULL)
    if (!is.null(dd)) save_plot(dd + ggtitle(sprintf("%s: delta per label (orange = pruned / low confidence)", nm)),
                                file.path(plot_dir, sprintf("02_delta_%s.png", nm)), 11, 2.4 * ceiling(length(unique(res$labels)) / 4) + 1)
  }
  obj <- r$obj
  md <- obj[[]]
  if ("umap" %in% Reductions(obj)) {
    dp <- function(col, title, legend = TRUE) {
      g <- DimPlot(obj, group.by = col, label = TRUE, repel = TRUE, label.size = 3) + ggtitle(title)
      if (legend) g else g + NoLegend()
    }
    lab_col <- final_label_col(md)
    pl <- list(dp("seurat_clusters", "Clusters (this sample)", legend = FALSE),
               dp(lab_col, if (lab_col == "cell_label") "Final label (cell_label)" else "Per-cell label (SingleR, pruned)"))
    if ("cell_lineage" %in% colnames(md)) pl <- c(pl, list(dp("cell_lineage", "Final label as broad lineage (cell_lineage)")))
    if ("label_agreement" %in% colnames(md)) pl <- c(pl, list(dp("label_agreement", "SingleR vs markers (same broad lineage?)")))
    primary <- if ("singler_primary_reference" %in% colnames(md)) md$singler_primary_reference[1] else NA
    if (!is.na(primary)) {
      if (lab_col == "cell_label") pl <- c(pl, list(dp("singler_label", sprintf("SingleR per-cell label (%s, pruned)", primary))))
      cl_col <- paste0("singler_", primary, "_cluster_label")
      if (cl_col %in% colnames(md)) pl <- c(pl, list(dp(cl_col, sprintf("SingleR cluster-level label (%s)", primary))))
      for (nm in setdiff(names(r$results), primary))
        pl <- c(pl, list(dp(paste0("singler_", nm, "_pruned"), sprintf("SingleR per-cell label (%s, pruned)", nm))))
    }
    if ("marker_pruned" %in% colnames(md)) pl <- c(pl, list(dp("marker_pruned", "Marker per-cell label (pruned)")))
    if ("marker_cluster_label" %in% colnames(md)) pl <- c(pl, list(dp("marker_cluster_label", "Marker cluster-level label")))
    save_plot(wrap_plots(pl, ncol = 2), file.path(plot_dir, "03_umap_labels.png"), 15, 6 * ceiling(length(pl) / 2))
    if ("cnv_status" %in% colnames(obj[[]])) save_plot(plot_cnv(obj), file.path(plot_dir, "04_cnv_status.png"), 14, 6)
  }
  # 06: the marker genes per final label, to check the marker calls (and the SingleR ones)
  if (length(r$markers)) {
    g <- tryCatch({
      o <- obj
      if (!"data" %in% Layers(o[[DefaultAssay(o)]])) o <- NormalizeData(o, verbose = FALSE)
      o$label_for_plot <- ifelse(is.na(md[[final_label_col(md)]]), "(unlabelled)", md[[final_label_col(md)]])
      DotPlot(o, features = r$markers, group.by = "label_for_plot") + RotatedAxis() +
        labs(title = "Marker genes by final label", x = NULL, y = NULL) +
        theme(axis.text.x = element_text(size = 7), strip.text.x = element_text(size = 7, angle = 90))
    }, error = function(e) { log_msg("[%s] marker dot plot skipped: %s", r$sample_id, conditionMessage(e)); NULL })
    if (!is.null(g)) save_plot(g, file.path(plot_dir, "06_marker_dotplot.png"),
                               max(10, 0.22 * length(unique(unlist(r$markers))) + 4), max(4, 0.35 * length(unique(g$data$id)) + 2.5))
  }
  # Keep copies of the tools' own heatmaps next to the other plots
  for (h in unlist(lapply(names(MALIGNANT_METHODS), function(m) r$cnv[[m]]$heatmap)))
    if (!is.null(h) && !is.na(h) && file.exists(h)) file.copy(h, file.path(plot_dir, paste0("05_", basename(h))), overwrite = TRUE)
}

CNV_COLOURS <- c(malignant = "#C0392B", normal = "#4C9BE8", "normal (reference)" = "grey70",
                 uncertain = "#E69F00", unresolved = "grey40", "not tested" = "grey90")

plot_cnv <- function(obj) {
  g1 <- DimPlot(obj, group.by = "cnv_status", cols = CNV_COLOURS) + ggtitle("Malignant calls from inferred CNVs")
  md <- obj[[]]
  if (all(is.na(md$cnv_score))) return(g1)
  g2 <- ggplot(md, aes(cnv_score, cnv_cor, colour = cnv_status)) + geom_point(size = 0.6, alpha = 0.7) +
    scale_colour_manual(values = CNV_COLOURS) + theme_bw(base_size = 11) +
    labs(title = "inferCNV: CNV score vs correlation with the tumour profile", x = "CNV score", y = "correlation", colour = NULL)
  g1 + g2
}

# What an earlier annotation added to an object (stale once the labels are made again)
ANNOTATION_COLUMN_RX <- "^(singler_|marker_)|^(cell_label|cell_lineage|cell_label_source|label_agreement|sample_cluster|seurat_clusters)$|_snn_res\\."
annotation_columns <- function(cols) unique(c(grep(ANNOTATION_COLUMN_RX, cols, value = TRUE), intersect(malignant_columns(), cols)))
strip_annotation <- function(obj) {
  old <- annotation_columns(colnames(obj@meta.data))
  if (length(old)) obj@meta.data <- obj@meta.data[setdiff(colnames(obj@meta.data), old)]
  obj@misc$annotation <- NULL
  obj@reductions <- list(); obj@graphs <- list()
  obj
}

# The scaled data made for clustering is large and not needed afterwards: not saved
drop_scaled <- function(obj) {
  for (a in names(obj@assays)) obj <- tryCatch({
    x <- obj[[a]]
    if (inherits(x, "Assay5")) { if ("scale.data" %in% SeuratObject::Layers(x)) { x$scale.data <- NULL; obj[[a]] <- x } }
    else if (length(x@scale.data)) { x@scale.data <- matrix(numeric(0), 0, 0); obj[[a]] <- x }
    obj
  }, error = function(e) obj)
  obj
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

# Short text for summary.csv and the log: the commonest final labels, and the malignant calls
count_text <- function(x, top = Inf) {
  t <- sort(table(ifelse(is.na(x), "(unlabelled)", as.character(x))), decreasing = TRUE)
  out <- paste(sprintf("%s %d", names(t)[seq_len(min(top, length(t)))], t[seq_len(min(top, length(t)))]), collapse = ", ")
  if (length(t) > top) paste0(out, sprintf(", +%d more", length(t) - top)) else out
}
malignant_text <- function(r) {
  md <- r$obj@meta.data
  if ("cnv_status" %in% colnames(md)) return(count_text(md$cnv_status))
  if (!length(r$cnv$methods)) "not run (malignant_method = none)" else paste("skipped:", r$cnv$skipped %||% "no result")
}

# Everything for one annotated sample, in the QC output folder: the labelled object replaces
# <sample>.rds, the label columns go into <sample>_cells.csv (or _spots.csv), the plots into
# annotation/, the settings into samplesheet.csv (annot_* columns), the labels into summary.csv.
# The QC app calls this too.
write_annotation_outputs <- function(r, out_dir, params) {
  sdir <- file.path(out_dir, r$dataset, r$sample_id)
  dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
  rds <- file.path(sdir, paste0(r$sample_id, ".rds"))
  r$obj@misc$annotation$params <- params
  save_rds_safely(drop_scaled(r$obj), rds)

  # label columns into the per-cell table (cells QC removed get blanks)
  lab <- r$obj@meta.data[, setdiff(label_columns(r$obj), c("sample", "dataset", "patient")), drop = FALSE]
  tab_file <- file.path(sdir, paste0(r$sample_id, c("_cells.csv", "_spots.csv")))
  tab_file <- if (any(file.exists(tab_file))) tab_file[file.exists(tab_file)][1] else
    file.path(sdir, paste0(r$sample_id, if (identical(r$type, "visium")) "_spots.csv" else "_cells.csv"))
  tab <- if (file.exists(tab_file)) read.csv(tab_file, check.names = FALSE, stringsAsFactors = FALSE) else data.frame(cell = rownames(lab))
  tab <- tab[setdiff(names(tab), annotation_columns(names(tab)))]
  i <- match(tab$cell, rownames(lab))
  for (cl in names(lab)) tab[[cl]] <- lab[[cl]][i]
  write.csv(tab, tab_file, row.names = FALSE)

  plot_dir <- file.path(sdir, "annotation")
  unlink(plot_dir, recursive = TRUE)
  tryCatch(make_sample_plots(r, plot_dir),
           error = function(e) log_msg("[%s] WARNING: plotting failed (%s)", r$sample_id, conditionMessage(e)))

  sheet <- file.path(out_dir, "samplesheet.csv")
  if (file.exists(sheet) && r$sample_id %in% read.csv(sheet, colClasses = "character", check.names = FALSE)$sample_id)
    update_csv_row(sheet, r$sample_id, annot_sheet_values(params))
  labels <- r$obj@meta.data[[final_label_col(r$obj@meta.data)]]
  update_csv_row(file.path(out_dir, "summary.csv"), r$sample_id,
                 list(dataset = r$dataset, annotated_at = format(Sys.time(), "%Y-%m-%d %H:%M"),
                      cell_labels = count_text(labels, 8), malignant = malignant_text(r)))
  log_msg("[%s] annotation saved in %s (labels: %s)", r$sample_id, sdir, count_text(labels, 8))
  invisible(list(dir = sdir, rds = rds, table = tab_file, plots = plot_dir))
}

# Settings that differ from PARAMS, as text (for the log)
param_changes <- function(params) {
  fmt <- function(v) if (!length(v) || all(is.na(v))) "blank" else paste(v, collapse = ";")
  nm <- names(PARAMS)[vapply(names(PARAMS), function(n) !identical(fmt(PARAMS[[n]]), fmt(params[[n]])), TRUE)]
  if (length(nm)) paste(sprintf("%s = %s", nm, vapply(nm, function(n) fmt(params[[n]]), "")), collapse = ", ") else "none"
}

# References take a while to prepare: once per set of reference settings
.refs_cache <- new.env()
references_for <- function(params) {
  key <- paste(deparse(params[c("label_methods", "label_priority", "references", "label_level", "custom_reference",
                                "custom_label_col", "custom_max_per_type", "custom_min_per_type", "primary", "seed")]), collapse = "")
  if (is.null(.refs_cache[[key]])) {
    refs <- prepare_references(params)
    log_msg("Labelling: %s%s", paste(label_methods(params), collapse = " > "),
            if (length(refs)) sprintf("; SingleR references: %s (primary: %s)",
                                      paste(sprintf("%s [%s]", names(refs), vapply(refs, `[[`, "", "source")), collapse = ", "),
                                      attr(refs, "primary")) else "")
    .refs_cache[[key]] <- refs
  }
  .refs_cache[[key]]
}

## =============================================================================
## main(): the whole pipeline, in order
## =============================================================================
main <- function(input = INPUT, params = PARAMS) {
  samples <- list_inputs(input)
  out_dir <- attr(samples, "out_dir")
  start_log(out_dir, input)
  log_msg("%d sample(s) in %d dataset(s); Seurat %s, SingleR %s", nrow(samples), length(unique(samples$dataset)),
          as.character(packageVersion("Seurat")), tryCatch(as.character(packageVersion("SingleR")), error = function(e) "not installed"))

  # Each sample's settings (PARAMS + its annot_* columns); check them all before starting
  rp <- lapply(seq_len(nrow(samples)), function(i) row_params(samples[i, , drop = FALSE], params))
  problems <- unlist(lapply(rp, `[[`, "problems"))
  if (length(problems)) stop(paste(c("Sample sheet problems (nothing was run):", problems), collapse = "\n  - "))
  for (p in unique(lapply(rp, `[[`, "params"))) check_packages(p)
  uses_copykat <- vapply(rp, function(r) "copykat" %in% malignant_methods(r$params), TRUE)
  if (any(uses_copykat)) log_msg("NOTE: copyKAT chosen for %d sample(s). %s", sum(uses_copykat), COPYKAT_WARNING)

  ok <- character(0); failures <- list()
  for (i in seq_len(nrow(samples))) {
    row <- samples[i, , drop = FALSE]; p <- rp[[i]]$params
    log_msg("[%s] settings changed from PARAMS: %s", row$sample_id, param_changes(p))
    res <- tryCatch({
      refs <- references_for(p)
      # <sample>.rds doubles as the checkpoint: the labels are saved there before the malignant-cell step
      r <- annotate_one_sample(row, refs, p, work_dir = file.path(out_dir, "working_files", row$dataset, row$sample_id),
                               checkpoint = row$path)
      write_annotation_outputs(r, out_dir, p)
      TRUE
    }, error = function(e) { log_msg("[%s] FAILED: %s", row$sample_id, conditionMessage(e)); conditionMessage(e) })
    if (isTRUE(res)) ok <- c(ok, row$sample_id) else {
      failures[[row$sample_id]] <- res
      update_csv_row(file.path(out_dir, "summary.csv"), row$sample_id,
                     list(annotated_at = paste("failed", format(Sys.time(), "%Y-%m-%d %H:%M")), cell_labels = res))
    }
    r <- NULL; gc(verbose = FALSE)
  }
  log_msg("Done: %d sample(s) annotated, %d failed. Outputs in %s", length(ok), length(failures), out_dir)
  invisible(list(out_dir = out_dir, ok = ok, failures = failures))
}

## =============================================================================
## Run
## =============================================================================
# Where this file is (Rscript), so relative paths in PARAMS can also be found next to qc_explorer/
.script_path <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
                        if (length(f)) normalizePath(f[1], winslash = "/") })
args <- commandArgs(trailingOnly = TRUE)
run <- main(input = if (length(args) >= 1) args[1] else INPUT)
if (length(run$failures) && !interactive()) quit(status = 1)
