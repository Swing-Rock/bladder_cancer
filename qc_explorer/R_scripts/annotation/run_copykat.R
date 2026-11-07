#!/usr/bin/env Rscript
###############################################################################
# run_copykat.R — copyKAT malignant-cell calls for QC'd scRNA-seq samples, as a batch
# job (Expanse or any Linux machine)
#
# Usage (from the qc_explorer folder):
#   Rscript R_scripts/annotation/run_copykat.R <qc output folder | its samplesheet.csv> [sample_id ...]
#   Rscript R_scripts/annotation/run_copykat.R --prefetch
#       checks the packages and downloads the SingleR references, then stops. Run it once on
#       a login node, in case the compute nodes can't reach ExperimentHub.
#
# A thin wrapper around annotate_cells.R. It loads that script's functions (without running
# it) and calls its main(), with two settings forced for every sample, whatever PARAMS or the
# sample sheet's annot_* columns say:
#   malignant_method = "copykat"
#   n_cores          = $SLURM_CPUS_PER_TASK (else $N_CORES, else PARAMS$n_cores)
# Everything else (references, cnv_* settings, cnv_max_cells, reuse_previous, ...) comes from
# PARAMS + each sample's annot_* columns, exactly as in annotate_cells.R and the QC app.
#
# Why force malignant_method: after a quick local run with malignant_method = "none",
# samplesheet.csv records annot_malignant_method = none for each sample, and that would
# otherwise make this run skip copyKAT.
#
# The usual split (annotation/README.md, "Running copyKAT on Expanse"):
#   1. locally:   annotate with malignant_method = none   -> labels saved in each <sample>.rds
#   2. copy the QC output folder to Expanse (Globus)
#   3. on Expanse: sbatch expanse.sbatch -> this script. The labels are reused (same cells
#      and settings), so only copyKAT runs. Results go into the same folder as always.
#   4. copy the folder back
#
# Sample IDs after the folder restrict the run to those samples.
# Exit status 1 if any sample failed or ended up without copyKAT calls (see log.txt).
###############################################################################

cli <- commandArgs(trailingOnly = TRUE)
prefetch <- "--prefetch" %in% cli
cli <- setdiff(cli, "--prefetch")

## ---- find and load annotate_cells.R ---------------------------------------
# It sits next to this file; ANNOTATION_SCRIPT overrides that.
this_file <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) normalizePath(f[1], winslash = "/")
  else normalizePath("R_scripts/annotation/run_copykat.R", winslash = "/", mustWork = FALSE)   # sourced in RStudio from qc_explorer/
})
annot_script <- Sys.getenv("ANNOTATION_SCRIPT", file.path(dirname(this_file), "annotate_cells.R"))
if (!file.exists(annot_script)) stop("annotate_cells.R not found: ", annot_script, " (set ANNOTATION_SCRIPT to its path)")

# Its definitions only, not its "Run" block: the same rule as the QC app's source_definitions()
# (skip what reads command-line arguments, quits, or runs main()).
an <- new.env(parent = globalenv())
for (e in parse(annot_script, keep.source = FALSE)) {
  nm <- all.names(e)
  if (any(c("commandArgs", "quit") %in% nm)) next
  if (is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name("run"))) next
  eval(e, an)
}
an$.script_path <- normalizePath(annot_script, winslash = "/")   # so relative paths in PARAMS resolve as usual

## ---- settings ---------------------------------------------------------------
env_int <- function(v) {
  x <- suppressWarnings(as.integer(Sys.getenv(v, "")))
  if (length(x) == 1 && !is.na(x) && x > 0) x else NA_integer_
}
n_cores <- env_int("SLURM_CPUS_PER_TASK")
if (is.na(n_cores)) n_cores <- env_int("N_CORES")
if (is.na(n_cores)) n_cores <- an$PARAMS$n_cores

FORCED <- list(malignant_method = "copykat", n_cores = n_cores)
params <- an$PARAMS
params[names(FORCED)] <- FORCED

# Sample-sheet annot_* columns are applied per sample by row_params(); put the forced
# settings back on top of them. (main() looks row_params up in `an`, so it uses this one.)
orig_row_params <- an$row_params
an$row_params <- function(row, params = an$PARAMS) {
  r <- orig_row_params(row, params)
  r$params[names(FORCED)] <- FORCED
  r
}

## ---- --prefetch: check packages + cache references, then stop ---------------
if (prefetch) {
  an$check_packages(params)
  refs <- an$prepare_references(params)
  message("Packages OK; SingleR references cached: ", if (length(refs)) paste(names(refs), collapse = ", ") else "none needed")
  quit(status = 0)
}

if (!length(cli)) stop("Usage: Rscript run_copykat.R <qc output folder | its samplesheet.csv> [sample_id ...]\n",
                       "       Rscript run_copykat.R --prefetch")
input <- cli[1]
sample_ids <- cli[-1]

# Optional: only the samples named on the command line
if (length(sample_ids)) {
  orig_list_inputs <- an$list_inputs
  an$list_inputs <- function(input) {
    s <- orig_list_inputs(input)
    miss <- setdiff(sample_ids, s$sample_id)
    if (length(miss)) stop("Not in ", input, " (or not QC'd yet): ", paste(miss, collapse = ", "))
    keep <- s$sample_id %in% sample_ids
    structure(s[keep, , drop = FALSE], out_dir = attr(s, "out_dir"), sheet = attr(s, "sheet"))
  }
}

message(sprintf("run_copykat.R | input: %s | samples: %s | cores: %d | R %s | copykat %s",
                input, if (length(sample_ids)) paste(sample_ids, collapse = ", ") else "all",
                n_cores, as.character(getRversion()),
                tryCatch(as.character(packageVersion("copykat")), error = function(e) "NOT INSTALLED")))
if (.Platform$OS.type == "windows") message("Note: on Windows copyKAT runs on 1 core whatever n_cores is")

## ---- run ----------------------------------------------------------------------
t0 <- Sys.time()
run <- an$main(input = input, params = params)

## ---- check every sample actually got copyKAT calls ---------------------------
# A copyKAT error inside annotate_cells.R is logged as a WARNING and the sample is still
# saved (labels only), so look at what was written rather than trusting "Done".
samples <- an$list_inputs(input)
samples <- samples[samples$sample_id %in% run$ok, , drop = FALSE]
no_calls <- character(0)
for (i in seq_len(nrow(samples))) {
  md <- tryCatch(readRDS(samples$path[i])@meta.data, error = function(e) NULL)
  n <- if (!is.null(md) && "copykat_call" %in% names(md)) sum(!is.na(md$copykat_call)) else 0L
  status <- if (!is.null(md) && "cnv_status" %in% names(md)) {
    t <- table(md$cnv_status); paste(sprintf("%s %d", names(t), t), collapse = ", ")
  } else "no cnv_status"
  message(sprintf("  %-25s copyKAT calls: %6d | %s", samples$sample_id[i], n, status))
  if (n == 0) no_calls <- c(no_calls, samples$sample_id[i])
}

message(sprintf("Finished in %.1f min. Outputs and log.txt in %s",
                as.numeric(difftime(Sys.time(), t0, units = "mins")), run$out_dir))
if (length(run$failures)) message("FAILED: ", paste(names(run$failures), collapse = ", "))
if (length(no_calls)) message("No copyKAT calls (skipped or copyKAT failed; see log.txt): ", paste(no_calls, collapse = ", "))
if (length(run$failures) || length(no_calls)) quit(status = 1)
