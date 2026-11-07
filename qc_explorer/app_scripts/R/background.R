# background.R — dataset-level runs in a separate R process (callr), so the app stays
# responsive: you can click through samples while the others are processed.
#
# Each job: the app writes the inputs to a temp file, a fresh R process sources the same
# scripts (engine.R / annotation.R and the lab's QC / annotation script) and runs the step,
# and the app picks up the result when the process finishes. Without callr, jobs run inside
# the app instead (it then can't respond until each sample is done).

HAS_CALLR <- requireNamespace("callr", quietly = TRUE)

# --- run in the child process: self-contained (callr sends only the function) ----------------
bg_setup <- function(app_dir) {
  setwd(app_dir)
  suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(ggplot2); library(patchwork) })
  assign("%||%", function(a, b) if (is.null(a)) b else a, envir = globalenv())
  for (f in c("engine.R", "export.R", "annotation.R")) sys.source(file.path(app_dir, "R", f), envir = globalenv())
}

bg_qc_fun <- function(app_dir, qc_script, modality, in_rds, out_rds, setup) {
  setup(app_dir)
  env <- source_definitions(qc_script)
  x <- readRDS(in_rds)
  saveRDS(prepare_qc(modality, env, x$d, x$p), out_rds, compress = FALSE)
  TRUE
}

bg_extra_fun <- function(app_dir, qc_script, in_rds, out_rds, setup) {
  setup(app_dir)
  env <- source_definitions(qc_script)
  x <- readRDS(in_rds)
  res <- run_extra_steps(env, x$d, x$meta, x$p)
  saveRDS(res, out_rds, compress = FALSE)
  TRUE
}

bg_annot_fun <- function(app_dir, annot_script, qc_rds, sample_id, dataset, sample_meta, params, work_dir, checkpoint,
                         out_rds, setup) {
  setup(app_dir)
  aenv <- source_definitions(annot_script, ANNOTATION_FUNCTIONS, "PARAMS")
  # work_dir and checkpoint are in the app's folders (this child process has its own tempdir)
  r <- annotate_saved_sample(aenv, qc_rds, sample_id, dataset, sample_meta, params, work_dir, checkpoint = checkpoint)
  # keep only what the app needs from the references (they are large)
  r$refs <- structure(lapply(r$refs, function(x) list(source = x$source)), primary = attr(r$refs, "primary"))
  saveRDS(r, out_rds, compress = FALSE)
  TRUE
}

# --- in the app -------------------------------------------------------------------------------
bg_job_dir <- function() { d <- tempfile("qcjob_"); dir.create(d); normalizePath(d, winslash = "/") }

bg_start <- function(fun, args, dir) {
  log <- file.path(dir, "log.txt")
  proc <- callr::r_bg(fun, args = c(args, list(setup = bg_setup)), stdout = log, stderr = "2>&1",
                      supervise = TRUE, package = FALSE)
  list(proc = proc, log = log, dir = dir, started = Sys.time())
}

# Finished? Returns NULL while running, else list(ok, error, log_lines)
bg_poll <- function(job) {
  if (job$proc$is_alive()) return(NULL)
  lines <- if (file.exists(job$log)) readLines(job$log, warn = FALSE) else character(0)
  err <- tryCatch({ job$proc$get_result(); NULL }, error = function(e) {
    m <- conditionMessage(e)
    if (!is.null(e$parent)) m <- conditionMessage(e$parent)   # the child's own error
    m
  })
  list(ok = is.null(err), error = err, log_lines = lines)
}

# ---- jobs that outlive the page ------------------------------------------------------------
# A job's process is killed when nothing refers to it any more, which happens when the page
# that started it disconnects. Annotation jobs are worth finishing (their labels and malignant-cell
# results are cached, so annotating the same cells again reuses them), so the app keeps them
# here, shared by every page, until they end; a page that connects later is told about them.
.app_jobs <- new.env()
bg_keep <- function(job, sample_id, task, label) {
  key <- paste0(sample_id, "_", format(Sys.time(), "%H%M%OS3"))
  assign(key, list(job = job, sample_id = sample_id, task = task, label = label, started = job$started, kept = Sys.time()), envir = .app_jobs)
  key
}
# Kept jobs, with whether each is still running; finished ones are forgotten after `keep_h` hours
bg_kept <- function(keep_h = 6) {
  out <- lapply(ls(.app_jobs), function(k) {
    x <- get(k, envir = .app_jobs)
    alive <- tryCatch(x$job$proc$is_alive(), error = function(e) FALSE)
    if (!alive && is.null(x$ended)) { x$ended <- Sys.time(); assign(k, x, envir = .app_jobs) }
    if (!alive && difftime(Sys.time(), x$ended, units = "hours") > keep_h) { rm(list = k, envir = .app_jobs); return(NULL) }
    c(x[c("sample_id", "task", "label", "started")], list(alive = alive, ended = x$ended))
  })
  Filter(Negate(is.null), out)
}
bg_running_for <- function(sample_id) any(vapply(bg_kept(), function(x) isTRUE(x$alive) && identical(x$sample_id, sample_id), TRUE))

# Kill a background job and anything it started (e.g. copyKAT's worker processes)
bg_kill <- function(job) {
  if (is.null(job$proc) || !job$proc$is_alive()) return(invisible())
  killed <- tryCatch({ job$proc$kill_tree(); TRUE }, error = function(e) FALSE)
  if (!killed || job$proc$is_alive()) job$proc$kill()
  invisible()
}
