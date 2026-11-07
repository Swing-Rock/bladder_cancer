# export.R — saves one sample from the app, with the QC script's own writer
# (write_sample_outputs()), so the app and the scripts write the same output folder:
#
#   <out_dir>/
#   ├── samplesheet.csv        one row per sample: every setting used, cutoffs set by hand as
#   │                          cutoff_* columns. Load it again (app or script) to reproduce the run.
#   ├── summary.csv            one row per sample: cells in / removed / kept, cutoffs, labels
#   ├── log.txt                added to, never wiped
#   └── <dataset>/<sample>/
#       ├── <sample>.rds        the QC'd object (with the annotation's labels when it was saved too)
#       ├── <sample>_cells.csv  one row per cell (Visium: _spots.csv)
#       ├── qc/  doublets_ambient/  annotation/
#
# Saving writes what the app has for the sample and runs nothing: scDblFinder + DecontX and the
# annotation are included when they were run on the cells being saved, and left out otherwise.

# Every parameter and cutoff that differs from the script's defaults / MAD rule, as readable text
fmt_val <- function(v) if (is.null(v) || length(v) == 0 || (length(v) == 1 && is.na(v))) "blank" else
  if (is.numeric(v)) format(v, digits = 6) else as.character(v)
same_val <- function(a, b) {
  na_a <- is.null(a) || (length(a) == 1 && is.na(a)); na_b <- is.null(b) || (length(b) == 1 && is.na(b))
  if (na_a || na_b) return(na_a && na_b)
  na <- suppressWarnings(as.numeric(a)); nb <- suppressWarnings(as.numeric(b))
  if (!is.na(na) && !is.na(nb)) return(isTRUE(all.equal(na, nb)))
  identical(toupper(as.character(a)), toupper(as.character(b)))
}
param_deviations <- function(env, p, th, th_auto, sources) {
  out <- character(0)
  for (nm in setdiff(names(env$DEFAULTS), "dataset")) {
    d <- env$DEFAULTS[[nm]]
    if (nm == "format" && identical(d, "auto")) next   # "auto" resolves to the detected format; not a choice
    if (!same_val(d, p[[nm]])) out <- c(out, sprintf("%s: default %s -> %s", nm, fmt_val(d), fmt_val(p[[nm]])))
  }
  for (nm in names(sources)[sources == "manual"])
    out <- c(out, sprintf("%s cutoff: MAD rule %s -> set by hand %s", nm, fmt_val(th_auto[[nm]]), fmt_val(th[[nm]])))
  out
}

# annot: list(r = the annotation result for exactly these cells, params = its settings), or NULL
export_sample <- function(modality, env, d, meta, th, p, out_dir, th_auto, sources, extra = NULL,
                          progress = function(...) NULL, sample_meta = list(), annot = NULL, aenv = NULL) {
  mod <- MODALITIES[[modality]]
  out_dir <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(out_dir)) stop("Could not create the output folder: ", out_dir)

  # The scripts' log_msg() writes to log.txt while we save
  log_file <- file.path(out_dir, "log.txt")
  env$.qc_log$file <- log_file
  on.exit(env$.qc_log$file <- NULL, add = TRUE)
  cat(sprintf("\n==== %s  QC app: saving %s ====\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), p$sample_id), file = log_file, append = TRUE)
  cut_txt <- paste(sprintf("%s=%s%s", names(sources), vapply(names(sources), function(nm) format(th[[nm]], digits = 6), ""),
                           ifelse(sources == "manual", sprintf(" (set by hand; MAD rule %s)", vapply(names(sources), function(nm) format(th_auto[[nm]], digits = 6), "")), "")),
                   collapse = ", ")
  env$log_msg("[%s] QC app (%s, script %s): cutoffs %s", p$sample_id, mod$label, basename(env$.script_path), cut_txt)
  dev <- param_deviations(env, p, th, th_auto, sources)
  if (length(dev)) env$log_msg("[%s] QC app: changes from the script's defaults: %s", p$sample_id, paste(dev, collapse = "; "))

  # Cutoffs set by hand go into samplesheet.csv as cutoff_* columns, so the sheet reproduces them
  manual <- names(sources)[sources == "manual"]
  q <- p
  for (nm in manual) q[[paste0("cutoff_", nm)]] <- th[[nm]]
  q$patient <- sample_meta$patient %||% NA_character_

  progress("Saving the Seurat object, table and plots")
  meta <- meta[setdiff(names(meta), ANNOT_MARKER_COLS)]   # labels read with the input: not this QC's (the annotation, if saved, adds its own)
  r <- result_object(modality, d, meta, th, p, extra, sample_meta)
  w <- env$write_sample_outputs(r, q, out_dir, manual = manual, saved_by = "QC app")
  qc_dir <- file.path(w$dir, "qc")
  if (!length(list.files(qc_dir))) fallback_plots(modality, env, meta, th, th_auto, sources, qc_dir, p$sample_id)

  ann <- NULL
  if (!is.null(annot)) {
    progress("Saving the annotation")
    aenv$.an_log$file <- log_file
    on.exit(aenv$.an_log$file <- NULL, add = TRUE)
    a <- annot$r
    a$dataset <- p$dataset; a$obj$dataset <- p$dataset
    for (nm in names(r$sample_meta)) a$obj@meta.data[[nm]] <- r$sample_meta[[nm]]
    ann <- aenv$write_annotation_outputs(a, out_dir, annot$params)
    ann$summary <- annotation_summary(a$obj)
  }
  env$log_msg("[%s] QC app: saved %d of %d %s%s", p$sample_id, sum(meta$keep), nrow(meta), mod$unit,
              if (is.null(ann)) "" else " with the annotation")

  list(dir = w$dir, rds = w$rds, out_dir = out_dir,
       files = c(file.path(out_dir, c("samplesheet.csv", "summary.csv", "log.txt")),
                 list.files(w$dir, recursive = TRUE, full.names = TRUE)),
       deviations = if (length(dev)) paste(dev, collapse = "; ") else "none (all script defaults)", n_deviations = length(dev),
       manual = manual, n_kept = sum(meta$keep), n_in = nrow(meta), annotation = ann)
}

# When the script's make_plots() can't draw a sample (e.g. a Visium .h5ad with no tissue image),
# save the app's own versions of the main QC figures instead, so the sample still has plots
fallback_plots <- function(modality, env, meta, th, th_auto, sources, plot_dir, id) {
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  specs <- MODALITIES[[modality]]$sliders
  srcs <- if (length(sources)) sources else setNames(rep("auto", 0), character(0))
  figs <- list(
    "01_distributions_app.png" = list(function() plot_distributions(meta, specs, th, srcs), 10, 3 * ceiling(length(specs) / 2)),
    "02_counts_vs_genes_app.png" = list(function() plot_scatter(meta, specs, th, reason_colours(env), 60000), 10, 4.5))
  if (isTRUE(MODALITIES[[modality]]$spatial) && all(c("x", "y") %in% names(meta)))
    figs[["03_spatial_qc_app.png"]] <- list(function() plot_spatial(meta, "qc_reason", reason_colours(env), NULL, 150000,
                                                                     if (modality == "xenium") 0.3 else 1, FALSE), 8, 7)
  done <- character(0)
  for (nm in names(figs)) {
    f <- figs[[nm]]
    g <- tryCatch(f[[1]](), error = function(e) NULL)
    if (!is.null(g) && !inherits(tryCatch(ggplot2::ggsave(file.path(plot_dir, nm), g, width = f[[2]], height = f[[3]], dpi = 120),
                                          error = function(e) e), "error")) done <- c(done, nm)
  }
  env$log_msg("[%s] saved the app's own QC plots instead (%s)", id, if (length(done)) paste(done, collapse = ", ") else "none")
  invisible(done)
}
