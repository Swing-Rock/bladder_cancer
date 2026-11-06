# export.R — writes one sample's QC results in the same layout as the scripts:
#
#   <out_dir>/
#   ├── all_samples_summary.csv      one row per sample (added/updated each time you save)
#   ├── samplesheet_resolved.csv     every parameter value used, per sample (same columns as the scripts')
#   ├── qc_parameters_used.csv       every parameter + final cutoffs, and what differs from the script defaults
#   ├── app_samplesheet_rows.csv     sample-sheet rows that reproduce the app's cutoffs with the script
#   ├── qc_log.txt                   appended, never wiped
#   └── <dataset>/
#       ├── <sample_id>_qc.rds       written by the script's own save_sample()
#       ├── qc_thresholds.csv        cutoffs used (+ the MAD cutoffs and which ones you moved by hand)
#       ├── <sample_id>_<cell|spot>_qc_metadata.csv
#       └── plots/<sample_id>/       the script's own make_plots() figures
#
# The script writes one combined metadata table and one plots/ folder per dataset;
# the app saves one sample at a time, so those are kept per sample here instead.

bind_fill <- function(a, b) {
  cols <- union(names(a), names(b))
  for (cl in setdiff(cols, names(a))) a[[cl]] <- if (nrow(a)) NA else logical(0)
  for (cl in setdiff(cols, names(b))) b[[cl]] <- NA
  rbind(a[cols], b[cols])
}

# Replace this sample's row(s) in a CSV (or create it)
upsert_csv <- function(file, df, key) {
  if (file.exists(file)) {
    old <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character")
    if (key %in% names(old)) old <- old[old[[key]] != as.character(df[[key]][1]), , drop = FALSE]
    df[] <- lapply(df, function(x) if (is.factor(x)) as.character(x) else x)
    df <- bind_fill(old, as.data.frame(lapply(df, as.character), stringsAsFactors = FALSE, check.names = FALSE))
  }
  write.csv(df, file, row.names = FALSE, na = "")
  invisible(file)
}

params_row <- function(env, p) {
  keep <- c("sample_id", "path", names(env$DEFAULTS))
  as.data.frame(lapply(p[keep], function(v) if (is.null(v) || length(v) == 0) NA else v),
                stringsAsFactors = FALSE, check.names = FALSE)
}

# Which sheet parameters set each cutoff when it is "pinned" to a value
PIN_RULES <- list(
  min_counts      = function(p, v) { p$min_counts <- v; p },
  max_counts      = function(p, v) { p$max_counts <- if (is.finite(v)) v else NA; p },
  min_features    = function(p, v) { p$min_features <- v; p },
  max_features    = function(p, v) { p$max_features <- if (is.finite(v)) v else NA; p },
  max_mito        = function(p, v) { p$mito_floor <- v; p$mito_ceiling <- v; p },
  min_cell_area   = function(p, v) { p$min_cell_area <- v; p },
  max_cell_area   = function(p, v) { p$max_cell_area <- if (is.finite(v)) v else NA; p },
  max_pct_control = function(p, v) { p$max_pct_control <- if (is.finite(v)) v else NA; p })

# A sample-sheet row that makes the script reproduce exactly the cutoffs used in the app.
# No manual changes -> the row is just the parameters. Manual changes -> every cutoff is
# pinned to its final value (MAD widened to 99 so it never binds), then checked by
# re-running the script's compute_thresholds().
reproducible_sheet_row <- function(modality, env, d, p, th_used, sources) {
  manual <- names(sources)[sources == "manual"]
  row <- params_row(env, p)
  if (!length(manual)) {
    row$notes <- "QC app: automatic (MAD) cutoffs, no manual changes"
    row$reproduces_app_cutoffs <- TRUE
    return(row)
  }
  q <- p
  for (nm in intersect(c("nmads_low", "nmads_high", "nmads_area"), names(q))) q[[nm]] <- 99
  for (nm in names(sources)) if (!is.null(PIN_RULES[[nm]])) q <- PIN_RULES[[nm]](q, th_used[[nm]])
  chk <- auto_thresholds(modality, env, d, q)
  same <- vapply(names(sources), function(nm) {
    a <- chk[[nm]]; b <- th_used[[nm]]
    (is.infinite(a) && is.infinite(b)) || isTRUE(abs(a - b) <= 1e-6 * max(1, abs(b)))
  }, TRUE)
  row <- params_row(env, q)
  row$notes <- sprintf("QC app: cutoffs pinned (manual: %s). MAD cutoffs disabled (nmads = 99) so the script uses these exact values.",
                       paste(manual, collapse = ", "))
  row$reproduces_app_cutoffs <- all(same)
  row
}

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

# One documentation row: all parameters, the final cutoffs and their source, and the deviations
parameters_used_row <- function(env, p, th, th_auto, sources, modality) {
  row <- params_row(env, p)
  row$data_type <- MODALITIES[[modality]]$label
  row$qc_script <- env$.script_path
  row$qc_script_modified <- format(env$.script_mtime, "%Y-%m-%d %H:%M")
  for (nm in names(sources)) {
    row[[paste0("cutoff_", nm)]] <- th[[nm]]
    row[[paste0("cutoff_", nm, "_source")]] <- if (sources[[nm]] == "manual") "manual" else "MAD rule"
  }
  dev <- param_deviations(env, p, th, th_auto, sources)
  row$n_deviations_from_default <- length(dev)
  row$deviations_from_default <- if (length(dev)) paste(dev, collapse = "; ") else "none (all script defaults)"
  row$saved_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  row
}

export_sample <- function(modality, env, d, meta, th, p, out_dir, th_auto, sources, extra = NULL,
                          progress = function(...) NULL, sample_meta = list()) {
  mod <- MODALITIES[[modality]]
  out_dir <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
  ds_dir  <- file.path(out_dir, p$dataset)
  dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(ds_dir)) stop("Could not create the output folder: ", ds_dir)

  # The script's log_msg() appends to this file while we save
  log_file <- file.path(out_dir, "qc_log.txt")
  env$.qc_log$file <- log_file
  on.exit(env$.qc_log$file <- NULL, add = TRUE)
  cut_txt <- paste(sprintf("%s=%s%s", names(sources), vapply(names(sources), function(nm) format(th[[nm]], digits = 6), ""),
                           ifelse(sources == "manual", sprintf(" (manual; MAD %s)", vapply(names(sources), function(nm) format(th_auto[[nm]], digits = 6), "")), "")),
                   collapse = ", ")
  env$log_msg("[%s] QC app (%s, script %s): cutoffs %s", p$sample_id, mod$label, basename(env$.script_path), cut_txt)

  # Thresholds table: script columns + what the app adds for documentation
  th_out <- th
  for (nm in names(sources)) th_out[[paste0("mad_", nm)]] <- th_auto[[nm]]
  th_out$manual_cutoffs <- paste(names(sources)[sources == "manual"], collapse = ";")
  th_out$saved_by <- "qc_explorer app"
  th_out$saved_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  progress("Saving the Seurat object")
  r <- result_object(modality, d, meta, th, p, extra)
  env$save_sample(r, ds_dir)
  rds <- file.path(ds_dir, paste0(p$sample_id, "_qc.rds"))
  sample_meta <- Filter(function(v) length(v) == 1 && !is.na(v) && nzchar(as.character(v)), sample_meta)
  if (length(sample_meta)) {   # patient ID etc.: one value per sample, added to every cell
    obj <- readRDS(rds)
    for (nm in names(sample_meta)) obj@meta.data[[nm]] <- sample_meta[[nm]]
    saveRDS(obj, rds)
    env$log_msg("[%s] QC app: sample metadata %s", p$sample_id,
                paste(sprintf("%s=%s", names(sample_meta), unlist(sample_meta)), collapse = ", "))
    for (nm in names(sample_meta)) meta[[nm]] <- sample_meta[[nm]]
  }

  progress("Writing tables")
  unit <- if (modality == "visium") "spot" else "cell"
  meta_file <- file.path(ds_dir, sprintf("%s_%s_qc_metadata.csv", p$sample_id, unit))
  write.csv(meta, meta_file, row.names = FALSE)
  upsert_csv(file.path(ds_dir, "qc_thresholds.csv"), th_out, "sample")
  sm <- env$summarise_sample(r); sm$manual_cutoffs <- th_out$manual_cutoffs
  sm$patient <- sample_meta$patient %||% NA_character_
  upsert_csv(file.path(out_dir, "all_samples_summary.csv"), sm, "sample")
  upsert_csv(file.path(out_dir, "samplesheet_resolved.csv"), params_row(env, p), "sample_id")
  used_row <- parameters_used_row(env, p, th, th_auto, sources, modality)
  used_row$patient <- sample_meta$patient %||% NA_character_
  for (nm in setdiff(names(sample_meta), "patient")) used_row[[paste0("meta_", nm)]] <- sample_meta[[nm]]
  used_row <- used_row[c("sample_id", "patient", setdiff(names(used_row), c("sample_id", "patient")))]
  upsert_csv(file.path(out_dir, "qc_parameters_used.csv"), used_row, "sample_id")
  if (used_row$n_deviations_from_default > 0)
    env$log_msg("[%s] QC app: deviations from defaults: %s", p$sample_id, used_row$deviations_from_default)
  sheet_row <- reproducible_sheet_row(modality, env, d, p, th, sources)
  if (identical(p$format, "h5ad")) sheet_row$notes <- paste(sheet_row$notes,
    "| Read from .h5ad by the QC app; the QC script can't read .h5ad, so save the counts as .rds (or 10x) to rerun it with the script.")
  upsert_csv(file.path(out_dir, "app_samplesheet_rows.csv"), sheet_row, "sample_id")

  progress("Drawing the QC plots")
  plot_dir <- file.path(ds_dir, "plots", p$sample_id)
  ok <- tryCatch({ env$make_plots(meta, th, plot_dir); TRUE },
                 error = function(e) { env$log_msg("[%s] WARNING: plotting failed (%s)", p$sample_id, conditionMessage(e)); FALSE })
  if (!ok) fallback_plots(modality, env, meta, th, th_auto, sources, plot_dir, p$sample_id)
  env$log_msg("[%s] QC app: saved %d of %d %s", p$sample_id, sum(meta$keep), nrow(meta), mod$unit)

  list(files = c(rds, meta_file, file.path(ds_dir, "qc_thresholds.csv"),
                 file.path(out_dir, c("all_samples_summary.csv", "samplesheet_resolved.csv", "qc_parameters_used.csv",
                                      "app_samplesheet_rows.csv", "qc_log.txt")),
                 list.files(plot_dir, full.names = TRUE)),
       reproducible = isTRUE(sheet_row$reproduces_app_cutoffs),
       deviations = used_row$deviations_from_default, n_deviations = used_row$n_deviations_from_default,
       n_kept = sum(meta$keep), n_in = nrow(meta), rds = rds)
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
