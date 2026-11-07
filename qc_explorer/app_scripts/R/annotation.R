# annotation.R — runs the annotation script (R_scripts/annotation/annotate_cells.R)
# on QC'd samples from inside the app. As with QC, the logic lives in the script; this file
# only turns the app's inputs into the script's PARAMS and draws the Annotation tab.
# Saving a sample writes its annotation with the script's own write_annotation_outputs()
# (labels into <sample>.rds and <sample>_cells.csv, plots into annotation/, settings into
# samplesheet.csv as annot_* columns), the same as the script does.

ANNOT_CHOICES <- list(label_level = c("main", "fine"), label_priority = c("SingleR", "markers"), marker_list_mode = c("add", "replace"))

# One input per annotation parameter, typed from the script's PARAMS
annotation_param_input <- function(nm, value, env, help) {
  id <- paste0("a_", nm)
  d <- env$PARAMS[[nm]]
  el <- if (nm == "references") {
    selectizeInput(id, nm, choices = names(env$REFERENCES), selected = value, multiple = TRUE,
                   options = list(placeholder = "none"))
  } else if (nm == "label_methods") {
    selectizeInput(id, nm, choices = env$LABEL_METHODS %||% c("SingleR", "markers"), selected = value, multiple = TRUE,
                   options = list(placeholder = "choose at least one"))
  } else if (nm == "malignant_method") {
    # one choice per method in the script's MALIGNANT_METHODS (+ its combinations), plus "none"
    ch <- env$malignant_choices()
    if (identical(value, "both")) value <- "copykat+infercnv"
    tagList(selectInput(id, nm, ch, selected = if (value %in% ch) value else "none"),
            conditionalPanel("input.a_malignant_method && input.a_malignant_method.indexOf('copykat') > -1",
                             div(class = "alert alert-warning py-2 px-2 mb-2 small", icon("triangle-exclamation"), " ",
                                 env$COPYKAT_WARNING %||% "copyKAT is slow and memory-hungry.")))
  } else if (!is.null(ANNOT_CHOICES[[nm]])) {
    selectInput(id, nm, ANNOT_CHOICES[[nm]], selected = value)
  } else if (is.logical(d) && length(d) == 1 && !is.na(d)) {
    checkboxInput(id, nm, value = isTRUE(value))
  } else if (is.numeric(d)) {
    numericInput(id, nm, value = if (is.null(value) || is.na(value)) NA else value)
  } else {
    textInput(id, nm, value = if (is.null(value) || length(value) == 0 || is.na(value)) "" else as.character(value))
  }
  div(class = "param-input", el, if (!is.null(help)) div(class = "form-text param-help", help))
}

# How the Annotation panel is laid out. Each section shows its `main` settings and folds the
# rarely changed ones into its own "Advanced" part. `when` is a JavaScript condition on the
# inputs: the section (or advanced part) only shows when it holds, e.g. only for the labelling
# or malignant-cell methods that are selected. Settings the script adds later, not listed
# here, go into "Other settings".
ANNOT_SECTIONS <- list(
  list(id = "methods", title = "Labelling methods", main = c("label_methods", "label_priority")),
  list(id = "singler", title = "SingleR", when = "input.a_label_methods && input.a_label_methods.indexOf('SingleR') > -1",
       main = c("references", "label_level"),
       advanced = list(list(title = "Custom reference",
                            params = c("custom_reference", "custom_label_col", "primary", "custom_max_per_type", "custom_min_per_type")))),
  list(id = "markers", title = "Marker genes", when = "input.a_label_methods && input.a_label_methods.indexOf('markers') > -1",
       note = "Uses the built-in marker list for the major bladder tumour cell types; add your own under Advanced.",
       main = character(0),
       advanced = list(list(title = "Custom markers & cutoffs",
                            params = c("custom_markers", "marker_list_mode", "marker_min_score", "marker_min_delta")))),
  list(id = "clustering", title = "Clustering", main = c("cluster_level", "cluster_resolution"),
       advanced = list(list(params = c("n_pcs", "n_hvg")))),
  list(id = "malignant", title = "Malignant cells (scRNA-seq)", main = "malignant_method",
       advanced = list(
         list(title = "Normal reference & cell limit", when = "input.a_malignant_method && input.a_malignant_method != 'none'",
              params = c("cnv_reference_types", "cnv_min_reference", "cnv_max_cells")),
         list(title = "inferCNV", when = "input.a_malignant_method && input.a_malignant_method.indexOf('infercnv') > -1",
              params = c("infercnv_gene_order", "cnv_score_quantile", "cnv_cor_cutoff")))),
  list(id = "other", title = "Other settings", main = character(0),
       advanced = list(list(params = c("assay", "min_shared_genes", "n_cores", "seed", "reuse_previous"))))
)
# Inside a section: settings that only make sense with another one on
ANNOT_PARAM_WHEN <- list(label_priority = "input.a_label_methods && input.a_label_methods.length > 1",
                      cluster_resolution = "input.a_cluster_level")

annot_settings_ui <- function(env) {
  sections <- ANNOT_SECTIONS
  listed <- unlist(lapply(sections, function(x) c(x$main, unlist(lapply(x$advanced, `[[`, "params")))))
  extra <- setdiff(names(env$PARAMS), listed)
  if (length(extra)) {
    k <- which(vapply(sections, `[[`, "", "id") == "other")
    sections[[k]]$advanced[[1]]$params <- c(sections[[k]]$advanced[[1]]$params, extra)
  }
  inputs <- function(nms) lapply(intersect(nms, names(env$PARAMS)), function(nm) {
    el <- annotation_param_input(nm, env$PARAMS[[nm]], env, env$.param_help[[nm]])
    if (!is.null(ANNOT_PARAM_WHEN[[nm]])) conditionalPanel(ANNOT_PARAM_WHEN[[nm]], el) else el
  })
  lapply(sections, function(sec) {
    adv <- Filter(function(a) length(intersect(a$params, names(env$PARAMS))), sec$advanced %||% list())
    if (!length(intersect(sec$main, names(env$PARAMS))) && !length(adv) && is.null(sec$note)) return(NULL)
    body_id <- paste0("annot-adv-", sec$id)
    adv_ui <- if (length(adv)) tagList(
      tags$button(type = "button", class = "btn btn-link annot-adv-toggle collapsed text-decoration-none ps-0",
                  `data-bs-toggle` = "collapse", `data-bs-target` = paste0("#", body_id), `aria-expanded` = "false",
                  icon("sliders"), " Advanced ", textOutput(paste0("annot_adv_n_", sec$id), inline = TRUE),
                  span(class = "qc-chevron", icon("chevron-down"))),
      div(id = body_id, class = "collapse",
          div(class = "annot-adv-body",
              lapply(adv, function(a) {
                part <- tagList(if (!is.null(a$title) && length(adv) > 1) h6(a$title, class = "mt-2 mb-2"), inputs(a$params))
                if (!is.null(a$when)) conditionalPanel(a$when, part) else part
              }))))
    body <- div(class = "annot-section",
                h6(sec$title, class = "mb-2"),
                if (!is.null(sec$note)) div(class = "form-text mb-2", sec$note),
                inputs(sec$main), adv_ui)
    if (!is.null(sec$when)) conditionalPanel(sec$when, body) else body
  })
}

# Settings in a section's Advanced part that differ from the script's defaults (so a hidden
# change still shows on the closed section)
annot_adv_changed <- function(env, params, sec_id) {
  sec <- Filter(function(x) identical(x$id, sec_id), ANNOT_SECTIONS)[[1]]
  nms <- unlist(lapply(sec$advanced, `[[`, "params"))
  if (sec_id == "other") nms <- c(nms, setdiff(names(env$PARAMS), unlist(lapply(ANNOT_SECTIONS, function(x) c(x$main, unlist(lapply(x$advanced, `[[`, "params")))))))
  fmt <- function(v) if (is.null(v) || !length(v) || all(is.na(v))) "" else paste(v, collapse = ";")
  nms <- intersect(nms, names(env$PARAMS))
  nms[vapply(nms, function(nm) !identical(fmt(params[[nm]]), fmt(env$PARAMS[[nm]])), TRUE)]
}

# The script's PARAMS with the app's inputs applied
collect_annotation_params <- function(input, env) {
  p <- env$PARAMS
  drawn <- !is.null(input$a_label_level)                  # the settings UI exists
  for (nm in names(p)) {
    v <- input[[paste0("a_", nm)]]
    # an emptied multi-select is NULL, which here means "no built-in references"
    if (nm %in% c("references", "label_methods")) { if (drawn) p[[nm]] <- as.character(v %||% character(0)); next }
    if (is.null(v)) next                                  # UI not drawn yet: keep the default
    d <- env$PARAMS[[nm]]
    if (is.numeric(d)) { p[[nm]] <- if (is.na(v)) NA_real_ else as.numeric(v); next }
    if (is.logical(d) && !is.na(d)) { p[[nm]] <- isTRUE(v); next }
    v <- trimws(as.character(v))
    p[[nm]] <- if (!nzchar(v)) NA else gsub('^["\']|["\']$', "", v)
  }
  p
}

# References take a while to load (download / normalize), so keep them between samples
.ref_cache <- new.env()
get_references <- function(env, params) {
  key <- paste(deparse(list(params[c("label_methods", "label_priority", "references", "label_level", "custom_reference", "custom_label_col",
                                     "custom_max_per_type", "custom_min_per_type", "primary", "seed")],
                            if (!is.na(params$custom_reference)) file.mtime(params$custom_reference),
                            env$.script_mtime)), collapse = "")
  if (!identical(.ref_cache$key, key)) {
    .ref_cache$refs <- env$prepare_references(params)
    .ref_cache$key <- key
  }
  .ref_cache$refs
}

# What a preview depends on: the cells kept and every annotation parameter
annotation_key <- function(cells, params) list(cells = sort(cells), params = params)

# Per-sample working folder (for as long as this R session runs, shared by every page). The
# labels checkpoint, each malignant-cell method's cached result and the copyKAT / inferCNV
# working files live here, for previews and saving alike: changing only the malignant-cell
# settings reuses the labels, switching back to a method reuses its result, and saving after a
# preview (or after the page reconnected) reuses what the preview computed. The tools' heatmaps
# are copied into the saved plots/ folder.
preview_dir <- function(id) file.path(normalizePath(tempdir(), winslash = "/"), "qcapp_preview", id)
labels_checkpoint <- function(id) file.path(preview_dir(id), paste0(id, "_labels_checkpoint.rds"))

# Annotate a QC'd object saved at rds_path. Returns the script's result list.
# checkpoint: where the script keeps (and looks for) this sample's labels; see annotate_one_sample()
annotate_saved_sample <- function(env, rds_path, sample_id, dataset, sample_meta, params, work_dir,
                                  checkpoint = labels_checkpoint(sample_id)) {
  env$check_packages(params)
  refs <- get_references(env, params)
  row <- data.frame(sample_id = sample_id, dataset = dataset, path = rds_path,
                    patient = sample_meta$patient %||% NA_character_, stringsAsFactors = FALSE)
  r <- env$annotate_one_sample(row, refs, params, work_dir = work_dir, checkpoint = checkpoint)
  for (nm in setdiff(names(sample_meta), "patient")) r$obj@meta.data[[nm]] <- sample_meta[[nm]]
  r$refs <- refs
  r
}

# Short text summary for the app
# The final label column (older annotated objects only have singler_label)
final_label_of <- function(obj) if ("cell_label" %in% colnames(obj@meta.data)) obj$cell_label else obj$singler_label

annotation_summary <- function(obj) {
  fl <- final_label_of(obj)
  lab <- table(ifelse(is.na(fl), "(unlabelled)", fl))
  lab <- sort(lab, decreasing = TRUE)
  out <- paste(sprintf("%s %d", names(lab), lab), collapse = ", ")
  if ("cnv_status" %in% colnames(obj[[]])) {
    cs <- table(obj$cnv_status)
    out <- paste0(out, " | CNV: ", paste(sprintf("%s %d", names(cs), cs), collapse = ", "))
  }
  out
}

## ---- plots for the Annotation tab ------------------------------------------------
# One UMAP coloured by a metadata column. No text on the plot: the legend lists every
# label with its cell count (largest first), so many small labels can't cover the cells.
umap_panel <- function(obj, col, title, colours = NULL) {
  emb <- Embeddings(obj, "umap")
  g <- as.character(obj[[col, drop = TRUE]]); g[is.na(g)] <- "(pruned)"
  n <- sort(table(g), decreasing = TRUE)
  lev <- c(setdiff(names(n), "(pruned)"), intersect("(pruned)", names(n)))
  labs_n <- sprintf("%s (%s)", lev, format(as.integer(n[lev]), big.mark = ",", trim = TRUE))
  set.seed(1)
  ord <- sample(nrow(emb))                                   # random draw order: no type hides another
  df <- data.frame(x = emb[ord, 1], y = emb[ord, 2], g = factor(g[ord], levels = lev, labels = labs_n))
  pt <- if (nrow(df) > 30000) 0.15 else if (nrow(df) > 8000) 0.3 else 0.6
  pal <- if (!is.null(colours)) {
    unname(ifelse(lev %in% names(colours), colours[lev], "grey60"))
  } else {
    p <- scales::hue_pal()(sum(lev != "(pruned)")); c(p, if ("(pruned)" %in% lev) "grey80")
  }
  ggplot(df, aes(x, y, colour = g)) +
    geom_point(size = pt, shape = 16, stroke = 0) +
    scale_colour_manual(values = setNames(pal, labs_n)) +
    coord_equal() +
    guides(colour = guide_legend(override.aes = list(size = 2.5), ncol = if (length(lev) > 12) 2 else 1)) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2", colour = NULL) +
    theme_classic(base_size = 10) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          plot.title = element_text(face = "bold", size = 11), plot.title.position = "plot",
          legend.text = element_text(size = 8), legend.key.height = unit(9, "pt"), legend.key.width = unit(9, "pt"),
          legend.box.spacing = unit(4, "pt"))
}

# The Annotation tab's panels: final label, its broad lineage, cluster-level label(s), SingleR vs
# markers agreement, malignant calls. Used by the plot and by the plot height in app.R.
annotation_panels <- function(r) {
  obj <- r$obj; md <- obj[[]]
  lab_col <- if ("cell_label" %in% colnames(md)) "cell_label" else "singler_label"
  methods <- if ("cell_label_source" %in% colnames(md)) sort(unique(na.omit(md$cell_label_source))) else "SingleR"
  primary <- if ("singler_primary_reference" %in% colnames(md)) md$singler_primary_reference[1] else NA
  if (!"umap" %in% Reductions(obj)) {
    fl <- final_label_of(obj)
    d <- as.data.frame(sort(table(label = ifelse(is.na(fl), "(unlabelled)", fl)), decreasing = TRUE))
    return(list(ggplot(d, aes(reorder(label, Freq), Freq)) + geom_col(fill = "#5B7DB1") + coord_flip() +
                  labs(title = "Final labels", x = NULL, y = "cells") + theme_qc()))
  }
  pl <- list(umap_panel(obj, lab_col, if (lab_col == "cell_label") sprintf("Final label (%s)", paste(methods, collapse = " + "))
                                       else sprintf("Per-cell label (%s, pruned)", primary)))
  if ("cell_lineage" %in% colnames(md) && !identical(md$cell_lineage, md$cell_label))
    pl <- c(pl, list(umap_panel(obj, "cell_lineage", "Final label as broad lineage")))
  cl_cols <- c(SingleR = if (!is.na(primary)) paste0("singler_", primary, "_cluster_label"), markers = "marker_cluster_label")
  for (m in names(cl_cols)) if (cl_cols[[m]] %in% colnames(md))
    pl <- c(pl, list(umap_panel(obj, cl_cols[[m]], sprintf("Cluster-level label (%s)", m))))
  if ("label_agreement" %in% colnames(md))
    pl <- c(pl, list(umap_panel(obj, "label_agreement", "SingleR vs markers",
                                c(agree = "#4C9BE8", disagree = "#C0392B", "one method only" = "#E69F00", "(pruned)" = "grey85"))))
  if ("cnv_status" %in% colnames(md))
    pl <- c(pl, list(umap_panel(obj, "cnv_status", "Malignant calls (inferred CNVs)", APP_CNV_COLOURS)))
  if ("cnv_status" %in% colnames(md) && any(!is.na(obj$cnv_score)))
    pl <- c(pl, list(ggplot(md, aes(cnv_score, cnv_cor, colour = cnv_status)) + geom_point(size = 0.5, alpha = 0.7, shape = 16) +
                       scale_colour_manual(values = APP_CNV_COLOURS) + theme_qc() +
                       guides(colour = guide_legend(override.aes = list(size = 2.5))) +
                       theme(plot.title.position = "plot", plot.title = element_text(size = 11),
                             legend.text = element_text(size = 8)) +
                       labs(title = "inferCNV score vs tumour-profile correlation", x = "CNV score", y = "correlation", colour = NULL)))
  pl
}
plot_annotation_result <- function(r) wrap_plots(annotation_panels(r), ncol = 2)

APP_CNV_COLOURS <- c(malignant = "#C0392B", normal = "#4C9BE8", "normal (reference)" = "grey70",
                     uncertain = "#E69F00", unresolved = "grey40", "not tested" = "grey90")

annotation_table <- function(r) {
  md <- r$obj[[]]
  fl <- final_label_of(r$obj)
  lab <- ifelse(is.na(fl), "(unlabelled)", fl)
  if ("cnv_status" %in% names(md)) {
    tb <- as.data.frame.matrix(table(lab, md$cnv_status))
    data.frame(Label = rownames(tb), Cells = as.integer(rowSums(tb)), tb, check.names = FALSE)[order(-rowSums(tb)), ]
  } else {
    tb <- sort(table(lab), decreasing = TRUE)
    data.frame(Label = names(tb), Cells = as.integer(tb), `%` = sprintf("%.1f", 100 * as.integer(tb) / nrow(md)), check.names = FALSE)
  }
}

# The script's check reports a package as "missing" whenever it can't be loaded. Say why:
# not installed in this R's libraries, or installed but failing to load (and the error).
explain_missing <- function(msg) {
  if (!grepl("Missing packages:", msg)) return(msg)
  pkgs <- trimws(strsplit(sub("(?s)^.*Missing packages: ([^\n]+).*$", "\\1", msg, perl = TRUE), ",")[[1]])
  why <- vapply(pkgs, function(p) {
    loc <- tryCatch(find.package(p, quiet = TRUE), error = function(e) character(0))
    if (!length(loc)) return(sprintf("- %s: not installed in the libraries this R uses:\n    %s", p, paste(.libPaths(), collapse = "\n    ")))
    v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) "?")
    err <- tryCatch({ loadNamespace(p); NULL }, error = function(e) conditionMessage(e))
    if (is.null(err)) return(sprintf("- %s %s: installed and loads now; restart the app.", p, v))
    sprintf("- %s %s is installed (%s) but fails to load:\n    %s", p, v, loc, gsub("\n", "\n    ", err))
  }, "")
  paste0(msg, "\n\nWhy:\n", paste(why, collapse = "\n"), "\n\nThe app is running ", R.version.string, ".")
}
