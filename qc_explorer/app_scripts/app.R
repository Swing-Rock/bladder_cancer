# QC explorer — interactive front end for the lab's QC scripts
# (R_scripts/: scRNA/scrna_qc.R, visium/visium_qc.R, xenium/xenium_qc.R, annotation/annotate_cells.R).
#
# Run it from RStudio (open this file, click "Run App"), or open ../run_app.R and click Source, or:
#   Rscript run_app.R   (from the qc_explorer folder)
#
# Files in R/ are loaded automatically by Shiny before this file.

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
  library(Seurat)
})
HAS_SHINYFILES <- requireNamespace("shinyFiles", quietly = TRUE)
`%||%` <- function(a, b) if (is.null(a)) b else a
is_manual <- function(src, nm) isTRUE(unname(src[nm]) == "manual")
MAX_BATCH_DEFAULT <- 8

clean_path <- function(x) {
  x <- trimws(x %||% "")
  gsub('^["\']|["\']$', "", x)   # Windows "Copy as path" adds quotes
}

# shinyFiles buttons each carry its script; included more than once, every click opens two dialogs.
# So the buttons are used without it, and the page loads it once (shinyfiles_head below).
drop_head <- function(x) x[!vapply(x, function(t) inherits(t, "shiny.tag") && identical(t$name, "head"), TRUE)]
browse_button <- function(id, label, type = c("dir", "file"), title = label) {
  if (!HAS_SHINYFILES) return(NULL)
  type <- match.arg(type)
  b <- if (type == "dir") shinyFiles::shinyDirButton(id, label, title, class = "btn-sm btn-outline-secondary")
       else shinyFiles::shinyFilesButton(id, label, title, multiple = FALSE, class = "btn-sm btn-outline-secondary")
  drop_head(b)
}
shinyfiles_head <- function() {
  if (!HAS_SHINYFILES) return(NULL)
  b <- shinyFiles::shinyDirButton("qc_sf_unused", "", "")
  b[vapply(b, function(t) inherits(t, "shiny.tag") && identical(t$name, "head"), TRUE)]
}

# One input per sample-sheet parameter, typed from the script's DEFAULTS
param_input <- function(nm, value, env, help) {
  id <- paste0("p_", nm)
  num <- exists("NUMERIC_PARAMS", envir = env) && nm %in% env$NUMERIC_PARAMS
  lgl <- exists("LOGICAL_PARAMS", envir = env) && nm %in% env$LOGICAL_PARAMS
  el <- if (nm == "format") {
    selectInput(id, nm, c("auto", env$KNOWN_FORMATS, "h5ad"), selected = value %||% "auto")
  } else if (nm == "is_filtered") {
    selectInput(id, nm, c("auto", "TRUE", "FALSE"), selected = as.character(value %||% "auto"))
  } else if (nm == "gene_column") {
    selectInput(id, nm, c("2 = gene symbol" = "2", "1 = Ensembl ID" = "1"), selected = as.character(value %||% 2))
  } else if (lgl) {
    checkboxInput(id, nm, value = isTRUE(as.logical(value)))
  } else if (num) {
    numericInput(id, nm, value = if (is.null(value) || is.na(value)) NA else as.numeric(value))
  } else {
    textInput(id, nm, value = if (is.null(value) || is.na(value)) "" else as.character(value))
  }
  div(class = "param-input", el, if (!is.null(help)) div(class = "form-text param-help", help))
}

# A button that shows a spinner and a message the moment it's clicked, until the app is idle again
busy_button <- function(id, label, busy_text, ...) {
  b <- actionButton(id, label, ...)
  b <- tagAppendAttributes(b, `data-busy-text` = busy_text)
  b$attribs$class <- paste(b$attribs$class, "qc-busy-btn")
  b
}

html_table <- function(df) {
  tags$table(class = "table table-sm table-striped small",
             tags$thead(tags$tr(lapply(names(df), tags$th))),
             tags$tbody(lapply(seq_len(nrow(df)), function(i) tags$tr(lapply(df[i, ], function(v) tags$td(as.character(v)))))))
}

# A full-window dialog. Not a Bootstrap modal, so the file browser's own dialog can open on top of it.
overlay <- function(id, title, ..., footer = NULL, show = FALSE) {
  div(id = id, class = paste("qc-overlay", if (show) "show"),
      div(class = "qc-dialog",
          div(class = "qc-dialog-header", h5(title),
              tags$button(type = "button", class = "btn-close", `aria-label` = "Close", `data-qc-close` = id)),
          div(class = "qc-dialog-body", ...),
          if (!is.null(footer)) div(class = "qc-dialog-footer", footer)))
}

STATUS_STYLE <- list(
  pending   = list(icon = "circle",               cls = "text-secondary", label = "not loaded"),
  loaded    = list(icon = "circle-dot",           cls = "text-primary",   label = "in memory, not saved"),
  saved     = list(icon = "circle-check",         cls = "text-success",   label = "saved"),
  on_disk   = list(icon = "circle-check",         cls = "text-success opacity-50", label = "already in output folder"),
  unsaved   = list(icon = "circle-exclamation",   cls = "text-warning",   label = "unloaded without saving"),
  failed    = list(icon = "circle-xmark",         cls = "text-danger",    label = "failed"))

## =============================================================================
## UI
## =============================================================================
ui <- page_sidebar(
  title = div(class = "d-flex align-items-center gap-2", span("QC explorer", class = "fw-semibold"),
              span(class = "text-muted small", "· scRNA-seq · Visium · Xenium")),
  theme = bs_theme(version = 5, primary = "#2C5F8A", "font-size-base" = "0.92rem"),
  fillable = FALSE,
  shinyfiles_head(),
  tags$head(tags$style(HTML("
    .param-help { font-size: 0.75rem; margin-top: -0.6rem; margin-bottom: 0.6rem; }
    .param-input .form-group { margin-bottom: 0.5rem; }
    .slider-readout { font-size: 0.8rem; margin-top: -0.6rem; margin-bottom: 0.9rem; color: var(--bs-secondary-color); }
    .slider-readout .manual { color: #B8590B; font-weight: 600; }
    .qc-cut-boxes { display: flex; gap: .5rem; margin-top: -0.5rem; margin-bottom: 0.35rem; }
    .qc-cut-boxes label { display: flex; align-items: center; gap: .3rem; font-size: .75rem; color: var(--bs-secondary-color); margin: 0; flex: 1; }
    .qc-cut-boxes input { width: 100%; min-width: 0; padding: .1rem .35rem; font-size: .8rem; height: auto; }
    .annot-section { border-top: 1px solid var(--bs-border-color); padding-top: .5rem; margin-top: .6rem; }
    .annot-adv-toggle { font-size: .78rem; padding: .1rem .4rem; }
    .annot-adv-toggle .qc-chevron { display: inline-block; transition: transform .2s; }
    .annot-adv-toggle:not(.collapsed) .qc-chevron { transform: rotate(180deg); }
    .annot-adv-body { border-left: 2px solid var(--bs-border-color); padding-left: .6rem; margin: .4rem 0 .2rem; }
    .annot-adv-toggle .shiny-text-output { color: #B8590B; font-weight: 600; }
    .status-banner { font-size: 0.85rem; }
    .sidebar h6 { text-transform: uppercase; letter-spacing: .04em; font-size: .72rem; color: var(--bs-secondary-color); }
    .bslib-value-box .value-box-value { font-size: 1.6rem; }
    .sample-col { width: 250px; flex-shrink: 0; position: sticky; top: 0.5rem; }
    .sample-list { max-height: 62vh; overflow-y: auto; }
    .sample-list .list-group-item { padding: .45rem .7rem; font-size: .85rem; }
    .sample-list .list-group-item.active .text-secondary, .sample-list .list-group-item.active .text-success,
    .sample-list .list-group-item.active .text-primary, .sample-list .list-group-item.active .small { color: #fff !important; }
    .sample-list .list-group-item.disabled { opacity: .75; }
    .sample-meta { font-size: .72rem; }
    #qc-adv-toggle .qc-chevron { display: inline-block; transition: transform .2s; }
    #qc-adv-toggle:not(.collapsed) .qc-chevron { transform: rotate(180deg); }
    #qc-busy-banner { display: none; position: fixed; top: 12px; left: 50%; transform: translateX(-50%); z-index: 2000;
                      background: var(--bs-primary); color: #fff; padding: .55rem 1.1rem; border-radius: 2rem;
                      box-shadow: 0 4px 14px rgba(0,0,0,.18); font-size: .9rem; }
    body.qc-busy #qc-busy-banner, body.qc-queue #qc-busy-banner { display: flex; align-items: center; gap: .6rem; }
    #qc-busy-banner .qc-queue-only { display: none; }
    body.qc-queue #qc-busy-banner .qc-queue-only { display: inline; }
    #qc-busy-banner .qc-queue-only a { color: #fff; text-decoration: underline; }
    .sample-list .list-group-item.qc-job-queued { opacity: .5; }
    .sample-list .list-group-item.qc-job-running { background: var(--bs-primary-bg-subtle); }
    .sample-list .list-group-item.active.qc-job-running { background: var(--bs-primary); }
    body.qc-busy .qc-busy-btn { pointer-events: none; }
    .qc-overlay { position: fixed; inset: 0; z-index: 1040; display: none; background: rgba(15, 23, 42, .45);
                  align-items: flex-start; justify-content: center; overflow-y: auto; padding: 5vh 1rem; }
    .qc-overlay.show { display: flex; }
    body.qc-overlay-open { overflow: hidden; }
    .qc-dialog { background: var(--bs-body-bg); border-radius: .8rem; box-shadow: 0 24px 60px rgba(0,0,0,.28); width: 100%; max-width: 620px; }
    .qc-dialog-header { display: flex; align-items: center; justify-content: space-between; padding: 1rem 1.25rem .5rem; }
    .qc-dialog-header h5 { margin: 0; font-weight: 600; }
    .qc-dialog-body { padding: .5rem 1.25rem 1rem; }
    .qc-dialog-footer { display: flex; gap: .5rem; justify-content: flex-end; padding: .75rem 1.25rem 1rem; border-top: 1px solid var(--bs-border-color); }
    .qc-collapse-toggle .qc-chevron { display: inline-block; transition: transform .2s; }
    .qc-collapse-toggle:not(.collapsed) .qc-chevron { transform: rotate(180deg); }
    .qc-data-box { font-size: .85rem; }
  ")),
  tags$script(HTML("
    // Number boxes under the QC sliders: a typed value is sent when you press Enter or leave the box
    $(document).on('change', '.qc-cut-box', function() {
      Shiny.setInputValue('cut_box', {slider: $(this).data('slider'), side: $(this).data('side'), value: this.value,
                                      nonce: Math.random()}, {priority: 'event'});
    });
    $(document).on('keydown', '.qc-cut-box', function(e) { if (e.key === 'Enter') { $(this).blur(); } });
    // ...and they follow the sliders (except the box being typed in)
    Shiny.addCustomMessageHandler('qc-cut-boxes', function(vals) {
      [].concat(vals).forEach(function(v) {
        var el = $('.qc-cut-box[data-slider=\"' + v.slider + '\"][data-side=\"' + v.side + '\"]');
        if (el.length && !el.is(':focus')) { el.val(v.value); }
      });
    });
    // Log-scale sliders hold log10 values; show the real numbers on the handles and ends
    Shiny.addCustomMessageHandler('qc-log-sliders', function(ids) {
      function apply(id, tries) {
        var s = $('#' + id).data('ionRangeSlider');
        if (!s) { if (tries < 50) setTimeout(function() { apply(id, tries + 1); }, 100); return; }
        s.update({prettify_enabled: true, prettify: function(n) {
          if (n >= s.options.max - 0.005) return 'no max';
          return Math.round(Math.pow(10, n)).toLocaleString();
        }});
      }
      [].concat(ids).forEach(function(id) { apply(id, 0); });
    });
    // Busy feedback: the clicked button turns into a spinner and a banner appears straight away
    // (before the server has even started), and both go when the server is idle again.
    var qcBusy = {on: false, sawBusy: false, btn: null, html: null, timer: null};
    function qcBusyEnd() {
      if (!qcBusy.on) return;
      if (qcBusy.btn) { qcBusy.btn.html(qcBusy.html).prop('disabled', false); }
      $('body').removeClass('qc-busy'); clearTimeout(qcBusy.timer);
      qcBusy = {on: false, sawBusy: false, btn: null, html: null, timer: null};
    }
    $(document).on('click', '.qc-busy-btn', function() {
      var b = $(this), text = b.data('busy-text') || 'Working…';
      qcBusyEnd();
      qcBusy = {on: true, sawBusy: false, btn: b, html: b.html(), timer: null};
      setTimeout(function() {   // after Shiny has read the click
        if (!qcBusy.on) return;
        b.empty().append($('<span>').addClass('spinner-border spinner-border-sm me-2'), document.createTextNode(text)).prop('disabled', true);
      }, 0);
      $('#qc-busy-text').text(text);
      $('body').addClass('qc-busy');
      qcBusy.timer = setTimeout(function() { if (!qcBusy.sawBusy) qcBusyEnd(); }, 5000);   // nothing happened: give the button back
    });
    $(document).on('shiny:busy', function() { if (qcBusy.on) qcBusy.sawBusy = true; });
    $(document).on('shiny:idle', function() { if (qcBusy.on && qcBusy.sawBusy) qcBusyEnd(); });
    $(document).on('shown.bs.modal', function() { qcBusyEnd(); });   // a question for the user: stop the spinner
    // Dataset-level runs: the server keeps the banner up (with progress) until the queue is done
    Shiny.addCustomMessageHandler('qc-queue', function(m) {
      if (m.on) { $('body').addClass('qc-queue'); $('#qc-busy-text').text(m.text); }
      else { $('body').removeClass('qc-queue'); }
    });

    // Load dialog: shown when the app opens, reopened from the sidebar
    function qcOverlay(id, show) {
      $('#' + id).toggleClass('show', !!show);
      $('body').toggleClass('qc-overlay-open', $('.qc-overlay.show').length > 0);
      if (show) qcBusyEnd();
    }
    Shiny.addCustomMessageHandler('qc-overlay', function(m) { qcOverlay(m.id, m.show); });
    $(document).on('click', '[data-qc-close]', function() { qcOverlay($(this).data('qc-close'), false); });
    $(document).on('click', '[data-qc-open]', function() { qcOverlay($(this).data('qc-open'), true); });
    $(document).on('keydown', function(e) {
      if (e.key === 'Escape' && !$('.modal.show').length) $('.qc-overlay.show').each(function() { qcOverlay(this.id, false); });
    });

    $(document).on('click', '.sample-pick', function() {
      if ($(this).hasClass('disabled')) return;
      Shiny.setInputValue('pick_sample', $(this).data('sample'), {priority: 'event'});
    });
  "))),
  sidebar = sidebar(
    width = 390,
    div(class = "card card-body p-2 qc-data-box",
        uiOutput("data_box"),
        tags$button(type = "button", class = "btn btn-sm btn-outline-primary mt-2", `data-qc-open` = "qc-load",
                    icon("folder-open"), " Load data…")),
    accordion(
      id = "side", open = c("QC cutoffs"), multiple = TRUE,
      accordion_panel(
        "Patient & sample info", icon = icon("user"),
        uiOutput("patient_ui"),
        conditionalPanel("input.input_mode == 'dataset'",
          div(class = "form-text mt-2 mb-1", "Or import a table for the whole dataset: a .csv with ",
              code("sample_id"), " and ", code("patient_id"), " columns. Any other columns (tissue, stage, treatment, …) are added too."),
          textInput("meta_path", NULL, placeholder = "path to sample_metadata.csv"),
          browse_button("browse_meta", "Browse…", "file", "Choose a sample metadata table"))
      ),
      accordion_panel(
        "QC cutoffs", icon = icon("sliders"),
        uiOutput("cutoff_heading"),
        uiOutput("sliders_ui"),
        uiOutput("reset_auto_ui"),
        # Every QC parameter, folded away under the sliders
        tags$button(id = "qc-adv-toggle", type = "button", class = "btn btn-sm btn-outline-secondary w-100 mt-2 collapsed",
                    `data-bs-toggle` = "collapse", `data-bs-target` = "#qc-advanced", `aria-expanded` = "false",
                    icon("gear"), " Advanced settings ", span(class = "qc-chevron", icon("chevron-down"))),
        div(id = "qc-advanced", class = "collapse",
          div(class = "pt-2",
            div(class = "form-text mb-2",
                "The QC script's filtering parameters, same defaults as its sample sheet, applied to every sample. ",
                "Cutoff rules update the sliders and plots live."),
            uiOutput("advanced_qc_ui")))
      ),
      accordion_panel(
        "Annotation", icon = icon("tags"),
        uiOutput("annotation_ui")
      ),
      accordion_panel(
        "Advanced settings", icon = icon("gear"),
        div(class = "form-text mb-2",
            "General settings, same defaults as the script's sample sheet. 'Reading the data' settings apply when you load again."),
        h6("Pipeline", class = "mt-2 mb-2"),
        div(class = "param-input",
            input_switch("auto_pipeline", "Run the next steps automatically after Apply QC", value = FALSE),
            div(class = "form-text param-help",
                "After Apply QC, runs scDblFinder + DecontX (scRNA-seq) and then the annotation preview, in that order, on the ",
                "cells that pass the MAD cutoffs. They run in the background; the banner at the top shows progress. With ",
                "dataset-level processing, each step runs on every sample in memory before the next one starts. Moving a slider ",
                "afterwards makes the results out of date: run the steps again before saving (saving doesn't run them). The annotation ",
                "uses the Annotation panel's settings, including the malignant-cell method (copyKAT can take an hour; choose ",
                "'none' there for a quick pass).")),
        conditionalPanel("input.input_mode == 'dataset'",
          h6("Dataset mode", class = "mt-2 mb-2"),
          div(class = "param-input",
              input_switch("dataset_level", "Dataset-level processing", value = FALSE),
              div(class = "form-text param-help",
                  "When on, 'Run scDblFinder + DecontX' and 'Annotate (preview)' run on every sample in memory, ",
                  "one after another. Waiting samples are greyed out with a spinner in the sample list; you can still click them ",
                  "to see their earlier results.")),
          div(class = "param-input",
              numericInput("max_batch", "Max samples in memory at once", value = MAX_BATCH_DEFAULT, min = 1, max = 50, step = 1),
              div(class = "form-text param-help",
                  sprintf("Samples are loaded in batches of up to this many; fewer if free RAM would drop below %g GB.", RAM_RESERVE_GB)))),
        uiOutput("advanced_ui")
      ),
      accordion_panel(
        "Save", icon = icon("floppy-disk"),
        textInput("out_dir", "Output folder", placeholder = "e.g. .../datasets/GSE.../qc_output"),
        browse_button("browse_out", "Choose folder…", "dir", "Choose the output folder"),
        div(class = "form-text mb-2", "Each sample goes to <output folder>/<dataset>/<sample>/, with samplesheet.csv, summary.csv and ",
            "log.txt at the top: the same layout as the scripts. Saving writes what is there: the QC result, plus ",
            "scDblFinder + DecontX and the annotation if they were run on the cells being saved. It runs nothing itself."),
        busy_button("save", "Save this sample", "Saving…", icon = icon("floppy-disk"), class = "btn-success w-100 mb-2"),
        conditionalPanel("input.input_mode == 'dataset'",
          busy_button("save_all", "Save all samples in memory", "Saving all samples in memory…", icon = icon("layer-group"), class = "btn-outline-success w-100"))
      )
    ),
    div(class = "form-text mt-2", uiOutput("script_info"))
  ),

  div(id = "qc-busy-banner", span(class = "spinner-border spinner-border-sm", role = "status"),
      span(id = "qc-busy-text", "Working…"),
      span(class = "qc-queue-only", actionLink("queue_stop_banner", "Stop after this one")),
      span(class = "qc-queue-only", "·"),
      span(class = "qc-queue-only", actionLink("queue_terminate_banner", "Terminate now"))),
  overlay("qc-load", "Load data", show = TRUE,
    radioButtons("modality", "Data type", inline = TRUE,
                 choices = c("scRNA-seq" = "scrna", "Visium" = "visium", "Xenium" = "xenium")),
    radioButtons("input_mode", "What to load", inline = TRUE,
                 choices = c("One sample" = "single", "Dataset (a folder, or a sample sheet)" = "dataset")),
    textInput("path", "Path", placeholder = "paste a folder or file path (dataset: a folder, or a sample sheet .csv)", width = "100%"),
    div(class = "d-flex gap-2 mb-2",
        browse_button("browse_in_dir", "Folder…", "dir", "Choose the folder"),
        browse_button("browse_in_file", "File…", "file", "Choose the sample file (or, for a dataset, its sample sheet .csv)")),
    div(class = "form-text mb-3", textOutput("path_hint", inline = TRUE)),
    layout_columns(col_widths = c(6, 6),
                   conditionalPanel("input.input_mode == 'single'", textInput("sample_id", "Sample ID")),
                   textInput("dataset", "Dataset", value = "dataset1")),
    conditionalPanel("input.input_mode == 'dataset'",
      checkboxInput("skip_done", "Skip samples already saved in the output folder", TRUE)),
    div(class = "form-text", "Loading reads the data and measures every cell or spot. Nothing is filtered until you click ",
        strong("Apply QC"), ", so data that has already been QC'd can be loaded and looked at as it is."),
    conditionalPanel("input.input_mode == 'single'",
      tags$button(type = "button", class = "btn btn-sm btn-outline-secondary w-100 mt-3 collapsed qc-collapse-toggle",
                  `data-bs-toggle` = "collapse", `data-bs-target` = "#qc-fill-sheet", `aria-expanded` = "false",
                  icon("table"), " Fill in from a sample sheet ", span(class = "qc-chevron", icon("chevron-down"))),
      div(id = "qc-fill-sheet", class = "collapse",
          div(class = "pt-2",
              div(class = "form-text mb-2", "One sample: copy its path and settings from an existing sample sheet. ",
                  "(To load every row of a sample sheet, choose 'Dataset' above and give the sheet as the path.)"),
              textInput("sheet_path", "Sample sheet (.csv)", placeholder = "path to samplesheet.csv"),
              browse_button("browse_sheet", "Browse…", "file", "Choose a sample sheet"),
              uiOutput("sheet_rows_ui")))),
    footer = tagList(tags$button(type = "button", class = "btn btn-outline-secondary", `data-qc-close` = "qc-load", "Close"),
                     uiOutput("load_button", inline = TRUE))),
  div(class = "d-flex gap-3 align-items-start",
      uiOutput("sample_column"),
      div(class = "flex-grow-1", style = "min-width: 0;",
          uiOutput("status"),
          uiOutput("value_boxes"),
          navset_card_tab(
            id = "tabs",
            nav_panel("Distributions", plotOutput("dist_plot", height = "auto")),
            nav_panel("Scatter",
                      plotOutput("scatter_plot", height = "460px"),
                      div(class = "form-text", textOutput("scatter_note", inline = TRUE))),
            nav_panel("Spatial",
                      div(class = "d-flex flex-wrap gap-3 align-items-end",
                          selectInput("sp_colour", "Colour by", choices = c("QC result" = "qc_reason"), width = "220px"),
                          conditionalPanel("input.modality == 'visium'", checkboxInput("sp_image", "Show tissue image", TRUE)),
                          checkboxInput("sp_removed", "Removed only", FALSE),
                          sliderInput("sp_size", "Point size", min = 0.1, max = 3, value = 1, step = 0.1, width = "180px", ticks = FALSE)),
                      plotOutput("spatial_plot", height = "640px"),
                      div(class = "form-text", textOutput("spatial_note", inline = TRUE))),
            nav_panel("Doublets & ambient",
                      div(class = "d-flex gap-2 align-items-center mb-2",
                          busy_button("run_extra", "Run scDblFinder + DecontX on the cells that pass", "Running scDblFinder + DecontX…", icon = icon("flask"), class = "btn-outline-primary"),
                          span(class = "form-text", "Runs the script's detect_doublets() and correct_ambient(). Saving includes their result if it is for the cells being saved.")),
                      plotOutput("extra_plot", height = "380px")),
            nav_panel("Annotation", value = "Annotation",
                      div(class = "d-flex flex-wrap gap-2 align-items-center mb-2",
                          busy_button("run_annot", "Annotate the cells that pass (preview)", "Annotating…", icon = icon("tags"), class = "btn-outline-primary"),
                          span(class = "form-text", "Uses the settings in the Annotation panel (and the sample's sample-sheet row). Saving includes this result if it was made for the cells being saved.")),
                      uiOutput("annot_status"),
                      plotOutput("annot_plot", height = "auto"),
                      tableOutput("annot_table")),
            nav_panel("Summary",
                      layout_columns(col_widths = c(5, 7),
                                     card(card_header("QC result (why each was removed)"), tableOutput("reason_table")),
                                     card(card_header("Cutoffs"), tableOutput("cutoff_table"))),
                      card(card_header("Medians"), tableOutput("median_table"))),
            nav_panel("Dataset", value = "Dataset",
                      div(class = "form-text mb-2", "Every sample in the dataset folder or sample sheet. Kept counts update when a sample is saved or viewed."),
                      tableOutput("dataset_table")),
            nav_panel("Log", verbatimTextOutput("log_text"))
          )))
)

## =============================================================================
## Server
## =============================================================================
server <- function(input, output, session) {
  # Active sample (the one on screen)
  rv <- reactiveValues(data = NULL, p_loaded = NULL, mod = NULL, slider_cfg = NULL, slider_set = list(), slider_exact = list(),
                       extra = NULL, active = NULL, load_n = 0, qc = FALSE, looks = NULL,
                       # dataset
                       samples = NULL, batch_n = 0, batch_info = NULL, pending_action = NULL,
                       # misc
                       log = character(0), prefill = NULL, sheet = NULL, auto_dataset = NULL,
                       # patient / sample metadata, annotation
                       sample_meta = NULL, prefill_patient = NULL, annot_n = 0,
                       # dataset-level runs
                       queue = NULL, queue_label = NULL, queue_go = 0, chain = NULL)
  # Loaded samples' heavy data live here (not reactive, so switching never copies them):
  # bank[[id]] = list(data, p_loaded, qc, looks, slider_cfg, slider_values, slider_set, extra)
  # qc = FALSE until you click Apply QC: the data are only read and measured until then.
  bank <- new.env()
  # A dataset loaded from a sample sheet: rows (by sample_id), the sheet's folder and file
  ds_sheet <- new.env()

  add_log <- function(x) isolate(rv$log <- tail(c(rv$log, x), 3000))
  logged <- function(expr) withCallingHandlers(expr,
    message = function(m) add_log(sub("\n$", "", conditionMessage(m))),
    warning = function(w) { add_log(paste("Warning:", conditionMessage(w))); invokeRestart("muffleWarning") })
  fail <- function(e, what) {
    msg <- if (inherits(e, "condition")) conditionMessage(e) else as.character(e)
    add_log(paste0(what, " failed: ", msg))
    showNotification(div(strong(what, " failed"), br(), pre(style = "white-space: pre-wrap; margin: 0", msg)),
                     type = "error", duration = NULL)
  }

  engine <- reactive(tryCatch(get_engine(input$modality), error = function(e) e))
  mod <- reactive(MODALITIES[[input$modality]])
  loaded_ok <- reactive(!is.null(rv$data) && identical(rv$mod, input$modality))
  qc_ok <- reactive(loaded_ok() && isTRUE(rv$qc))       # QC applied to the sample on screen
  has_qc <- function(id) exists(id, envir = bank, inherits = FALSE) && isTRUE(bank[[id]]$qc)
  need_qc <- function(what) showNotification(sprintf("Apply QC first (the button in the QC cutoffs panel); %s needs the cutoffs.", what),
                                             type = "warning")
  in_memory <- function() ls(bank)

  set_sample <- function(id, ...) {
    s <- isolate(rv$samples); if (is.null(s)) return()
    i <- match(id, s$sample_id); if (is.na(i)) return()
    v <- list(...)
    for (nm in names(v)) s[[nm]][i] <- v[[nm]]
    rv$samples <- s
  }

  clear_all <- function() {
    rm(list = ls(bank), envir = bank)
    rm(list = ls(ds_sheet), envir = ds_sheet)
    rv$data <- NULL; rv$extra <- NULL; rv$slider_cfg <- NULL; rv$active <- NULL; rv$p_loaded <- NULL; rv$qc <- FALSE
    rv$samples <- NULL; rv$batch_n <- 0; rv$batch_info <- NULL
    gc(verbose = FALSE)
  }

  ## ---- modality / mode switches ----------------------------------------------------
  observeEvent(input$modality, {
    clear_all(); rv$prefill <- NULL; rv$sheet <- NULL
    m <- mod()
    if (m$spatial) nav_show("tabs", "Spatial") else nav_hide("tabs", "Spatial")
    if (isTRUE(m$extra_step)) nav_show("tabs", "Doublets & ambient") else nav_hide("tabs", "Doublets & ambient")
  })
  observe({ if (identical(input$input_mode, "dataset")) nav_show("tabs", "Dataset") else nav_hide("tabs", "Dataset") })
  output$path_hint <- renderText({
    if (identical(input$input_mode, "dataset"))
      paste("A folder with one sub-folder per sample (for scRNA, one count file per sample also works; sub-folder names become sample IDs),",
            "or a sample sheet .csv: every row is loaded, each with its own settings (QC settings, cutoff_* columns, patient,",
            "annot_* annotation settings). An output folder's samplesheet.csv reproduces that run.")
    else mod()$path_hint
  })
  output$load_button <- renderUI({
    if (identical(input$input_mode, "dataset"))
      busy_button("load", "Load the dataset (first batch)", "Finding the samples and loading the first batch…", icon = icon("play"), class = "btn-primary")
    else busy_button("load", "Load sample", "Loading the sample…", icon = icon("play"), class = "btn-primary")
  })
  output$data_box <- renderUI({
    s <- rv$samples
    if (is.null(s) || is.null(rv$mod)) return(div(class = "text-muted", icon("circle-info"), " No data loaded yet."))
    n_qc <- sum(vapply(in_memory(), has_qc, TRUE)); rv$qc
    tagList(div(class = "fw-semibold", MODALITIES[[rv$mod]]$label, " · ",
                if (nrow(s) > 1) sprintf("%s (%d samples)", paste(unique(s$dataset %||% input$dataset), collapse = ", "), nrow(s)) else s$sample_id[1]),
            div(class = "text-muted small", if (nrow(s) > 1) sprintf("%d in memory, QC applied to %d", length(in_memory()), n_qc)
                                            else if (n_qc) "QC applied" else "QC not applied yet"))
  })
  output$script_info <- renderUI({
    env <- engine()
    if (inherits(env, "error")) return(span(class = "text-danger", conditionMessage(env)))
    tagList("Using ", code(env$.script_path), br(), "modified ", format(env$.script_mtime, "%Y-%m-%d %H:%M"))
  })

  ## ---- file browsing ----------------------------------------------------------------
  if (HAS_SHINYFILES) {
    roots <- c("bladder cancer" = dirname(APP_DIR), Home = normalizePath("~", winslash = "/"),
               shinyFiles::getVolumes()())
    roots <- roots[!duplicated(roots)]
    shinyFiles::shinyDirChoose(input, "browse_in_dir", roots = roots, session = session)
    shinyFiles::shinyFileChoose(input, "browse_in_file", roots = roots, session = session)
    shinyFiles::shinyFileChoose(input, "browse_sheet", roots = roots, session = session, filetypes = c("csv"))
    shinyFiles::shinyDirChoose(input, "browse_out", roots = roots, session = session)
    picked_dir  <- function(x) if (is.list(x)) shinyFiles::parseDirPath(roots, x) else character(0)
    picked_file <- function(x) if (is.list(x)) shinyFiles::parseFilePaths(roots, x)$datapath else character(0)
    observeEvent(input$browse_in_dir,  { p <- picked_dir(input$browse_in_dir);   if (length(p)) updateTextInput(session, "path", value = p) })
    observeEvent(input$browse_in_file, { p <- picked_file(input$browse_in_file); if (length(p)) updateTextInput(session, "path", value = p) })
    observeEvent(input$browse_sheet,   { p <- picked_file(input$browse_sheet);   if (length(p)) updateTextInput(session, "sheet_path", value = p) })
    observeEvent(input$browse_out,     { p <- picked_dir(input$browse_out);      if (length(p)) updateTextInput(session, "out_dir", value = p) })
    shinyFiles::shinyFileChoose(input, "browse_meta", roots = roots, session = session, filetypes = c("csv"))
    observeEvent(input$browse_meta,    { p <- picked_file(input$browse_meta);    if (length(p)) updateTextInput(session, "meta_path", value = p) })
  }
  # Suggest a sample ID (one-sample mode) or dataset name (dataset mode) from the path
  observeEvent(input$path, {
    p <- gsub("\\\\", "/", clean_path(input$path))
    if (!nzchar(p)) return()
    if (identical(input$input_mode, "dataset")) {
      if (grepl("\\.csv$", p, ignore.case = TRUE)) return()   # a sample sheet: its rows have their own dataset
      if (input$dataset %in% c("", "dataset1") || identical(input$dataset, rv$auto_dataset)) {
        rv$auto_dataset <- dataset_name_from_dir(p)
        updateTextInput(session, "dataset", value = rv$auto_dataset)
      }
    } else if (!nzchar(input$sample_id %||% "")) {
      b <- basename(p)
      if (tolower(b) %in% c("outs", "filtered_feature_bc_matrix", "raw_feature_bc_matrix")) b <- basename(dirname(p))
      b <- sub("\\.(h5|h5ad|rds|csv|tsv|txt)(\\.gz)?$", "", b, ignore.case = TRUE)
      b <- sub("_(qc|annotated)$", "", b)
      updateTextInput(session, "sample_id", value = sanitize_id(b))
    }
  }, ignoreInit = TRUE)

  ## ---- sample sheet (one-sample mode) ------------------------------------------------
  observeEvent(input$sheet_path, {
    f <- clean_path(input$sheet_path)
    rv$sheet <- NULL
    if (!nzchar(f)) return()
    if (!file.exists(f)) { showNotification("Sample sheet not found", type = "warning"); return() }
    sh <- tryCatch({
      s <- read.csv(f, colClasses = "character", na.strings = c("", "NA"), check.names = FALSE, strip.white = TRUE)
      names(s) <- tolower(trimws(names(s)))
      s <- s[rowSums(!is.na(s)) > 0, , drop = FALSE]
      if (!all(c("sample_id", "path") %in% names(s))) stop("needs sample_id and path columns")
      s
    }, error = function(e) { fail(e, "Reading the sample sheet"); NULL })
    if (!is.null(sh)) rv$sheet <- list(df = sh, dir = dirname(normalizePath(f, winslash = "/")))
  })
  output$sheet_rows_ui <- renderUI({
    req(rv$sheet)
    tagList(
      selectInput("sheet_row", "Sample", choices = rv$sheet$df$sample_id),
      actionButton("use_row", "Fill in this sample's settings", class = "btn-sm btn-outline-primary w-100"))
  })
  observeEvent(input$use_row, {
    env <- engine(); req(!inherits(env, "error"), rv$sheet, input$sheet_row)
    i <- match(input$sheet_row, rv$sheet$df$sample_id)
    res <- resolve_row(env, rv$sheet$df[i, , drop = FALSE], i, rv$sheet$dir)
    if (length(res$problems)) showNotification(paste(res$problems, collapse = "\n"), type = "warning", duration = 15)
    p <- res$params
    updateRadioButtons(session, "input_mode", selected = "single")
    updateTextInput(session, "path", value = p$path)
    updateTextInput(session, "sample_id", value = p$sample_id)
    updateTextInput(session, "dataset", value = p$dataset)
    rv$prefill <- p
    rv$prefill_row <- list(row = rv$sheet$df[i, , drop = FALSE], dir = rv$sheet$dir)   # its cutoff_* and other columns, used on loading
    pcol <- intersect(c("patient_id", "patient"), names(rv$sheet$df))
    rv$prefill_patient <- if (length(pcol) && !is.na(rv$sheet$df[[pcol[1]]][i])) rv$sheet$df[[pcol[1]]][i] else NULL
    if (!nzchar(input$out_dir %||% "")) updateTextInput(session, "out_dir", value = file.path(rv$sheet$dir, "qc_output"))
    showNotification(sprintf("Filled in settings for %s. Click Load sample.", p$sample_id), type = "message")
  })

  ## ---- advanced settings --------------------------------------------------------------
  # The script's sample-sheet parameters, split in two: QC filtering parameters (under the
  # QC cutoffs sliders) and everything else (Advanced settings panel)
  param_groups_ui <- function(qc) {
    env <- engine()
    if (inherits(env, "error")) return(div(class = "text-danger", conditionMessage(env)))
    vals <- rv$prefill %||% env$DEFAULTS
    groups <- mod()$groups
    other <- setdiff(names(env$DEFAULTS), c(unlist(groups), "dataset"))
    if (length(other)) groups[["Other"]] <- other
    groups <- groups[if (qc) names(groups) %in% QC_GROUPS else !names(groups) %in% QC_GROUPS]
    tagList(lapply(names(groups), function(g) {
      nms <- intersect(groups[[g]], names(env$DEFAULTS))
      if (!length(nms)) return(NULL)
      tagList(h6(g, class = "mt-3 mb-2"),
              lapply(nms, function(nm) param_input(nm, vals[[nm]], env, env$.param_help[[nm]])))
    }))
  }
  output$advanced_qc_ui <- renderUI(param_groups_ui(qc = TRUE))
  output$advanced_ui    <- renderUI(param_groups_ui(qc = FALSE))
  # Render even while collapsed, so the inputs always exist
  outputOptions(output, "advanced_qc_ui", suspendWhenHidden = FALSE)
  outputOptions(output, "advanced_ui", suspendWhenHidden = FALSE)

  # A one-row sample sheet from the inputs, parsed by the script's own resolve_sample().
  # sheet: list(row, dir, pin) for a sample from a sample sheet. pin = TRUE (a dataset loaded from a
  # sheet): the row's values win over the inputs, and its QC settings stay with the sample
  # (attr(p, "row_set")). pin = FALSE (one sample filled in from a sheet; the inputs hold its
  # settings): only the row's other columns are used (cutoff_*, patient, notes, annot_*, ...).
  params_for <- function(sample_id, path, dataset = input$dataset, sheet = NULL) {
    env <- engine()
    row <- list(sample_id = trimws(sample_id %||% ""), path = clean_path(path), dataset = trimws(dataset %||% ""))
    for (nm in setdiff(names(env$DEFAULTS), "dataset")) {
      v <- input[[paste0("p_", nm)]]
      if (is.null(v)) v <- env$DEFAULTS[[nm]]
      v <- as.character(v)
      row[[nm]] <- if (!length(v) || is.na(v) || !nzchar(trimws(v))) NA_character_ else v
    }
    for (k in c("sample_id", "path", "dataset")) if (!nzchar(row[[k]])) row[[k]] <- NA_character_
    dir <- APP_DIR
    if (!is.null(sheet)) {
      dir <- sheet$dir
      sr <- sheet$row
      cols <- if (isTRUE(sheet$pin)) names(sr) else setdiff(names(sr), c("sample_id", "path", "dataset", names(env$DEFAULTS)))
      for (nm in cols) { v <- sr[[nm]]; if (length(v) == 1 && !is.na(v) && nzchar(trimws(v))) row[[nm]] <- as.character(v) }
      if (isTRUE(sheet$pin)) row$path <- gsub("\\\\", "/", sr$path)   # relative to the sheet's folder
    }
    res <- resolve_row(env, as.data.frame(row, stringsAsFactors = FALSE, check.names = FALSE), 1, dir)
    if (isTRUE(sheet$pin)) attr(res$params, "row_set") <- c(sheet_set_params(env, sheet$row), "dataset")
    res
  }
  # The QC settings a sample's sheet row sets, for messages: the ones that differ from the defaults
  row_set_text <- function(env, p) {
    rs <- setdiff(attr(p, "row_set"), c("dataset", "format"))
    dif <- rs[!vapply(rs, function(nm) same_val(env$DEFAULTS[[nm]], p[[nm]]), TRUE)]
    txt <- function(nms) paste(sprintf("%s = %s", nms, vapply(nms, function(nm) fmt_val(p[[nm]]), "")), collapse = ", ")
    if (length(rs) <= 6) return(txt(rs))
    sprintf("%d QC settings (%s)", length(rs), if (length(dif)) paste("not the defaults:", txt(dif)) else "all the defaults")
  }
  # The sample-sheet row a sample was (or will be) loaded from, as params_for()'s `sheet`
  sheet_for <- function(id) {
    if (!is.null(ds_sheet$rows[[id]])) return(list(row = ds_sheet$rows[[id]], dir = ds_sheet$dir, pin = TRUE))
    pr <- isolate(rv$prefill_row)
    if (!is.null(pr) && identical(isolate(input$input_mode), "single") && identical(pr$row$sample_id, id))
      return(list(row = pr$row, dir = pr$dir, pin = FALSE))
    NULL
  }
  current_params_raw <- reactive({
    env <- engine(); req(!inherits(env, "error"))
    params_for(input$sample_id, input$path)
  })
  current_params <- debounce(current_params_raw, 400)

  load_keys <- reactive(c(if (identical(input$input_mode, "single")) c("sample_id", "path"), load_params(input$modality)))

  # A sample's loaded parameters + the live (non-loading) parameters from the inputs, except the
  # ones its sample-sheet row sets
  merge_live <- function(p_loaded, cur) {
    p <- p_loaded
    m <- isolate(input$modality)
    keep <- c("sample_id", "path", load_params(m), attr(p_loaded, "row_set"))
    for (nm in setdiff(intersect(names(cur), names(get_engine(m)$DEFAULTS)), keep)) p[[nm]] <- cur[[nm]]
    p
  }
  p_live <- reactive({ req(rv$p_loaded); merge_live(rv$p_loaded, current_params()$params) })


  ## ---- patient ID and sample metadata ------------------------------------------------------
  output$patient_ui <- renderUI({
    id <- rv$active
    if (is.null(id)) return(div(class = "form-text", "Load a sample to set its patient ID."))
    s <- isolate(rv$samples); val <- s$patient[match(id, s$sample_id)]
    tagList(textInput("patient_id", sprintf("Patient ID for %s", id), value = if (is.na(val)) "" else val, placeholder = "e.g. PT01"),
            div(class = "form-text param-help mt-0",
                "Added to every cell as ", code("patient"), ", so samples from one patient (e.g. tumour + adjacent normal) ",
                "can be grouped later or used as the Harmony batch."))
  })
  outputOptions(output, "patient_ui", suspendWhenHidden = FALSE)
  observeEvent(input$patient_id, {
    id <- rv$active; req(id)
    set_sample(id, patient = trimws(input$patient_id))
  }, ignoreInit = TRUE)

  # A table of per-sample metadata: sample_id + patient_id (+ any other columns)
  observeEvent(input$meta_path, {
    f <- clean_path(input$meta_path)
    if (!nzchar(f)) { rv$sample_meta <- NULL; return() }
    if (!file.exists(f)) { showNotification("Sample metadata file not found", type = "warning"); return() }
    mt <- tryCatch({
      d <- read.csv(f, colClasses = "character", na.strings = c("", "NA"), check.names = FALSE, strip.white = TRUE)
      names(d) <- tolower(trimws(names(d)))
      if (!"sample_id" %in% names(d)) stop("needs a sample_id column")
      names(d)[names(d) == "patient_id"] <- "patient"
      d
    }, error = function(e) { fail(e, "Reading the sample metadata"); NULL })
    rv$sample_meta <- mt
    apply_sample_meta()
  })
  apply_sample_meta <- function(quiet = FALSE) {
    mt <- isolate(rv$sample_meta); s <- isolate(rv$samples)
    if (is.null(mt) || is.null(s)) return()
    i <- match(s$sample_id, mt$sample_id)
    if ("patient" %in% names(mt)) s$patient[!is.na(i)] <- ifelse(is.na(mt$patient[i[!is.na(i)]]), "", mt$patient[i[!is.na(i)]])
    rv$samples <- s
    if (!quiet) showNotification(sprintf("Sample metadata matched %d of %d samples%s.", sum(!is.na(i)), nrow(s),
                                         if (any(is.na(i))) paste0(" (no row for: ", paste(head(s$sample_id[is.na(i)], 8), collapse = ", "), ")") else ""),
                                 type = if (any(is.na(i))) "warning" else "message", duration = 10)
    id <- isolate(rv$active)
    if (!is.null(id)) updateTextInput(session, "patient_id", value = s$patient[match(id, s$sample_id)])
  }
  # Everything known about one sample, added to every cell when it is saved
  sample_meta_for <- function(id) {
    s <- isolate(rv$samples); i <- match(id, s$sample_id)
    out <- list(patient = if (!is.na(i) && nzchar(s$patient[i] %||% "")) s$patient[i] else NA_character_)
    mt <- isolate(rv$sample_meta)
    if (!is.null(mt) && id %in% mt$sample_id)
      for (nm in setdiff(names(mt), c("sample_id", "patient", "path", "dataset"))) out[[make.names(nm)]] <- mt[[nm]][match(id, mt$sample_id)]
    out
  }

  ## ---- annotation settings ----------------------------------------------------------------
  aengine <- reactive(tryCatch(get_engine("annotation"), error = function(e) e))
  output$annotation_ui <- renderUI({
    env <- aengine()
    if (inherits(env, "error")) return(div(class = "text-danger small", conditionMessage(env)))
    tagList(
      div(class = "form-text mb-2", "Same parameters and defaults as ", code(basename(env$.script_path)),
          ". Each section shows the main settings; the rest are under its Advanced."),
      annot_settings_ui(env))
  })
  outputOptions(output, "annotation_ui", suspendWhenHidden = FALSE)
  # "(2 changed)" on an Advanced toggle when settings hidden there differ from the defaults
  for (sec in ANNOT_SECTIONS) local({
    sid <- sec$id
    output[[paste0("annot_adv_n_", sid)]] <- renderText({
      env <- aengine(); req(!inherits(env, "error"))
      p <- collect_annotation_params(input, env)
      ch <- annot_adv_changed(env, p, sid)
      if (length(ch)) sprintf("(%d changed)", length(ch)) else ""
    })
    outputOptions(output, paste0("annot_adv_n_", sid), suspendWhenHidden = FALSE)
  })
  annot_params <- function() collect_annotation_params(input, isolate(aengine()))

  ## ---- loading -------------------------------------------------------------------------
  # Read one sample into the bank. Returns TRUE/FALSE.
  load_one <- function(id, path, dataset) {
    m <- isolate(input$modality)
    env <- tryCatch(get_engine(m), error = function(e) { fail(e, "Loading the QC script"); NULL })
    if (is.null(env)) return(FALSE)
    res <- params_for(id, path, dataset, sheet_for(id))
    if (length(res$problems)) {
      add_log(paste(res$problems, collapse = "\n")); set_sample(id, status = "failed", note = paste(res$problems, collapse = "; "))
      return(FALSE)
    }
    p <- res$params
    add_log(sprintf("---- loading %s (%s) ----", id, MODALITIES[[m]]$label))
    d <- tryCatch(logged(load_data(m, env, p)), error = function(e) e)
    if (inherits(d, "error")) {
      add_log(paste("Loading failed:", conditionMessage(d)))
      set_sample(id, status = "failed", note = conditionMessage(d)); return(FALSE)
    }
    add_log(sprintf("[%s] input looks like %s", id, d$looks$text))
    tc <- take_cutoffs(env, p)            # cutoff_* columns of its sample-sheet row: the sliders start there
    if (length(setdiff(attr(p, "row_set"), "dataset")) || length(tc$cutoffs))
      add_log(sprintf("[%s] from its sample-sheet row: %s%s", id, row_set_text(env, p),
                      if (length(tc$cutoffs)) paste0("; cutoffs ", paste(sprintf("%s = %s", names(tc$cutoffs), vapply(tc$cutoffs, fmt_val, "")), collapse = ", ")) else ""))
    bank[[id]] <- list(data = d, p_loaded = tc$p, row_cutoffs = tc$cutoffs, qc = FALSE, looks = d$looks,
                       slider_cfg = NULL, slider_values = NULL, slider_set = NULL, extra = NULL)
    set_sample(id, status = "loaded", n_in = nrow(d$meta), n_kept = NA_integer_, looks = d$looks$short, note = "")
    TRUE
  }

  # Apply QC to one loaded sample: cell calling for scRNA (prepare_qc), then the MAD cutoffs.
  # Uses the settings the sample was loaded with (change those, and you load again).
  apply_qc <- function(id) {
    if (!exists(id, envir = bank, inherits = FALSE) || has_qc(id)) return(TRUE)
    m <- isolate(input$modality); env <- get_engine(m); e <- bank[[id]]
    add_log(sprintf("---- applying QC to %s ----", id))
    d <- tryCatch(logged(prepare_qc(m, env, e$data, e$p_loaded)), error = function(err) err)
    if (inherits(d, "error")) { fail(d, paste("Applying QC to", id)); return(FALSE) }
    finish_qc(id, d)
  }
  # The part after cell calling (also used when cell calling ran in the background)
  finish_qc <- function(id, d) {
    e <- bank[[id]]
    e$data <- d; e$qc <- TRUE
    e$slider_cfg <- NULL; e$slider_values <- NULL; e$slider_set <- NULL; e$slider_exact <- NULL; e$extra <- NULL; e$annot <- NULL
    bank[[id]] <- e
    st <- tryCatch(sample_state(id), error = function(err) { fail(err, paste("Computing cutoffs for", id)); NULL })
    set_sample(id, n_in = nrow(d$meta), n_kept = if (is.null(st)) NA_integer_ else st$n_kept)
    if (identical(id, isolate(rv$active))) activate(id, force = TRUE)
    TRUE
  }

  # Remember the on-screen sample's slider positions and doublet results
  stash_active <- function() {
    id <- isolate(rv$active)
    if (is.null(id) || !exists(id, envir = bank, inherits = FALSE)) return()
    e <- bank[[id]]
    cfg <- isolate(rv$slider_cfg)
    e$slider_cfg <- cfg
    e$slider_values <- lapply(cfg, function(c) isolate(input[[c$input_id]]) %||% c$value)
    e$slider_set <- isolate(rv$slider_set)
    e$slider_exact <- isolate(rv$slider_exact)
    e$extra <- isolate(rv$extra)
    bank[[id]] <- e
    if (!isTRUE(e$qc)) return()
    st <- tryCatch(sample_state(id), error = function(e) NULL)
    if (!is.null(st)) set_sample(id, n_kept = st$n_kept, manual = paste(names(st$sources)[st$sources == "manual"], collapse = ", "))
  }

  # Put a loaded sample on screen
  activate <- function(id, force = FALSE) {
    if (!exists(id, envir = bank, inherits = FALSE)) return()
    if (!identical(isolate(rv$active), id)) stash_active()
    m <- isolate(input$modality); env <- get_engine(m); e <- bank[[id]]
    if (isTRUE(e$qc)) {
      p <- merge_live(e$p_loaded, isolate(current_params()$params))
      th <- auto_thresholds(m, env, e$data, p)
      rv$load_n <- isolate(rv$load_n) + 1
      sl <- restore_sliders(MODALITIES[[m]]$sliders, e$data$meta, th, isolate(rv$load_n), e$slider_cfg, e$slider_values, e$slider_set)
      sl$exact <- e$slider_exact %||% list()
      if (is.null(e$slider_cfg) && length(e$row_cutoffs)) {   # first time on screen: the sample sheet's cutoffs
        sc <- sheet_cutoff_sliders(sl$cfg, e$row_cutoffs); sl$cfg <- sc$cfg; sl$exact <- sc$exact
      }
    } else sl <- list(cfg = NULL, set = list(), exact = list())     # QC not applied: no sliders yet
    rv$qc <- isTRUE(e$qc)
    rv$looks <- e$looks
    rv$mod <- m
    rv$p_loaded <- e$p_loaded
    rv$extra <- e$extra
    rv$slider_set <- sl$set
    rv$slider_exact <- sl$exact
    rv$slider_cfg <- sl$cfg
    rv$data <- e$data
    rv$active <- id
    if (MODALITIES[[m]]$spatial) {
      cols <- intersect(c("nCount_Spatial", "nFeature_Spatial", "percent.mt", "percent.ribo", "percent.hb",
                          "nCount_Xenium", "nFeature_Xenium", "cell_area", "nucleus_area", "pct_control"), names(e$data$meta))
      updateSelectInput(session, "sp_colour", choices = c("QC result" = "qc_reason", setNames(cols, cols)),
                        selected = isolate(input$sp_colour) %||% "qc_reason")
    }
  }

  # Unload the current batch and load the next one: samples one at a time until the
  # batch cap is reached or free RAM would drop below RAM_RESERVE_GB.
  load_next_batch <- function() {
    if (!is.null(isolate(rv$queue))) { showNotification("Wait for the dataset-level run to finish (or stop it) first.", type = "warning"); return() }
    s <- isolate(rv$samples); req(s)
    for (id in in_memory()) if (s$status[match(id, s$sample_id)] == "loaded") set_sample(id, status = "unsaved")
    rm(list = ls(bank), envir = bank)
    rv$data <- NULL; rv$active <- NULL; rv$slider_cfg <- NULL; rv$extra <- NULL
    gc(verbose = FALSE)
    s <- isolate(rv$samples)
    todo <- c(s$sample_id[s$status == "pending"], s$sample_id[s$status == "unsaved"])
    if (!length(todo)) { showNotification("No samples left to load.", type = "message"); return() }
    cap <- max(1, as.integer(isolate(input$max_batch) %||% MAX_BATCH_DEFAULT))
    loaded <- character(0); biggest <- 0; stop_reason <- NULL
    withProgress(message = "Loading samples", value = 0, {
      for (id in todo) {
        if (length(loaded) >= cap) { stop_reason <- sprintf("batch size %d reached", cap); break }
        if (length(loaded) >= 1) {
          free <- free_ram_gb()
          if (!is.na(free) && free - 1.5 * biggest < RAM_RESERVE_GB) {
            stop_reason <- sprintf("free RAM %.1f GB; the next sample needs ~%.1f GB", free, 1.5 * biggest); break
          }
        }
        incProgress(1 / min(cap, length(todo)), detail = id)
        before <- r_mem_used_gb()
        ok <- load_one(id, s$path[match(id, s$sample_id)], s$dataset[match(id, s$sample_id)])
        biggest <- max(biggest, r_mem_used_gb() - before)
        if (ok) loaded <- c(loaded, id)
      }
    })
    rv$batch_n <- isolate(rv$batch_n) + 1
    left <- sum(isolate(rv$samples)$status %in% c("pending", "unsaved"))
    rv$batch_info <- sprintf("Batch %d: %d sample(s) in memory (R is using %.1f GB)%s. %d left to load.",
                             isolate(rv$batch_n), length(loaded), r_mem_used_gb(),
                             if (!is.null(stop_reason) && left) paste0("; stopped at ", stop_reason) else "", left)
    add_log(isolate(rv$batch_info))
    if (length(loaded)) activate(loaded[1])
    else showNotification("No samples in this batch could be loaded; see the Dataset tab and the Log.", type = "warning")
  }

  observeEvent(input$load, {
    env <- engine()
    if (inherits(env, "error")) { fail(env, "Loading the QC script"); return() }
    # Visium data that isn't in Space Ranger folders is turned into them first, written next to
    # the data, then loaded as a dataset:
    #   - a combined .h5ad (several sections in one file, e.g. GSE171351) is split into one folder
    #     per section (R/h5ad.R)
    #   - GEO loose files (counts .csv + spot_coordinates .csv + image per sample, e.g. GSE246011)
    #     become one folder per sample (R/geo_visium.R)
    # Picked as a file: that file's sections / sample only. A folder: everything in it is converted,
    # and the folder scan then finds the new folders (a one-sample folder works in either mode).
    mode <- input$input_mode
    split_sheet <- NULL; split_out <- NULL; sheet_shown <- TRUE
    if (identical(input$modality, "visium")) {
      p0 <- clean_path(input$path)
      is_file <- nzchar(p0) && file.exists(p0) && !dir.exists(p0)
      is_dir <- nzchar(p0) && dir.exists(p0)
      targets <- if (is_file) p0 else if (identical(mode, "dataset") && is_dir)
        list.files(p0, pattern = "\\.h5ad$", full.names = TRUE, ignore.case = TRUE) else character(0)
      targets <- Filter(is_combined_visium_h5ad, targets)
      for (f in targets) {
        sp <- tryCatch(withProgress(message = sprintf("Splitting %s into sections", basename(f)), value = 0.05,
                                    split_visium_h5ad(f, progress = function(k, n, nm) setProgress(k / (n + 1), detail = nm))),
                       error = function(e) e)
        if (inherits(sp, "error")) { fail(sp, sprintf("Splitting %s", basename(f))); return() }
        add_log(sp$log)
        showNotification(sprintf("%s holds %d Visium sections: %s %s next to it. Loading them as dataset %s.", basename(f),
                                 nrow(sp$samples), if (sp$reused) "using the folders already split" else "split into one folder each",
                                 paste0("(", basename(sp$sheet), ")"), sp$dataset), type = "message", duration = 12)
        if (is_file) { split_sheet <- sp$sheet; split_out <- file.path(dirname(dirname(sp$sheet)), "qc_output") }
      }
      geo_dir <- if (is_file && !length(targets)) dirname(p0) else if (is_dir) p0 else ""
      if (nzchar(geo_dir)) {
        gv <- tryCatch(withProgress(message = "Converting GEO Visium files into Space Ranger folders", value = 0.05,
                                    convert_geo_visium_folder(geo_dir, only = if (is_file) p0,
                                                              progress = function(k, n, d) setProgress((k - 0.5) / n, detail = d))),
                       error = function(e) e)
        if (inherits(gv, "error")) { fail(gv, "Converting the GEO Visium files"); return() }
        if (!is.null(gv)) {
          add_log(gv$log)
          weak <- grep("WARNING", gv$log, value = TRUE)
          showNotification(div(sprintf("Converted %d GEO Visium sample(s) into Space Ranger folders next to the files (%s). The image scale was fitted to the spots: check spatial/alignment_check_app.png in each sample folder.",
                                       nrow(gv$samples), basename(gv$sheet)),
                               if (length(weak)) div(class = "text-danger", "Weak image fit: ", paste(sub("^\\[([^]]+)\\].*", "\\1", weak), collapse = ", "))),
                           type = if (length(weak)) "warning" else "message", duration = 20)
          if (is_file || identical(mode, "single")) {
            split_out <- file.path(dirname(normalizePath(geo_dir, winslash = "/")), "qc_output")
            if (is_file) {   # just this sample: a one-row sheet
              split_sheet <- tempfile("visium_converted_", fileext = ".csv"); sheet_shown <- FALSE
              one <- gv$samples; one$path <- file.path(normalizePath(geo_dir, winslash = "/"), one$path)
              write.csv(one, split_sheet, row.names = FALSE)
            } else split_sheet <- gv$sheet
          }
        }
      }
      if (!is.null(split_sheet)) {
        mode <- "dataset"
        updateRadioButtons(session, "input_mode", selected = "dataset")
        if (sheet_shown) updateTextInput(session, "path", value = split_sheet)
        if (!nzchar(clean_path(input$out_dir))) updateTextInput(session, "out_dir", value = split_out)
      }
    }
    if (identical(mode, "single")) {
      res <- current_params_raw()
      if (length(res$problems)) { showNotification(paste(res$problems, collapse = "\n"), type = "error", duration = 15); return() }
      p <- res$params
      if (input$modality %in% c("scrna", "xenium")) {
        chk <- tryCatch({ env$check_packages(list(p)); NULL }, error = function(e) e)
        if (!is.null(chk)) { fail(chk, "Loading"); return() }
      }
      clear_all()
      rv$samples <- data.frame(sample_id = p$sample_id, path = p$path, kind = "", dataset = p$dataset, looks = "", status = "pending",
                               n_in = NA_integer_, n_kept = NA_integer_, manual = "", note = "",
                               patient = isolate(rv$prefill_patient) %||% "", annotated = "", job = "", stringsAsFactors = FALSE)
      ok <- withProgress(message = sprintf("Loading %s", p$sample_id), detail = "reading counts and computing QC metrics…",
                         value = 0.3, load_one(p$sample_id, p$path, p$dataset))
      if (!ok) { fail(isolate(rv$samples)$note[1], "Loading"); return() }
      activate(p$sample_id)
      session$sendCustomMessage("qc-overlay", list(id = "qc-load", show = FALSE))
      if (!nzchar(input$out_dir %||% "")) updateTextInput(session, "out_dir", value = file.path(dirname(p$path), "qc_output"))
      return()
    }
    # Dataset: a folder (one sub-folder per sample) or a sample sheet (one row per sample)
    dir <- split_sheet %||% clean_path(input$path)
    from_sheet <- grepl("\\.csv$", dir, ignore.case = TRUE) && file.exists(dir) && !dir.exists(dir)
    rows <- NULL
    if (from_sheet) {
      ds0 <- trimws(input$dataset %||% ""); if (!nzchar(ds0)) ds0 <- "dataset1"   # for rows without a dataset
      sh <- tryCatch(read_dataset_sheet(dir, ds0), error = function(e) e)
      if (inherits(sh, "error")) { fail(sh, "Reading the sample sheet"); return() }
      if (length(sh$problems)) showNotification(paste("Sample sheet:", paste(sh$problems, collapse = "; ")), type = "warning", duration = 15)
      if (!nrow(sh$samples)) { fail("The sample sheet has no samples.", "Reading the sample sheet"); return() }
      found <- list(samples = sh$samples, skipped = data.frame(name = character(0), reason = character(0)))
      rows <- sh$rows
    } else {
      if (!dir.exists(dir)) { showNotification("Dataset folder (or sample sheet .csv) not found.", type = "error"); return() }
      found <- find_samples(dir, input$modality)
      if (!nrow(found$samples)) {
        fail(paste0("No ", mod()$label, " samples found in ", dir, ".",
                    if (nrow(found$skipped)) paste0("\nSkipped: ", paste(sprintf("%s (%s)", found$skipped$name, found$skipped$reason), collapse = "; "))),
             "Scanning the folder")
        return()
      }
      if (nrow(found$skipped)) showNotification(paste0("Skipped ", nrow(found$skipped), " folder(s): ",
                                                       paste(sprintf("%s (%s)", found$skipped$name, found$skipped$reason), collapse = "; ")),
                                                type = "warning", duration = 15)
      found$samples$dataset <- input$dataset
    }
    # check the packages the samples need (and, for a sheet, every row's settings) before loading anything
    ps <- lapply(seq_len(nrow(found$samples)), function(i) {
      id <- found$samples$sample_id[i]
      params_for(id, found$samples$path[i], found$samples$dataset[i],
                 if (from_sheet) list(row = rows[[id]], dir = dirname(normalizePath(dir, winslash = "/")), pin = TRUE))
    })
    bad <- Filter(function(r) length(r$problems), ps)
    if (length(bad)) showNotification(div(strong("Some rows can't be loaded:"), pre(style = "white-space: pre-wrap; margin: 0",
                                          paste(unlist(lapply(bad, `[[`, "problems")), collapse = "\n"))), type = "warning", duration = NULL)
    good <- lapply(Filter(function(r) !length(r$problems), ps), `[[`, "params")
    if (input$modality %in% c("scrna", "xenium") && length(good)) {
      chk <- tryCatch({ env$check_packages(good); NULL }, error = function(e) e)
      if (!is.null(chk)) { fail(chk, "Loading"); return() }
    }
    out_dir <- clean_path(input$out_dir)
    if (!nzchar(out_dir) && !is.null(split_out)) out_dir <- split_out
    if (!nzchar(out_dir)) {
      out_dir <- if (from_sheet) env$default_out_dir(dir) else file.path(dirname(normalizePath(dir, winslash = "/")), "qc_output")
      updateTextInput(session, "out_dir", value = out_dir)
    }
    clear_all()
    if (from_sheet) { ds_sheet$rows <- rows; ds_sheet$dir <- dirname(normalizePath(dir, winslash = "/")); ds_sheet$file <- dir }
    patient <- if (from_sheet) vapply(found$samples$sample_id, function(id) {
      r <- rows[[id]]; pc <- intersect(c("patient", "patient_id"), names(r))
      if (length(pc) && !is.na(r[[pc[1]]])) r[[pc[1]]] else ""
    }, "", USE.NAMES = FALSE) else ""
    s <- data.frame(found$samples, looks = "", status = "pending", n_in = NA_integer_, n_kept = NA_integer_, manual = "", note = "",
                    patient = patient, annotated = "", job = "", stringsAsFactors = FALSE)
    for (r in bad) {
      i <- match(r$params$sample_id, s$sample_id)
      if (!is.na(i)) { s$status[i] <- "failed"; s$note[i] <- paste(r$problems, collapse = "; ") }
    }
    if (isTRUE(input$skip_done)) {
      done <- s$status == "pending" & file.exists(rds_path(out_dir, s$dataset, s$sample_id))
      s$status[done] <- "on_disk"
    }
    rv$samples <- s
    apply_sample_meta(quiet = TRUE)
    add_log(sprintf("---- %s: %d samples in %s (%d already saved) ----", if (from_sheet) "sample sheet" else paste("dataset", input$dataset),
                    nrow(s), dir, sum(s$status == "on_disk")))
    load_next_batch()
    if (length(in_memory())) session$sendCustomMessage("qc-overlay", list(id = "qc-load", show = FALSE))
  })

  observeEvent(input$pick_sample, {
    id <- input$pick_sample
    if (exists(id, envir = bank, inherits = FALSE)) activate(id)
  })

  ## ---- sliders ---------------------------------------------------------------------------
  output$cutoff_heading <- renderUI({
    req(rv$active)
    rs <- setdiff(attr(rv$p_loaded, "row_set"), "dataset")
    cut <- if (exists(rv$active, envir = bank, inherits = FALSE)) bank[[rv$active]]$row_cutoffs
    multi <- nrow(rv$samples %||% data.frame()) > 1
    if (!multi && !length(rs) && !length(cut)) return(NULL)
    div(class = "form-text mb-2",
        if (multi) tagList("Cutoffs for ", strong(rv$active), " (each sample keeps its own). "),
        if (length(rs)) tagList("Its sample-sheet row sets ", code(row_set_text(engine(), rv$p_loaded)),
                                ": the Advanced settings don't change those for this sample. "),
        if (length(cut)) tagList("The row's ", code(paste(sprintf("cutoff_%s = %s", names(cut), vapply(cut, fmt_val, "")), collapse = ", ")),
                                 " set the starting position of its sliders."))
  })
  output$sliders_ui <- renderUI({
    cfg <- rv$slider_cfg
    if (loaded_ok() && !isTRUE(rv$qc)) return(apply_qc_ui())
    if (is.null(cfg)) return(div(class = "form-text", "Load a sample to see its cutoffs."))
    tagList(lapply(cfg, function(c) tagList(
      sliderInput(c$input_id, c$label, min = c$min, max = c$max, value = c$value, step = c$step, ticks = FALSE),
      # type a value instead of dragging (Enter, or click elsewhere, to apply)
      div(class = "qc-cut-boxes", lapply(names(slider_sides(c)), function(side) tags$label(
        if (side == "lo") "min" else "max",
        tags$input(type = "number", class = "form-control qc-cut-box", step = "any", min = 0,
                   placeholder = if (side == "hi") "no max" else NULL, `data-slider` = c$input_id, `data-side` = side,
                   title = "Type a cutoff and press Enter")))),
      div(class = "slider-readout", uiOutput(paste0("readout_", c$input_id), inline = TRUE)))))
  })
  outputOptions(output, "sliders_ui", suspendWhenHidden = FALSE)

  output$reset_auto_ui <- renderUI({
    req(qc_ok())
    actionButton("reset_auto", "Reset to MAD cutoffs", icon = icon("rotate-left"), class = "btn-sm btn-outline-secondary w-100")
  })
  outputOptions(output, "reset_auto_ui", suspendWhenHidden = FALSE)
  apply_qc_ui <- function() {
    tagList(
      div(class = "form-text mb-2", "QC isn't applied to ", strong(rv$active), " yet: the plots show the data as loaded. ",
          if (identical(rv$mod, "scrna")) "Applying QC calls cells (emptyDrops, for a raw matrix) and sets the MAD cutoffs; " else
            "Applying QC sets the MAD cutoffs; ", "you can then move them with sliders here."),
      busy_button("apply_qc", apply_qc_label(), "Applying QC…", icon = icon("play"), class = "btn-primary btn-sm w-100 mb-2"))
  }
  apply_qc_label <- function() {
    if (isTRUE(input$dataset_level) && identical(input$input_mode, "dataset")) "Apply QC to all samples in memory" else "Apply QC"
  }
  # With dataset-level processing on, Apply QC goes to every sample in memory without QC (in the background)
  apply_qc_click <- function() {
    id <- rv$active; req(id)
    if (auto_on() && !is.null(isolate(rv$queue))) { showNotification("A background run is already going; wait for it or stop it first.", type = "warning"); return() }
    if (dataset_level_on()) {
      stash_active()
      rv$chain <- if (auto_on()) auto_steps() else NULL
      started <- start_queue("qc", "Applying QC")
      # every sample already had QC: go straight to the next steps
      if (!started && !is.null(isolate(rv$chain))) run_chain(Filter(has_qc, in_memory()))
      return()
    }
    if (sample_job(id) %in% c("queued", "running")) { showNotification(sprintf("%s is already queued in the dataset-level run.", id), type = "message"); return() }
    ok <- apply_qc(id)
    if (isTRUE(ok) && auto_on()) { rv$chain <- auto_steps(); run_chain(id) }
  }

  ## ---- automatic pipeline (Advanced settings -> Pipeline) -------------------------------------
  # After Apply QC: scDblFinder + DecontX (scRNA-seq), then the annotation preview. Each step is a
  # run of the dataset-level queue (on one sample, or on every sample in memory), so it runs in
  # the background and the next step starts when the queue finishes, on the samples that succeeded.
  auto_on <- function() isTRUE(isolate(input$auto_pipeline))
  CHAIN_LABELS <- c(extra = "scDblFinder + DecontX", annot = "Annotating")
  auto_steps <- function() {
    steps <- c(if (isTRUE(MODALITIES[[isolate(rv$mod) %||% isolate(input$modality)]]$extra_step)) "extra", "annot")
    aenv <- isolate(aengine())
    chk <- if (inherits(aenv, "error")) conditionMessage(aenv) else
      tryCatch({ aenv$check_packages(annot_params()); NULL }, error = function(e) explain_missing(conditionMessage(e)))
    if (!is.null(chk)) {
      steps <- setdiff(steps, "annot")
      showNotification(paste("Automatic pipeline: the annotation step is skipped.", chk), type = "warning", duration = 15)
    }
    steps
  }
  run_chain <- function(ids) {          # start the next step on ids; stops the pipeline when there's nothing left
    steps <- isolate(rv$chain)
    if (!length(steps) || !length(ids)) { rv$chain <- NULL; return(invisible(FALSE)) }
    rv$chain <- steps[-1]
    ok <- start_queue(steps[1], CHAIN_LABELS[[steps[1]]], only = ids)
    if (!isTRUE(ok)) rv$chain <- NULL
    invisible(ok)
  }
  observeEvent(input$apply_qc, apply_qc_click())
  observeEvent(input$apply_qc_banner, apply_qc_click())

  observeEvent(rv$slider_cfg, {
    ids <- unname(vapply(Filter(function(c) c$scale == "log", rv$slider_cfg), function(c) c$input_id, ""))
    if (length(ids)) session$onFlushed(function() session$sendCustomMessage("qc-log-sliders", as.list(ids)), once = TRUE)
    for (c in rv$slider_cfg) local({
      cc <- c
      output[[paste0("readout_", cc$input_id)]] <- renderUI({
        st <- th_state()
        parts <- lapply(c(cc$lo, cc$hi), function(nm) {
          man <- is_manual(st$sources, nm)
          lab <- if (grepl("^min", nm)) "min" else "max"
          span(class = if (man) "manual" else NULL,
               sprintf("%s %s%s", lab, fmt_num(st$th[[nm]]), if (man) sprintf(" (manual; MAD %s)", fmt_num(st$auto[[nm]])) else " (MAD)"))
        })
        if (length(parts) == 2) tagList(parts[[1]], " · ", parts[[2]]) else parts[[1]]
      })
    })
  })

  # Moving a slider from the server resets its handle labels to the raw log10 values; show the
  # real numbers again on the log-scale ones
  relabel_log_sliders <- function(cfg) {
    ids <- unname(vapply(Filter(function(c) c$scale == "log", cfg), function(c) c$input_id, ""))
    if (length(ids)) session$onFlushed(function() session$sendCustomMessage("qc-log-sliders", as.list(ids)), once = TRUE)
  }
  # A value typed into a box under a slider: move the slider there, and keep the exact value
  observeEvent(input$cut_box, {
    b <- input$cut_box; cfg <- rv$slider_cfg; req(cfg)
    c <- Filter(function(x) identical(x$input_id, b$slider), cfg)
    if (!length(c)) return()
    c <- c[[1]]; side <- b$side; i <- slider_sides(c)[[side]]; nm <- c[[side]]
    v <- suppressWarnings(as.numeric(b$value))
    if (is.na(v) && nzchar(trimws(b$value %||% ""))) { showNotification("Type a number.", type = "warning"); rv$box_tick <- isolate(rv$box_tick %||% 0) + 1; return() }
    if (is.na(v)) v <- if (side == "hi") Inf else 0                 # emptied box: no max / no min
    if (v < 0) v <- 0
    pos <- if (is.infinite(v)) c$max else slider_pos(c, v)
    cur <- isolate(input[[c$input_id]]) %||% c$value
    new <- cur; new[i] <- round(pos / c$step) * c$step
    if (length(new) == 2 && new[1] > new[2]) {
      showNotification("The min can't be above the max.", type = "warning"); rv$box_tick <- isolate(rv$box_tick %||% 0) + 1; return()
    }
    ex <- isolate(rv$slider_exact)
    ex[[nm]] <- if (is.finite(v) && pos > c$min && pos < c$max) v else NULL   # at the ends the slider itself says 0 / no max
    rv$slider_exact <- ex
    updateSliderInput(session, c$input_id, value = new)
    relabel_log_sliders(list(c))
  })
  # Keep the boxes showing the cutoffs in use
  observe({
    rv$box_tick
    st <- tryCatch(th_state(), error = function(e) NULL); cfg <- rv$slider_cfg
    if (is.null(st) || is.null(cfg)) return()
    vals <- unname(unlist(lapply(unname(cfg), function(c) lapply(names(slider_sides(c)), function(side) {
      x <- st$th[[c[[side]]]]
      list(slider = c$input_id, side = side,
           value = if (is.null(x) || !is.finite(x)) "" else if (c$scale == "log") round(x) else signif(x, 4))
    })), recursive = FALSE))
    session$onFlushed(function() session$sendCustomMessage("qc-cut-boxes", vals), once = TRUE)
  })

  cut_inputs <- debounce(reactive({
    cfg <- rv$slider_cfg; req(cfg)
    lapply(cfg, function(c) input[[c$input_id]])
  }), 250)

  auto_th <- reactive({
    req(qc_ok())
    tryCatch(auto_thresholds(rv$mod, engine(), rv$data, p_live()), error = function(e) { fail(e, "Computing cutoffs"); req(FALSE) })
  })
  th_state <- reactive(resolve_cutoffs(auto_th(), rv$slider_cfg, cut_inputs(), rv$slider_set, rv$slider_exact))

  # When a cutoff rule changes, sliders still on "auto" follow the new MAD values
  live_rule_values <- reactive({
    p <- current_params()$params
    p[intersect(live_params(input$modality), names(p))]
  })
  observeEvent(live_rule_values(), {
    req(loaded_ok(), rv$slider_cfg)
    src <- tryCatch(isolate(th_state()$sources), error = function(e) character(0))
    new <- slider_config(rv$slider_cfg, rv$data$meta, auto_th())
    set <- rv$slider_set
    for (c in new) {
      cur <- isolate(input[[c$input_id]]) %||% c$value
      sides <- slider_sides(c)
      val <- cur
      for (side in names(sides)) if (!is_manual(src, c[[side]])) val[sides[[side]]] <- c$value[sides[[side]]]
      set[[c$id]] <- c$value
      updateSliderInput(session, c$input_id, value = val)
    }
    relabel_log_sliders(new)
    rv$slider_set <- set
  }, ignoreInit = TRUE)

  observeEvent(input$reset_auto, {
    cfg <- rv$slider_cfg; req(cfg)
    new <- slider_config(cfg, rv$data$meta, auto_th())
    for (c in new) updateSliderInput(session, c$input_id, value = c$value)
    relabel_log_sliders(new)
    rv$slider_set <- lapply(new, `[[`, "value")
    rv$slider_exact <- list()
  })

  ## ---- one sample's full state (on screen or not) ----------------------------------------
  # Used for saving: the same cutoff rule as the sliders, from the stored slider positions.
  sample_state <- function(id) {
    m <- isolate(rv$mod) %||% isolate(input$modality); env <- get_engine(m); e <- bank[[id]]
    p <- merge_live(e$p_loaded, isolate(current_params()$params))
    th <- auto_thresholds(m, env, e$data, p)
    if (identical(id, isolate(rv$active))) {
      cfg <- isolate(rv$slider_cfg); set <- isolate(rv$slider_set); exact <- isolate(rv$slider_exact)
      vals <- lapply(cfg, function(c) isolate(input[[c$input_id]]))
    } else if (!is.null(e$slider_cfg)) {
      cfg <- e$slider_cfg; set <- e$slider_set; vals <- e$slider_values; exact <- e$slider_exact
    } else {                 # never on screen: MAD cutoffs, or the ones its sample-sheet row fixes
      cfg <- slider_config(MODALITIES[[m]]$sliders, e$data$meta, th)
      set <- lapply(cfg, `[[`, "value"); vals <- set; exact <- NULL
      if (length(e$row_cutoffs)) { sc <- sheet_cutoff_sliders(cfg, e$row_cutoffs); vals <- lapply(sc$cfg, `[[`, "value"); exact <- sc$exact }
    }
    st <- resolve_cutoffs(th, cfg, vals, set, exact)
    meta <- suppressMessages(apply_filters(m, env, e$data, st$th, p))
    c(st, list(p = p, meta = meta, n_kept = sum(meta$keep), env = env, m = m))
  }

  ## ---- filtering (on-screen sample) ---------------------------------------------------------
  filtered_base <- reactive({
    if (!qc_ok()) { req(loaded_ok()); meta <- rv$data$meta; meta$qc_reason <- "pass"; meta$keep <- TRUE; return(meta) }
    st <- th_state()
    suppressMessages(apply_filters(rv$mod, engine(), rv$data, st$th, p_live()))
  })
  extra_valid <- reactive({
    if (!qc_ok()) return(FALSE)
    ex <- rv$extra
    !is.null(ex) && identical(ex$key, extra_key(filtered_base(), p_live()))
  })
  filtered <- reactive({
    meta <- filtered_base()
    if (isTRUE(mod()$extra_step) && extra_valid()) meta <- apply_contamination(rv$extra$meta, p_live())
    meta
  })

  ## ---- sample list (vertical tabs) -----------------------------------------------------------
  output$sample_column <- renderUI({
    s <- rv$samples
    if (is.null(s) || nrow(s) < 2 || !identical(input$input_mode, "dataset")) return(NULL)
    act <- rv$active
    mem <- in_memory()
    n_done <- sum(s$status %in% c("saved", "on_disk"))
    items <- lapply(seq_len(nrow(s)), function(i) {
      id <- s$sample_id[i]; st <- STATUS_STYLE[[s$status[i]]] %||% STATUS_STYLE$pending
      loaded <- id %in% mem
      meta <- if (!is.na(s$n_in[i]) && !is.na(s$n_kept[i])) sprintf("%s / %s kept", format(s$n_kept[i], big.mark = ","), format(s$n_in[i], big.mark = ","))
              else if (!is.na(s$n_in[i]) && loaded) sprintf("%s %s · QC not applied", format(s$n_in[i], big.mark = ","), mod()$unit) else st$label
      if (loaded && s$status[i] == "saved") meta <- paste(meta, "· saved")
      if (nzchar(s$patient[i] %||% "")) meta <- paste(meta, "·", s$patient[i])
      job <- s$job[i] %||% ""
      if (job %in% c("queued", "running")) meta <- paste(meta, "·", if (job == "running") paste0(rv$queue_label %||% "running", "…") else "queued")
      tags$button(type = "button", `data-sample` = id, title = paste0(st$label, if (nzchar(s$note[i] %||% "")) paste0(": ", s$note[i]) else ""),
                  class = paste("list-group-item list-group-item-action sample-pick d-flex gap-2 align-items-start",
                                if (identical(id, act)) "active", if (!loaded) "disabled",
                                if (job %in% c("queued", "running")) paste0("qc-job-", job)),
                  if (job %in% c("queued", "running")) span(class = "spinner-border spinner-border-sm text-primary flex-shrink-0 mt-1", role = "status")
                  else span(class = st$cls, icon(st$icon)),
                  div(class = "text-start", style = "min-width:0",
                      div(class = "text-truncate", style = "max-width: 180px", id),
                      div(class = "sample-meta small text-secondary", meta)))
    })
    div(class = "sample-col",
        card(
          card_header(div(class = "d-flex justify-content-between", span("Samples"),
                          span(class = "text-muted small", sprintf("%d / %d saved", n_done, nrow(s))))),
          if (!is.null(rv$queue)) div(class = "px-3 pt-2 small d-flex justify-content-between align-items-center",
              span(sprintf("%s: %d of %d done", rv$queue$label, rv$queue$done, rv$queue$total)),
              span(actionLink("queue_stop", "Stop after this one"), " · ",
                   actionLink("queue_terminate", "Terminate now", class = "text-danger"))),
          div(class = "list-group list-group-flush sample-list", items),
          card_footer(
            div(class = "form-text mb-2", rv$batch_info),
            busy_button("next_batch", "Load next batch", "Loading the next batch…", icon = icon("forward"), class = "btn-sm btn-outline-primary w-100 mb-2"),
            actionButton("save_rest", "Auto-QC the rest (MAD cutoffs)", icon = icon("wand-magic-sparkles"),
                         class = "btn-sm btn-outline-secondary w-100"),
            div(class = "form-text mt-1", "Loads, saves and unloads the samples not yet loaded, one at a time, with automatic cutoffs.")))
    )
  })

  output$dataset_table <- renderTable({
    s <- rv$samples; req(s)
    out <- data.frame(Sample = s$sample_id,
               Status = vapply(s$status, function(x) (STATUS_STYLE[[x]] %||% STATUS_STYLE$pending)$label, ""),
               `Looks like` = s$looks,
               In = ifelse(is.na(s$n_in), "", format(s$n_in, big.mark = ",")),
               Kept = ifelse(is.na(s$n_kept), "", format(s$n_kept, big.mark = ",")),
               `% kept` = ifelse(is.na(s$n_in), "", sprintf("%.1f", 100 * s$n_kept / pmax(s$n_in, 1))),
               `Manual cutoffs` = s$manual, Patient = s$patient, Annotated = s$annotated,
               Input = s$path, Note = s$note, check.names = FALSE)
    if (length(unique(s$dataset)) > 1) out <- cbind(out[1], Dataset = s$dataset, out[-1])
    out
  }, striped = TRUE, spacing = "s")

  ## ---- status + value boxes -------------------------------------------------------------------
  output$status <- renderUI({
    if (!loaded_ok()) {
      return(div(class = "alert alert-light status-banner",
                 icon("circle-info"), " Choose the data type, point to a sample or a dataset folder, and click ", strong("Load"), "."))
    }
    msgs <- list()
    lk <- rv$looks
    if (!isTRUE(rv$qc)) {
      msgs <- c(msgs, list(div(class = "alert alert-info status-banner py-2 d-flex align-items-center gap-3",
        span(icon("circle-info"), " ", strong(rv$active), paste0(" looks like ", lk$text, ". QC isn't applied: nothing is filtered yet.",
             if (identical(lk$kind, "raw")) " These are all barcodes, mostly empty droplets; Apply QC calls the cells first." else
             if (lk$kind %in% c("qc", "annotated")) " If it's already QC'd you may not need to filter it again." else "")),
        if (sample_job(rv$active) %in% c("queued", "running"))
          span(class = "ms-auto text-nowrap", span(class = "spinner-border spinner-border-sm me-1"),
               if (sample_job(rv$active) == "running") "Applying QC…" else "Queued")
        else busy_button("apply_qc_banner", apply_qc_label(), "Applying QC…", class = "btn-sm btn-primary ms-auto text-nowrap"))))
    } else if (isTRUE(lk$kind %in% c("qc", "annotated"))) {
      msgs <- c(msgs, list(div(class = "alert alert-warning status-banner py-2", icon("triangle-exclamation"),
        sprintf(" %s looks %s; QC cutoffs filter it a second time.", rv$active, lk$text))))
    }
    cur <- current_params()$params
    keys <- load_keys()
    if (identical(input$p_format %||% "auto", "auto")) keys <- setdiff(keys, "format")   # "auto" is worked out per sample
    changed <- Filter(function(nm) !identical(as.character(cur[[nm]]), as.character(rv$p_loaded[[nm]])), keys)
    if (length(changed)) msgs <- c(msgs, list(div(class = "alert alert-warning status-banner py-2",
      icon("triangle-exclamation"), sprintf(" Changed since loading: %s. Load again to apply.", paste(changed, collapse = ", ")))))
    if (isTRUE(mod()$extra_step) && !is.null(rv$extra) && !extra_valid()) msgs <- c(msgs, list(div(class = "alert alert-info status-banner py-2",
      icon("circle-info"), " The cell cutoffs or doublet settings changed after scDblFinder/DecontX ran; they'll be rerun when you save (or click Run again).")))
    tagList(msgs)
  })

  output$value_boxes <- renderUI({
    req(loaded_ok())
    meta <- filtered(); u <- mod()$unit
    n <- nrow(meta); k <- sum(meta$keep)
    if (!qc_ok()) return(layout_columns(
      fill = FALSE, col_widths = c(3, 3, 3, 3),
      value_box(title = paste(if (identical(rv$looks$kind, "raw")) "barcodes" else u, "loaded"), value = format(n, big.mark = ","), theme = "light", p(rv$active)),
      value_box(title = paste(u, "kept"), value = "—", theme = "primary", p("QC not applied")),
      value_box(title = "removed", value = "—", theme = "light", p("nothing filtered yet")),
      value_box(title = "input", value = rv$looks$short %||% "", theme = "light", p("as loaded"))))
    st <- th_state()
    layout_columns(
      fill = FALSE, col_widths = c(3, 3, 3, 3),
      value_box(title = paste(u, "in"), value = format(n, big.mark = ","), theme = "light", p(rv$active)),
      value_box(title = paste(u, "kept"), value = format(k, big.mark = ","), theme = "primary", p(sprintf("%.1f%%", 100 * k / max(n, 1)))),
      value_box(title = "removed", value = format(n - k, big.mark = ","), theme = "light",
                p(if (isTRUE(mod()$extra_step) && !extra_valid()) "before doublets / ambient" else "all filters")),
      value_box(title = "manual cutoffs", value = sum(st$sources == "manual"), theme = "light",
                p(if (any(st$sources == "manual")) paste(names(st$sources)[st$sources == "manual"], collapse = ", ") else "all from the MAD rule"))
    )
  })

  ## ---- plots ------------------------------------------------------------------------------------
  output$dist_plot <- renderPlot({
    req(loaded_ok())
    if (!qc_ok()) return(plot_distributions(filtered_base(), mod()$sliders, NULL, NULL,
                                            if (identical(rv$looks$kind, "raw")) "  (all barcodes; QC not applied)" else "  (QC not applied)"))
    st <- th_state()
    plot_distributions(filtered_base(), mod()$sliders, st$th, st$sources)
  }, height = function() 300 * ceiling(length(mod()$sliders) / 2), res = 96)

  SCATTER_MAX <- 60000
  output$scatter_plot <- renderPlot({
    req(loaded_ok())
    plot_scatter(filtered(), mod()$sliders, if (qc_ok()) th_state()$th, reason_colours(engine()), SCATTER_MAX)
  }, res = 96)
  output$scatter_note <- renderText({
    req(loaded_ok()); n <- nrow(rv$data$meta)
    if (n > SCATTER_MAX) sprintf("Showing a random %s of %s %s.", format(SCATTER_MAX, big.mark = ","), format(n, big.mark = ","), mod()$unit) else ""
  })

  SPATIAL_MAX <- 150000
  output$spatial_plot <- renderPlot({
    req(loaded_ok(), mod()$spatial)
    meta <- filtered()
    req(input$sp_colour %in% c("qc_reason", names(meta)))
    img <- if (isTRUE(input$sp_image)) rv$data$image else NULL
    plot_spatial(meta, input$sp_colour, reason_colours(engine()), img, SPATIAL_MAX,
                 input$sp_size * (if (rv$mod == "xenium") 0.3 else 1), isTRUE(input$sp_removed))
  }, res = 96)
  output$spatial_note <- renderText({
    req(loaded_ok(), mod()$spatial)
    n <- nrow(rv$data$meta)
    paste0(if (rv$mod == "visium" && is.null(rv$data$image)) "No tissue image found for this sample. " else "",
           if (n > SPATIAL_MAX) sprintf("Showing a random %s of %s %s.", format(SPATIAL_MAX, big.mark = ","), format(n, big.mark = ","), mod()$unit) else "")
  })

  ## ---- doublets / ambient (scRNA) ---------------------------------------------------------------
  run_extra_for <- function(env, d, meta, p, label) {
    if (!isTRUE(p$run_doublets) && !isTRUE(p$run_decontx)) return(run_extra_steps(env, d, meta, p))
    withProgress(message = sprintf("scDblFinder + DecontX: %s", label), detail = "this can take a few minutes", value = 0.3,
                 tryCatch(logged(run_extra_steps(env, d, meta, p)),
                          error = function(e) { fail(e, "Doublets / ambient RNA"); NULL }))
  }
  # scDblFinder + DecontX for one loaded sample, on the cells that pass its current cutoffs
  job_extra <- function(id) {
    st <- sample_state(id); e <- bank[[id]]
    ex <- run_extra_for(st$env, e$data, st$meta, st$p, id)
    if (is.null(ex)) return(FALSE)
    e$extra <- ex; bank[[id]] <- e
    if (identical(id, isolate(rv$active))) rv$extra <- ex
    TRUE
  }
  observeEvent(input$run_extra, {
    req(loaded_ok(), rv$mod == "scrna")
    if (!qc_ok() && !dataset_level_on()) { need_qc("scDblFinder + DecontX"); return() }
    if (dataset_level_on()) { stash_active(); start_queue("extra", "scDblFinder + DecontX"); return() }
    stash_active()
    job_extra(rv$active)
  })
  output$extra_plot <- renderPlot({
    req(loaded_ok(), rv$mod == "scrna")
    validate(need(qc_ok(), "Apply QC first: these run on the cells that pass the cutoffs."))
    validate(need(!is.null(rv$extra), "Not run yet. Click the button above (saving doesn't run it)."))
    th <- add_doublet_columns(th_state()$th, rv$extra$meta, rv$extra$doublet_threshold)
    g <- plot_extra(rv$extra$meta, th)
    validate(need(!is.null(g), "Both steps are switched off (run_doublets / run_decontx)."))
    g
  }, res = 96)

  ## ---- summary tables ---------------------------------------------------------------------------
  output$reason_table <- renderTable({
    req(loaded_ok())
    meta <- filtered()
    tb <- as.data.frame(table(reason = meta$qc_reason), stringsAsFactors = FALSE)
    tb <- tb[order(tb$reason != "pass", -tb$Freq), ]
    data.frame(Result = tb$reason, N = format(tb$Freq, big.mark = ","), `%` = sprintf("%.1f", 100 * tb$Freq / nrow(meta)), check.names = FALSE)
  }, striped = TRUE, spacing = "s")
  output$cutoff_table <- renderTable({
    req(loaded_ok())
    validate(need(qc_ok(), "QC not applied yet."))
    st <- th_state()
    nms <- names(st$sources)
    out <- data.frame(Cutoff = nms, `MAD rule` = fmt_num(unlist(st$auto[nms])), Used = fmt_num(unlist(st$th[nms])),
                      Source = unname(st$sources), check.names = FALSE)
    if (rv$mod %in% c("scrna", "visium")) out <- rbind(out, data.frame(Cutoff = "max_mito before clamping", `MAD rule` = fmt_num(st$auto$max_mito_mad),
                                                                        Used = "", Source = "", check.names = FALSE))
    out
  }, striped = TRUE, spacing = "s")
  output$median_table <- renderTable({
    req(loaded_ok())
    meta <- filtered()
    cols <- intersect(c(vapply(mod()$sliders, `[[`, "", "col"), "percent.ribo", "percent.hb", "nucleus_area", "decontX_contamination"), names(meta))
    data.frame(Metric = cols,
               `All` = vapply(cols, function(cl) fmt_num(median(meta[[cl]], na.rm = TRUE)), ""),
               `Kept` = vapply(cols, function(cl) fmt_num(median(meta[[cl]][meta$keep], na.rm = TRUE)), ""), check.names = FALSE)
  }, striped = TRUE, spacing = "s")
  output$log_text <- renderText(paste(rv$log, collapse = "\n"))

  ## ---- saving ---------------------------------------------------------------------------------------
  out_dir_ok <- function() {
    out_dir <- clean_path(input$out_dir)
    if (!nzchar(out_dir)) { showNotification("Choose an output folder first (Save panel).", type = "warning"); return(NULL) }
    out_dir
  }
  rds_path <- function(out_dir, dataset, id) file.path(out_dir, dataset, id, paste0(id, ".rds"))

  # A loaded sample's final QC result: cutoffs from its sliders, plus (scRNA) the scDblFinder / DecontX
  # results if they were computed for these cells (nothing is run here)
  qc_final <- function(id) {
    st <- tryCatch(sample_state(id), error = function(e) { fail(e, paste("Preparing", id)); NULL })
    if (is.null(st)) return(NULL)
    e <- bank[[id]]; p <- st$p; meta <- st$meta; th <- st$th; extra <- NULL
    if (isTRUE(MODALITIES[[st$m]]$extra_step)) {
      ex <- if (identical(id, isolate(rv$active))) isolate(rv$extra) else e$extra
      valid <- !is.null(ex) && identical(ex$key, extra_key(meta, p))
      if (valid) {
        e$extra <- ex; bank[[id]] <- e
        if (identical(id, isolate(rv$active))) rv$extra <- ex
        extra <- ex
        meta <- apply_contamination(ex$meta, p)
        th <- add_doublet_columns(th, meta, ex$doublet_threshold)
      } else {   # not run (or for other cells): the columns scrna_qc.R's save_sample() expects, left empty
        meta$scDblFinder.score <- NA_real_; meta$scDblFinder.class <- NA_character_; meta$decontX_contamination <- NA_real_
      }
    }
    meta$dataset <- p$dataset
    list(st = st, e = e, p = p, meta = meta, th = th, extra = extra, doublets_done = !is.null(extra))
  }

  # The settings a sample is annotated with: the Annotation panel's, with its sample-sheet row's
  # annot_* columns on top
  annot_params_for <- function(id) {
    base <- annot_params()
    ex <- if (exists(id, envir = bank, inherits = FALSE)) bank[[id]]$p_loaded$extra_columns
    if (!length(grep("^annot_", names(ex)))) return(base)
    rp <- isolate(aengine())$row_params(c(list(sample_id = id), ex), base)
    if (length(rp$problems)) add_log(paste(rp$problems, collapse = "\n"))
    rp$params
  }
  annot_row_text <- function(id) {
    ex <- if (exists(id, envir = bank, inherits = FALSE)) bank[[id]]$p_loaded$extra_columns
    ex <- ex[grep("^annot_", names(ex))]
    ex <- ex[!vapply(ex, function(v) is.na(v) || !nzchar(v), TRUE)]
    if (length(ex)) paste(sprintf("%s = %s", sub("^annot_", "", names(ex)), unlist(ex)), collapse = ", ")
  }

  # Save one loaded sample (on screen or not), as it is: QC, plus scDblFinder / DecontX and the
  # annotation when they were run on the cells being saved. Runs nothing. Returns export_sample()'s
  # result (with $left_out: what wasn't included and why), or NULL.
  save_one <- function(id, out_dir) {
    q <- qc_final(id)
    if (is.null(q)) return(NULL)
    st <- q$st; sm <- sample_meta_for(id)
    left_out <- character(0)
    if (isTRUE(MODALITIES[[st$m]]$extra_step) && !q$doublets_done && (isTRUE(q$p$run_doublets) || isTRUE(q$p$run_decontx)))
      left_out <- c(left_out, if (is.null(q$e$extra)) "scDblFinder + DecontX: not run" else "scDblFinder + DecontX: run on other cells (the cutoffs changed since); run them again")
    pre <- q$e$annot; kept <- rownames(q$meta)[q$meta$keep]
    annot <- NULL
    if (!is.null(pre) && identical(pre$key$cells, sort(kept))) annot <- list(r = pre, params = pre$key$params)
    else left_out <- c(left_out, if (is.null(pre)) "annotation: not run" else "annotation: made for other cells (the cutoffs or doublets changed since); annotate again")
    aenv <- if (!is.null(annot)) isolate(aengine())
    res <- withProgress(message = sprintf("Saving %s", id), value = 0.1,
      tryCatch(logged(export_sample(st$m, st$env, q$e$data, q$meta, q$th, q$p, out_dir, st$auto, st$sources, q$extra,
                                    progress = function(x) incProgress(0.3, detail = x), sample_meta = sm, annot = annot, aenv = aenv)),
               error = function(err) { fail(err, paste("Saving", id)); NULL }))
    if (is.null(res)) return(NULL)
    set_sample(id, status = "saved", n_in = res$n_in, n_kept = res$n_kept, manual = paste(res$manual, collapse = ", "), note = "",
               annotated = if (!is.null(annot)) "saved" else if (is.null(pre)) "" else "not saved")
    res$patient <- sm$patient
    res$left_out <- left_out
    res
  }

  # Annotate one loaded sample's passing cells (a preview: nothing is saved)
  job_annot <- function(id) {
    aenv <- isolate(aengine()); params <- annot_params_for(id)
    q <- qc_final(id)
    if (is.null(q)) return(FALSE)
    tmp <- preview_dir(id)
    dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
    r <- withProgress(message = sprintf("Annotating %s", id), detail = "SingleR (+ malignant-cell calling for scRNA-seq; can take a while)", value = 0.2,
      tryCatch(logged({
        rds <- q$st$env$save_sample(result_object(q$st$m, q$e$data, q$meta, q$th, q$p, q$extra, sample_meta_for(id)), file.path(tmp, "qc"))
        annotate_saved_sample(aenv, rds, id, q$p$dataset, sample_meta_for(id), params, work_dir = file.path(tmp, "cnv"))
      }), error = function(e) { fail(e, paste("Annotating", id)); NULL }))
    if (is.null(r)) return(FALSE)
    r$key <- annotation_key(rownames(q$meta)[q$meta$keep], params)
    r$doublets_done <- q$doublets_done
    e <- bank[[id]]; e$annot <- r; bank[[id]] <- e
    set_sample(id, annotated = "preview")
    rv$annot_n <- isolate(rv$annot_n) + 1
    TRUE
  }
  observeEvent(input$run_annot, {
    req(loaded_ok(), rv$active)
    if (!qc_ok() && !dataset_level_on()) { need_qc("annotation"); return() }
    aenv <- aengine()
    if (inherits(aenv, "error")) { fail(aenv, "Loading the annotation script"); return() }
    chk <- tryCatch({ aenv$check_packages(annot_params_for(rv$active)); NULL }, error = function(e) e)
    if (!is.null(chk)) { fail(explain_missing(conditionMessage(chk)), "Annotation"); return() }
    stash_active()
    if (dataset_level_on()) { start_queue("annot", "Annotating"); return() }
    # in a background process, so the app stays responsive (and connected) during a long copyKAT run
    if (HAS_CALLR) start_queue("annot", "Annotating", only = rv$active) else job_annot(rv$active)
  })

  ## ---- dataset-level runs: a queue, one sample per step ---------------------------------------
  # Between samples the app hands control back to the browser, so the sample list updates
  # and you can click samples (to look at earlier results) while the rest are still running.
  dataset_level_on <- function() isTRUE(input$dataset_level) && identical(input$input_mode, "dataset") && length(in_memory()) > 1
  JOBS <- list(qc = function(id) apply_qc(id), extra = function(id) job_extra(id), annot = function(id) job_annot(id))
  sample_job <- function(id) { s <- rv$samples; j <- s$job[match(id, s$sample_id)]; if (length(j) && !is.na(j)) j else "" }
  queue_banner <- function() {
    q <- isolate(rv$queue)
    if (is.null(q)) session$sendCustomMessage("qc-queue", list(on = FALSE))
    else session$sendCustomMessage("qc-queue", list(on = TRUE, text = sprintf("%s: %d of %d samples (now %s)…",
                                                                         q$label, q$done + 1, q$total, q$ids[1])))
  }
  # Run `task` on every sample in memory (or only on `only`, e.g. for the automatic pipeline).
  # Returns TRUE if a run started.
  start_queue <- function(task, label, only = NULL) {
    if (!is.null(isolate(rv$queue))) { showNotification("A dataset-level run is already going.", type = "warning"); return(FALSE) }
    s <- isolate(rv$samples)
    ids <- s$sample_id[s$sample_id %in% in_memory()]
    if (!is.null(only)) ids <- ids[ids %in% only]
    if (task == "qc") {                     # only the samples that still need it
      ids <- ids[!vapply(ids, has_qc, TRUE)]
      if (!length(ids)) { showNotification("QC is already applied to every sample in memory.", type = "message"); return(FALSE) }
    } else {
      skipped <- ids[!vapply(ids, has_qc, TRUE)]; ids <- setdiff(ids, skipped)
      if (length(skipped)) showNotification(sprintf("Skipping %d sample(s) without QC applied: %s", length(skipped), paste(skipped, collapse = ", ")),
                                            type = "warning", duration = 10)
      if (!length(ids)) return(FALSE)
    }
    active <- isolate(rv$active)
    ids <- c(intersect(active, ids), setdiff(ids, active))          # the sample on screen first
    for (id in ids) set_sample(id, job = "queued")
    set_sample(ids[1], job = "running")
    rv$queue_label <- label
    rv$queue <- list(task = task, label = label, ids = ids, total = length(ids), done = 0, ok = 0, ok_ids = character(0), started = FALSE)
    add_log(sprintf("---- %s%s on %d sample(s): %s ----", if (!is.null(isolate(rv$chain))) "automatic pipeline: " else "dataset-level ",
                    label, length(ids), paste(ids, collapse = ", ")))
    queue_banner()
    rv$queue_go <- isolate(rv$queue_go) + 1
    TRUE
  }
  # Background versions of job_extra / job_annot: prepare inputs here, run in a child R process
  start_bg_job <- function(task, id) {
    dir <- bg_job_dir(); out <- file.path(dir, "result.rds")
    if (task == "qc") {                     # cell calling (scRNA) in the child; cutoffs here when it's back
      m <- isolate(rv$mod); e <- bank[[id]]
      add_log(sprintf("---- applying QC to %s (background) ----", id))
      saveRDS(list(d = e$data, p = e$p_loaded), file.path(dir, "input.rds"), compress = FALSE)
      job <- bg_start(bg_qc_fun, list(app_dir = APP_DIR, qc_script = get_engine(m)$.script_path, modality = m,
                                      in_rds = file.path(dir, "input.rds"), out_rds = out), dir)
    } else if (task == "extra") {
      st <- sample_state(id); e <- bank[[id]]
      saveRDS(list(d = list(counts = e$data$counts, background = e$data$background), meta = st$meta, p = st$p),
              file.path(dir, "input.rds"), compress = FALSE)
      job <- bg_start(bg_extra_fun, list(app_dir = APP_DIR, qc_script = st$env$.script_path,
                                         in_rds = file.path(dir, "input.rds"), out_rds = out), dir)
    } else {
      if (bg_running_for(id)) stop(id, " is still being annotated by a run started before the page reconnected; ",
                                   "annotate again once it has finished (it then reuses that run's result)")
      aenv <- isolate(aengine()); params <- annot_params_for(id)
      q <- qc_final(id)
      if (is.null(q)) return(NULL)
      rds <- logged(q$st$env$save_sample(result_object(q$st$m, q$e$data, q$meta, q$th, q$p, q$extra, sample_meta_for(id)), file.path(dir, "qc")))
      job <- bg_start(bg_annot_fun, list(app_dir = APP_DIR, annot_script = aenv$.script_path,
                                         qc_rds = rds, sample_id = id, dataset = q$p$dataset,
                                         sample_meta = sample_meta_for(id), params = params,
                                         work_dir = file.path(preview_dir(id), "cnv"), checkpoint = labels_checkpoint(id),
                                         out_rds = out), dir)
      job$key <- annotation_key(rownames(q$meta)[q$meta$keep], params)
      job$doublets_done <- q$doublets_done
    }
    job$out <- out
    job
  }
  finish_bg_job <- function(task, id, job, res) {
    add_log(c(sprintf("---- %s: background log ----", id), res$log_lines))
    if (!res$ok) {
      fail(res$error, paste(switch(task, qc = "Applying QC to", extra = "scDblFinder + DecontX", "Annotating"), id)); return(FALSE)
    }
    if (!exists(id, envir = bank, inherits = FALSE)) return(FALSE)
    e <- bank[[id]]
    if (task == "qc") {
      ok <- finish_qc(id, readRDS(job$out))
      unlink(job$dir, recursive = TRUE)
      return(ok)
    }
    if (task == "extra") {
      ex <- readRDS(job$out); e$extra <- ex; bank[[id]] <- e
      if (identical(id, isolate(rv$active))) rv$extra <- ex
    } else {
      r <- readRDS(job$out); r$key <- job$key; r$doublets_done <- job$doublets_done
      e$annot <- r; bank[[id]] <- e
      set_sample(id, annotated = "preview")
      rv$annot_n <- isolate(rv$annot_n) + 1
    }
    unlink(job$dir, recursive = TRUE)
    TRUE
  }

  observe({
    rv$queue_go
    q <- isolate(rv$queue)
    if (is.null(q)) return()
    if (!q$started) { q$started <- TRUE; rv$queue <- q; invalidateLater(150); return() }   # let the list draw first
    id <- q$ids[1]
    # Applying QC to Visium / Xenium is only the cutoffs (instant), so it runs here
    if (HAS_CALLR && !(q$task == "qc" && !identical(isolate(rv$mod), "scrna"))) {
      if (is.null(q$job)) {                                   # start this sample's job
        job <- tryCatch(if (exists(id, envir = bank, inherits = FALSE)) start_bg_job(q$task, id),
                        error = function(e) { fail(e, paste(q$label, id)); NULL })
        if (!is.null(job)) { q$job <- job; rv$queue <- q; invalidateLater(500); return() }
        ok <- FALSE
      } else {
        res <- bg_poll(q$job)
        if (is.null(res)) { invalidateLater(500); return() }  # still running: check again soon
        ok <- isTRUE(tryCatch(finish_bg_job(q$task, id, q$job, res), error = function(e) { fail(e, paste(q$label, id)); FALSE }))
        q$job <- NULL
      }
    } else {
      ok <- isTRUE(tryCatch(if (exists(id, envir = bank, inherits = FALSE)) JOBS[[q$task]](id) else FALSE,
                            error = function(e) { fail(e, paste(q$label, id)); FALSE }))
    }
    set_sample(id, job = if (ok) "" else "failed", note = if (ok) "" else paste(q$label, "failed"))
    q$ids <- q$ids[-1]; q$done <- q$done + 1; q$ok <- q$ok + ok
    if (ok) q$ok_ids <- c(q$ok_ids, id)
    if (length(q$ids) && !isTRUE(q$stop)) {
      set_sample(q$ids[1], job = "running")
      rv$queue <- q; queue_banner(); invalidateLater(150)
    } else {
      for (rest in q$ids) set_sample(rest, job = "")
      rv$queue <- NULL; queue_banner()
      showNotification(sprintf("%s: %d of %d samples done%s.", q$label, q$ok, q$total,
                               if (length(q$ids)) sprintf(" (stopped; %d not run)", length(q$ids)) else ""),
                       type = if (q$ok < q$done) "warning" else "message", duration = 8)
      # automatic pipeline: the next step, on the samples this step succeeded on (none if stopped)
      if (!is.null(isolate(rv$chain))) {
        if (isTRUE(q$stop)) rv$chain <- NULL else run_chain(q$ok_ids)
        if (is.null(isolate(rv$chain)) && is.null(isolate(rv$queue)) && !isTRUE(q$stop) && length(q$ok_ids))
          showNotification("Automatic pipeline finished.", type = "message", duration = 8)
      }
    }
  })
  # When the page disconnects: an annotation job keeps running to the end (its labels and
  # malignant-cell results are cached), so a slow copyKAT run isn't lost; other jobs stop
  session$onSessionEnded(function() {
    q <- isolate(rv$queue)
    if (is.null(q$job)) return()
    if (q$task == "annot" && q$job$proc$is_alive()) bg_keep(q$job, q$ids[1], q$task, q$label)
    else bg_kill(q$job)
  })
  # A new page: say if an annotation started by an earlier page is still running, or has finished
  local({
    kept <- bg_kept()
    if (!length(kept)) return()
    msg <- lapply(kept, function(x) tags$li(sprintf("%s: %s (%s, started %s)", x$sample_id,
      if (x$alive) "still running" else sprintf("finished %s", format(x$ended, "%H:%M")),
      "annotation preview",
      format(x$started, "%H:%M"))))
    showNotification(duration = NULL, type = "message", tagList(
      strong("Annotation runs from before the page reconnected"), tags$ul(class = "mb-1 ps-3", msg),
      div(class = "small", "Annotating the same sample again with the same cells and settings reuses their labels and ",
          "malignant-cell results, so it's quick once they have finished. Then save it.")))
  })
  stop_queue <- function() {
    q <- rv$queue; req(q)
    q$stop <- TRUE; rv$queue <- q
    showNotification(if (!is.null(isolate(rv$chain))) "Stopping after the current sample; the automatic pipeline's later steps won't run."
                     else "Stopping after the current sample.", type = "message")
  }
  observeEvent(input$queue_stop, stop_queue())
  observeEvent(input$queue_stop_banner, stop_queue())

  # Terminate: kill the running background job now, drop the rest of the queue and of the
  # automatic pipeline, and remove the job's temporary files. Finished samples keep their results;
  # the sample being processed keeps whatever it had before this step (earlier steps' results,
  # the labels checkpoint, finished malignant-cell methods' caches).
  confirm_terminate <- function() {
    q <- rv$queue; req(q)
    showModal(modalDialog(
      title = "Terminate the current run?", easyClose = TRUE,
      p(sprintf("%s is running on %s. Terminating stops it now: that sample loses the progress of this step, and %s",
                q$label, q$ids[1],
                if (length(q$ids) > 1) sprintf("the %d sample(s) still waiting are not run.", length(q$ids) - 1) else "nothing else is waiting.")),
      if (length(isolate(rv$chain))) p("The automatic pipeline's later steps are cancelled too."),
      p(class = "text-muted small", "Samples already finished keep their results."),
      footer = tagList(modalButton("Keep running"),
                       actionButton("confirm_terminate", "Terminate now", class = "btn-danger"))))
  }
  observeEvent(input$queue_terminate, confirm_terminate())
  observeEvent(input$queue_terminate_banner, confirm_terminate())
  observeEvent(input$confirm_terminate, {
    removeModal()
    terminate_queue()
  })
  terminate_queue <- function() {
    q <- isolate(rv$queue)
    rv$chain <- NULL
    if (is.null(q)) return(invisible(FALSE))
    running <- q$ids[1]
    if (!is.null(q$job)) {
      bg_kill(q$job)
      unlink(q$job$dir, recursive = TRUE)
    }
    set_sample(running, job = "failed", note = paste(q$label, "terminated"))
    for (rest in q$ids[-1]) set_sample(rest, job = "")
    rv$queue <- NULL; queue_banner()
    add_log(sprintf("---- terminated: %s (%s; %d sample(s) not run) ----", q$label, running, length(q$ids) - 1))
    showNotification(sprintf("Terminated %s. %d of %d samples had finished.", q$label, q$ok, q$total), type = "warning", duration = 8)
    invisible(TRUE)
  }
  # Button labels say what they'll do
  observe({
    on <- isTRUE(input$dataset_level) && identical(input$input_mode, "dataset")
    updateActionButton(session, "run_extra", label = if (on) "Run scDblFinder + DecontX on all samples in memory" else "Run scDblFinder + DecontX on the cells that pass")
    updateActionButton(session, "run_annot", label = if (on) "Annotate all samples in memory (preview)" else "Annotate the cells that pass (preview)")
    updateActionButton(session, "apply_qc", label = if (on) "Apply QC to all samples in memory" else "Apply QC")
  })

  current_annot <- reactive({
    rv$annot_n; id <- rv$active
    if (is.null(id) || !exists(id, envir = bank, inherits = FALSE)) return(NULL)
    bank[[id]]$annot
  })
  output$annot_status <- renderUI({
    req(loaded_ok())
    r <- current_annot()
    aenv <- aengine()
    if (inherits(aenv, "error")) return(div(class = "alert alert-danger status-banner py-2", conditionMessage(aenv)))
    jb <- rv$samples$job[match(rv$active, rv$samples$sample_id)]
    if (isTRUE(jb %in% c("queued", "running"))) return(tagList(
      div(class = "alert alert-info status-banner py-2", span(class = "spinner-border spinner-border-sm me-2"),
          sprintf("This sample is %s in a background run (see the banner at the top).", if (jb == "running") "being processed" else "queued"),
          if (!is.null(r)) " Showing its previous annotation below."),
      if (!is.null(r)) div(class = "small mb-2", strong("Labels: "), annotation_summary(r$obj))))
    row_txt <- annot_row_text(rv$active)
    row_note <- if (!is.null(row_txt)) div(class = "form-text mb-2", "This sample's sample-sheet row sets ", code(row_txt),
                                           "; those override the Annotation panel for it.")
    if (is.null(r)) return(tagList(div(class = "alert alert-light status-banner py-2", icon("circle-info"),
                                       " Not annotated yet. Annotate it (button above) before saving to include the annotation."), row_note))
    kept <- sort(rownames(filtered())[filtered()$keep])
    cells_changed <- !identical(r$key$cells, kept)
    settings_changed <- !cells_changed && !identical(r$key$params, annot_params_for(rv$active))
    tagList(
      if (cells_changed) div(class = "alert alert-warning status-banner py-2", icon("triangle-exclamation"),
                             " Made for other cells: the cutoffs (or the doublet removal) changed since. Saving leaves it out; annotate again first."),
      if (settings_changed) div(class = "alert alert-info status-banner py-2", icon("circle-info"),
                                " The annotation settings changed since this was made. Saving keeps this result (and its settings); annotate again to use the new ones."),
      row_note,
      div(class = "small mb-2", strong("Labels: "), annotation_summary(r$obj)),
      if (!is.na(r$cnv$skipped %||% NA)) div(class = "form-text mb-2", "Malignant-cell step skipped: ", r$cnv$skipped))
  })
  output$annot_plot <- renderPlot({
    r <- current_annot(); req(r)
    plot_annotation_result(r)
  }, height = function() {
    r <- current_annot(); if (is.null(r)) return(1)
    if (!"umap" %in% Reductions(r$obj)) return(360)
    n <- length(annotation_panels(r))
    340 * ceiling(n / 2)
  }, res = 96)
  output$annot_table <- renderTable({ r <- current_annot(); req(r); annotation_table(r) }, striped = TRUE, spacing = "s")

  saved_modal <- function(results, ids) {
    ok <- Filter(Negate(is.null), results)
    if (!length(ok)) return()
    out_dir <- ok[[1]]$out_dir
    where <- p(class = "form-text", "Every setting used is in ", code("samplesheet.csv"), " (load it again, in the app or with the script, ",
               "to reproduce this; cutoffs you set by hand are its cutoff_* columns), the numbers per sample in ", code("summary.csv"),
               ", in ", code(out_dir), ".")
    if (length(ok) == 1) {
      r <- ok[[1]]; id <- names(ok)
      showModal(modalDialog(
        title = sprintf("Saved %s", id), easyClose = TRUE, size = "l",
        p(sprintf("%s of %s %s kept.", format(r$n_kept, big.mark = ","), format(r$n_in, big.mark = ","), mod()$unit)),
        p(strong("Patient: "), if (is.na(r$patient %||% NA)) span(class = "text-warning", "not set") else r$patient),
        p(strong("Changes from the QC defaults: "), r$deviations),
        if (!is.null(r$annotation)) p(strong("Annotation: "), r$annotation$summary),
        if (length(r$left_out)) div(class = "alert alert-light py-2", strong("Not included: "), tags$ul(class = "mb-0", lapply(r$left_out, tags$li))),
        where,
        tags$details(tags$summary(sprintf("Files (%s)", r$dir)),
                     tags$ul(lapply(unique(r$files[file.exists(r$files)]), function(f) tags$li(code(f)))))))
    } else {
      tab <- data.frame(Sample = names(ok), Patient = vapply(ok, function(r) if (is.na(r$patient %||% NA)) "—" else r$patient, ""),
                        Kept = vapply(ok, function(r) sprintf("%s / %s", format(r$n_kept, big.mark = ","), format(r$n_in, big.mark = ",")), ""),
                        `QC changes from defaults` = vapply(ok, `[[`, "", "deviations"),
                        Annotation = vapply(ok, function(r) if (is.null(r$annotation)) "—" else r$annotation$summary, ""),
                        `Not included` = vapply(ok, function(r) if (length(r$left_out)) paste(r$left_out, collapse = "; ") else "—", ""),
                        check.names = FALSE)
      showModal(modalDialog(title = sprintf("Saved %d samples", length(ok)), easyClose = TRUE, size = "l", html_table(tab), where))
    }
  }

  # Confirm before overwriting existing outputs, then run `action`
  confirm_then <- function(ids, action_name) {
    out_dir <- out_dir_ok(); if (is.null(out_dir)) return()
    s <- isolate(rv$samples)
    exist <- ids[file.exists(rds_path(out_dir, s$dataset[match(ids, s$sample_id)], ids))]
    if (length(exist)) {
      rv$pending_action <- action_name
      showModal(modalDialog(title = "Overwrite?",
        p(sprintf("%d sample(s) already have output in %s: %s", length(exist), out_dir,
                  paste(head(exist, 10), collapse = ", "), if (length(exist) > 10) ", …" else "")),
        p("Saving replaces their folders (an annotation saved there earlier too, unless it is saved again now) ",
          "and their rows in samplesheet.csv and summary.csv."),
        footer = tagList(modalButton("Cancel"), busy_button("confirm_action", "Overwrite", "Saving…", class = "btn-danger"))))
    } else run_action(action_name)
  }
  run_action <- function(action_name) {
    out_dir <- clean_path(input$out_dir)
    if (action_name == "save_active") {
      id <- rv$active
      saved_modal(setNames(list(save_one(id, out_dir)), id))
    } else if (action_name == "save_all") {
      stash_active()
      ids <- qc_ids()
      res <- setNames(lapply(ids, function(id) save_one(id, out_dir)), ids)
      saved_modal(res)
    } else if (action_name == "save_all_then_next") {
      stash_active()
      ids <- qc_ids()
      for (id in ids) save_one(id, out_dir)
      load_next_batch()
    } else if (action_name == "save_rest") {
      auto_qc_rest(out_dir)
    }
  }
  observeEvent(input$confirm_action, { removeModal(); act <- rv$pending_action; rv$pending_action <- NULL; run_action(act) })

  observeEvent(input$save, {
    req(loaded_ok(), rv$active)
    if (!qc_ok()) { need_qc("saving"); return() }
    confirm_then(rv$active, "save_active")
  })
  qc_ids <- function() Filter(has_qc, in_memory())
  observeEvent(input$save_all, {
    req(loaded_ok())
    ids <- qc_ids(); rest <- setdiff(in_memory(), ids)
    if (!length(ids)) { need_qc("saving"); return() }
    if (length(rest)) showNotification(sprintf("Saving the %d sample(s) with QC applied; skipping %s.", length(ids), paste(rest, collapse = ", ")),
                                       type = "warning", duration = 10)
    confirm_then(ids, "save_all")
  })

  observeEvent(input$next_batch, {
    req(rv$samples)
    s <- rv$samples
    unsaved <- intersect(in_memory(), s$sample_id[s$status == "loaded"])
    if (length(unsaved)) {
      showModal(modalDialog(title = "Unsaved samples",
        p(sprintf("%d sample(s) in memory haven't been saved: %s", length(unsaved), paste(unsaved, collapse = ", "))),
        p("Loading the next batch removes them from memory. Unsaved samples are loaded again after the rest."),
        footer = tagList(modalButton("Cancel"),
                         busy_button("next_discard", "Load next without saving", "Loading the next batch…", class = "btn-outline-danger"),
                         busy_button("next_save", "Save them, then load next", "Saving, then loading the next batch…", class = "btn-success"))))
    } else load_next_batch()
  })
  observeEvent(input$next_discard, { removeModal(); load_next_batch() })
  observeEvent(input$next_save, { removeModal(); confirm_then(qc_ids(), "save_all_then_next") })

  # Not-yet-loaded samples: load -> MAD cutoffs -> save -> unload, one at a time
  observeEvent(input$save_rest, {
    s <- rv$samples; req(s)
    rest <- s$sample_id[s$status %in% c("pending", "unsaved") & !s$sample_id %in% in_memory()]
    if (!length(rest)) { showNotification("Every sample is already loaded or saved.", type = "message"); return() }
    showModal(modalDialog(title = "Auto-QC the remaining samples?",
      p(sprintf("%d sample(s) will be loaded, QC applied with the MAD cutoffs and the current settings, saved, and unloaded, one at a time:", length(rest))),
      p(class = "small", paste(rest, collapse = ", ")),
      p("Samples in memory now are left alone. Only QC is saved for these (no scDblFinder + DecontX, no annotation); ",
        "each one's settings go into samplesheet.csv."),
      footer = tagList(modalButton("Cancel"), busy_button("save_rest_go", "Run", "Auto-QC running…", class = "btn-primary"))))
  })
  observeEvent(input$save_rest_go, {
    removeModal()
    s <- rv$samples
    confirm_then(s$sample_id[s$status %in% c("pending", "unsaved") & !s$sample_id %in% in_memory()], "save_rest")
  })
  auto_qc_rest <- function(out_dir) {
    stash_active()
    s <- isolate(rv$samples)
    rest <- s$sample_id[s$status %in% c("pending", "unsaved") & !s$sample_id %in% in_memory()]
    withProgress(message = "Auto-QC", value = 0, {
      for (k in seq_along(rest)) {
        id <- rest[k]
        incProgress(1 / length(rest), detail = sprintf("%s (%d of %d)", id, k, length(rest)))
        if (!load_one(id, s$path[match(id, s$sample_id)], s$dataset[match(id, s$sample_id)])) next
        if (!apply_qc(id)) next
        save_one(id, out_dir)
        rm(list = id, envir = bank); gc(verbose = FALSE)
      }
    })
    s <- isolate(rv$samples)
    showNotification(sprintf("Auto-QC finished: %d saved, %d failed (see the Dataset tab).",
                             sum(s$sample_id %in% rest & s$status == "saved"), sum(s$sample_id %in% rest & s$status == "failed")),
                     type = "message", duration = 10)
  }
}

shinyApp(ui, server)
