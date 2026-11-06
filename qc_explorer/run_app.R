# Start the QC app in your browser.
#   RStudio: open this file and click Source.
#   Terminal: Rscript run_app.R   (from this folder, or with the full path to this file)
# The app itself is in app_scripts/ next to this file.
here <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))      # Rscript run_app.R
  if (length(f)) return(dirname(normalizePath(f)))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {  # RStudio: Source
    p <- rstudioapi::getSourceEditorContext()$path
    if (nzchar(p)) return(dirname(normalizePath(p)))
  }
  getwd()
})
app_dir <- file.path(here, "app_scripts")
if (!file.exists(file.path(app_dir, "app.R")))
  stop("Can't find app_scripts/app.R in ", here, ". Open run_app.R from the qc_explorer folder (or setwd() there) and try again.")
options(shiny.maxRequestSize = 5 * 1024^3)
shiny::runApp(app_dir, launch.browser = TRUE)
