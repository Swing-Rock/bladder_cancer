# dataset.R — dataset-folder mode: find the samples in a folder, and decide how many
# fit in memory at once.

## =============================================================================
## Finding samples: one sub-folder (or, for scRNA, one count file) per sample
## =============================================================================
SCRNA_FILE_RX <- "\\.(h5|h5ad|rds|csv|tsv|txt)(\\.gz)?$"
SPATIAL_FILE_RX <- "\\.(h5ad|rds)$"

sanitize_id <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)

# A combined Visium .h5ad split into one folder per section (R/h5ad.R) lists its sections here
split_manifest <- function(path) sub("\\.h5ad$", "_sections.csv", path, ignore.case = TRUE)

# Dataset name from a dataset folder; skips generic names like datasets/GSE123/data
dataset_name_from_dir <- function(dir) {
  dir <- gsub("\\\\", "/", sub("[/\\\\]+$", "", dir))
  b <- basename(dir)
  if (tolower(b) %in% c("data", "raw", "raw_data", "samples", "input", "counts", "matrices", "outs")) b <- basename(dirname(dir))
  sanitize_id(b)
}

# A sample folder of a QC output folder: <sample>/<sample>.rds (next to <sample>_cells.csv, qc/, ...)
saved_object <- function(d) {
  f <- file.path(d, paste0(basename(d), ".rds"))
  if (file.exists(f)) list(path = f, note = "saved by the QC scripts / app") else NULL
}

# Pick the single file/folder in a sample folder that the scRNA script should read.
# Returns list(path, note) or list(path = NA, note = why not).
resolve_scrna_input <- function(d) {
  if (!is.null(r <- saved_object(d))) return(r)
  # 10x files directly in the folder (including GEO-style prefixed names)
  if (length(list.files(d, pattern = "matrix\\.mtx(\\.gz)?$", ignore.case = TRUE)))
    return(list(path = d, note = "10x folder"))
  pick <- function(f, what) {
    if (length(f) > 1 && any(grepl("filtered", f, ignore.case = TRUE))) f <- f[grepl("filtered", f, ignore.case = TRUE)]
    if (length(f) > 1) f <- f[nchar(f) == min(nchar(f))]
    if (length(f) == 1) list(path = f, note = what) else NULL
  }
  h5 <- list.files(d, pattern = "feature_bc_matrix\\.h5$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (!is.null(r <- pick(h5, "10x .h5"))) return(r)
  mtx <- list.files(d, pattern = "matrix\\.mtx(\\.gz)?$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (!is.null(r <- pick(dirname(mtx), "10x folder"))) return(r)
  files <- list.files(d, pattern = SCRNA_FILE_RX, full.names = TRUE, ignore.case = TRUE)
  files <- files[!grepl("barcodes|features|genes|readme|notes|metadata|annotation|_cells\\.csv$|_spots\\.csv$", basename(files), ignore.case = TRUE)]
  files <- files[file.size(files) > 50 * 1024]   # count tables are big; this skips notes and small tables
  if (length(files) == 1) return(list(path = files, note = "single file"))
  if (length(files) > 1) return(list(path = NA, note = sprintf("%d count files; put one per folder", length(files))))
  list(path = NA, note = "no count matrix found")
}

resolve_spatial_input <- function(d, modality) {
  if (!is.null(r <- saved_object(d))) return(r)
  rds <- list.files(d, pattern = SPATIAL_FILE_RX, full.names = TRUE, ignore.case = TRUE)
  marker <- if (modality == "visium") "tissue_positions(_list)?\\.csv(\\.gz)?$" else "(cell_feature_matrix(\\.h5)?|cells\\.csv(\\.gz)?|cells\\.parquet)$"
  hit <- c(list.files(d, pattern = marker, recursive = TRUE, ignore.case = TRUE),
           list.dirs(d, full.names = FALSE)[grepl(marker, basename(list.dirs(d, full.names = FALSE)), ignore.case = TRUE)])
  if (length(hit)) return(list(path = d, note = if (modality == "visium") "Space Ranger folder" else "Xenium folder"))
  if (length(rds) == 1) return(list(path = rds, note = if (grepl("h5ad$", rds, ignore.case = TRUE)) ".h5ad" else "Seurat .rds"))
  if (length(rds) > 1) return(list(path = NA, note = sprintf("%d .rds/.h5ad files; put one per folder", length(rds))))
  list(path = NA, note = if (modality == "visium") "no tissue_positions file found" else "no cell_feature_matrix / cells table found")
}

# Returns list(samples = data.frame(sample_id, path, kind), skipped = data.frame(name, reason)).
# Also reads a dataset folder of a QC output folder (<out>/<dataset>/<sample>/<sample>.rds).
find_samples <- function(dir, modality) {
  entries <- list.files(dir, full.names = TRUE, no.. = TRUE)
  entries <- entries[!grepl("^(qc_output|annotation|working_files|\\.)|^(samplesheet|summary)\\.csv$|^log\\.txt$", basename(entries))]
  dirs <- entries[dir.exists(entries)]
  rows <- list(); skipped <- list()
  for (d in dirs) {
    r <- if (modality == "scrna") resolve_scrna_input(d) else resolve_spatial_input(d, modality)
    if (is.na(r$path)) skipped[[length(skipped) + 1]] <- data.frame(name = basename(d), reason = r$note)
    else rows[[length(rows) + 1]] <- data.frame(sample_id = sanitize_id(basename(d)), path = normalizePath(r$path, winslash = "/"), kind = r$note)
  }
  # GEO often ships one file per sample in a single folder (count files; .rds / .h5ad for any data type)
  rx <- if (modality == "scrna") SCRNA_FILE_RX else SPATIAL_FILE_RX
  files <- entries[!dir.exists(entries) & grepl(rx, entries, ignore.case = TRUE)]
  files <- files[!grepl("readme|notes|metadata|samplesheet|annotation", basename(files), ignore.case = TRUE) & file.size(files) > 50 * 1024]
  files <- files[!(grepl("\\.h5ad$", files, ignore.case = TRUE) & file.exists(split_manifest(files)))]   # split into section folders already
  for (f in files) rows[[length(rows) + 1]] <- data.frame(
    sample_id = sanitize_id(sub("_(qc|annotated)$", "", sub(rx, "", basename(f), ignore.case = TRUE))),
    path = normalizePath(f, winslash = "/"), kind = "file")
  samples <- if (length(rows)) do.call(rbind, rows) else data.frame(sample_id = character(0), path = character(0), kind = character(0))
  dup <- duplicated(samples$sample_id)
  if (any(dup)) samples$sample_id[dup] <- make.unique(samples$sample_id)[dup]
  samples <- samples[order(samples$sample_id), , drop = FALSE]
  list(samples = samples,
       skipped = if (length(skipped)) do.call(rbind, skipped) else data.frame(name = character(0), reason = character(0)))
}

## =============================================================================
## A sample sheet: every row is a sample with its own settings
## =============================================================================
# Returns list(samples = data.frame(sample_id, path, kind, dataset), rows = the sheet's rows by
# sample_id, dir = the sheet's folder, problems). Paths are resolved by the QC script when each
# sample is loaded (relative paths are relative to the sheet's folder).
read_dataset_sheet <- function(file, default_dataset) {
  sh <- read.csv(file, colClasses = "character", na.strings = c("", "NA"), check.names = FALSE, strip.white = TRUE)
  names(sh) <- tolower(trimws(names(sh)))
  sh <- sh[rowSums(!is.na(sh)) > 0, , drop = FALSE]
  miss <- setdiff(c("sample_id", "path"), names(sh))
  if (length(miss)) stop("the sample sheet has no ", paste(miss, collapse = " or "), " column")
  probs <- character(0)
  dup <- unique(sh$sample_id[duplicated(sh$sample_id)])
  if (length(dup)) probs <- c(probs, sprintf("duplicate sample_id: %s", paste(dup, collapse = ", ")))
  if (anyNA(sh$sample_id)) probs <- c(probs, sprintf("%d row(s) without a sample_id", sum(is.na(sh$sample_id))))
  sh <- sh[!is.na(sh$sample_id) & !duplicated(sh$sample_id), , drop = FALSE]
  if (!"dataset" %in% names(sh)) sh$dataset <- NA_character_
  sh$dataset[is.na(sh$dataset)] <- default_dataset
  dir <- dirname(normalizePath(file, winslash = "/"))
  rows <- lapply(seq_len(nrow(sh)), function(i) sh[i, , drop = FALSE]); names(rows) <- sh$sample_id
  list(samples = data.frame(sample_id = sh$sample_id, path = ifelse(is.na(sh$path), "", sh$path), kind = "sample sheet row",
                            dataset = sh$dataset, stringsAsFactors = FALSE),
       rows = rows, dir = dir, problems = probs)
}

## =============================================================================
## Memory: how many samples to hold at once
## =============================================================================
# Rule used by the app ("x"): load samples one at a time into the batch until either
#   - MAX_BATCH samples are loaded (setting, default 8), or
#   - loading one more would leave less than RAM_RESERVE_GB of free RAM, predicting the
#     next sample's footprint as 1.5x the largest sample loaded so far.
# On a 16 GB laptop that's typically 4-8 scRNA samples, ~8 Visium sections, 1-3 Xenium regions.
RAM_RESERVE_GB <- 2

# Free physical RAM in GB, or NA if it can't be read on this OS
free_ram_gb <- function() {
  tryCatch({
    sys <- Sys.info()[["sysname"]]
    if (sys == "Windows") {
      out <- suppressWarnings(system2("powershell", c("-NoProfile", "-NonInteractive", "-Command",
               "(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory"), stdout = TRUE, stderr = FALSE))
      kb <- as.numeric(trimws(out[nzchar(trimws(out))][1]))
      if (is.na(kb)) {   # older Windows: wmic
        out <- suppressWarnings(system2("wmic", c("OS", "get", "FreePhysicalMemory", "/Value"), stdout = TRUE, stderr = FALSE))
        kb <- as.numeric(sub(".*=", "", grep("FreePhysicalMemory=", out, value = TRUE)[1]))
      }
      kb / 1024^2
    } else if (file.exists("/proc/meminfo")) {
      l <- readLines("/proc/meminfo")
      as.numeric(gsub("\\D", "", grep("^MemAvailable", l, value = TRUE))) / 1024^2
    } else if (sys == "Darwin") {
      vm <- system2("vm_stat", stdout = TRUE)
      page <- as.numeric(sub(".*page size of (\\d+).*", "\\1", vm[1]))
      pages <- function(k) as.numeric(gsub("\\D", "", grep(k, vm, value = TRUE)[1]))
      (pages("Pages free") + pages("Pages inactive")) * page / 1024^3
    } else NA_real_
  }, error = function(e) NA_real_)
}

# Memory R is using right now (GB)
r_mem_used_gb <- function() sum(gc(verbose = FALSE)[, 2]) / 1024
