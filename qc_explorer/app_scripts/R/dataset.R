# dataset.R — dataset-folder mode: find the samples in a folder, and decide how many
# fit in memory at once.

## =============================================================================
## Finding samples: one sub-folder (or, for scRNA, one count file) per sample
## =============================================================================
SCRNA_FILE_RX <- "\\.(h5|h5ad|rds|csv|tsv|txt)(\\.gz)?$"
SPATIAL_FILE_RX <- "\\.(h5ad|rds)$"

sanitize_id <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)

# Dataset name from a dataset folder; skips generic names like datasets/GSE123/data
dataset_name_from_dir <- function(dir) {
  dir <- gsub("\\\\", "/", sub("[/\\\\]+$", "", dir))
  b <- basename(dir)
  if (tolower(b) %in% c("data", "raw", "raw_data", "samples", "input", "counts", "matrices", "outs")) b <- basename(dirname(dir))
  sanitize_id(b)
}

# Pick the single file/folder in a sample folder that the scRNA script should read.
# Returns list(path, note) or list(path = NA, note = why not).
resolve_scrna_input <- function(d) {
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
  files <- files[!grepl("barcodes|features|genes|readme|notes|metadata|annotation", basename(files), ignore.case = TRUE)]
  files <- files[file.size(files) > 50 * 1024]   # count tables are big; this skips notes and small tables
  if (length(files) == 1) return(list(path = files, note = "single file"))
  if (length(files) > 1) return(list(path = NA, note = sprintf("%d count files; put one per folder", length(files))))
  list(path = NA, note = "no count matrix found")
}

resolve_spatial_input <- function(d, modality) {
  rds <- list.files(d, pattern = SPATIAL_FILE_RX, full.names = TRUE, ignore.case = TRUE)
  marker <- if (modality == "visium") "tissue_positions(_list)?\\.csv(\\.gz)?$" else "(cell_feature_matrix(\\.h5)?|cells\\.csv(\\.gz)?|cells\\.parquet)$"
  hit <- c(list.files(d, pattern = marker, recursive = TRUE, ignore.case = TRUE),
           list.dirs(d, full.names = FALSE)[grepl(marker, basename(list.dirs(d, full.names = FALSE)), ignore.case = TRUE)])
  if (length(hit)) return(list(path = d, note = if (modality == "visium") "Space Ranger folder" else "Xenium folder"))
  if (length(rds) == 1) return(list(path = rds, note = if (grepl("h5ad$", rds, ignore.case = TRUE)) ".h5ad" else "Seurat .rds"))
  if (length(rds) > 1) return(list(path = NA, note = sprintf("%d .rds/.h5ad files; put one per folder", length(rds))))
  list(path = NA, note = if (modality == "visium") "no tissue_positions file found" else "no cell_feature_matrix / cells table found")
}

# Returns list(samples = data.frame(sample_id, path, kind), skipped = data.frame(name, reason))
find_samples <- function(dir, modality) {
  entries <- list.files(dir, full.names = TRUE, no.. = TRUE)
  entries <- entries[!grepl("^(qc_output|annotation|\\.)", basename(entries))]
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
