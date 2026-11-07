# geo_visium.R — Visium samples shipped on GEO as loose files instead of Space Ranger folders
# (e.g. GSE246011), converted into Space-Ranger-style folders the Visium QC script can read.
#
# Per sample, files sharing a prefix (GSM7853987_BLCA-B1 below):
#   GSM7853987_BLCA-B1_spot_coordinates.csv.gz   Space Ranger tissue_positions (v1 or v2 layout)
#   GSM7853987_BLCA-B1-count.csv.gz              dense counts (spots x genes or genes x spots)
#   GSM7853987_BLCA-B1.tif.gz                    H&E image (.tif / .png / .jpg, gzipped or not)
#   GSM7853987_BLCA-B1-Metadata.csv.gz           optional, not used
# There is no scalefactors_json.json, so the image scale (and any crop offset) is found by fitting
# the spots to the tissue: in-tissue spots should land on stained tissue, the other spots on glass.
# Check <sample>/spatial/alignment_check_app.png after converting.
#
# Each sample becomes <folder>/<prefix>/ with filtered_feature_bc_matrix/ and spatial/; the folder's
# converted samples are listed in visium_converted.csv (a sample sheet). Converting again is
# skipped while the folder is newer than its source files.

GEO_COORD_RX <- "[-_.]spot_coordinates\\.csv(\\.gz)?$"

# The loose-file Visium samples in a folder (only the one `only` belongs to, if given):
# data.frame(prefix, coords, counts, image)
geo_visium_groups <- function(dir, only = NULL) {
  none <- data.frame(prefix = character(0), coords = character(0), counts = character(0), image = character(0))
  if (!dir.exists(dir)) return(none)
  files <- list.files(dir, full.names = TRUE)
  files <- files[!dir.exists(files)]
  coords <- files[grepl(GEO_COORD_RX, basename(files), ignore.case = TRUE)]
  rows <- lapply(coords, function(cf) {
    pre <- sub(GEO_COORD_RX, "", basename(cf), ignore.case = TRUE)
    mine <- files[startsWith(basename(files), pre)]
    rest <- substring(basename(mine), nchar(pre) + 1)
    cnt <- mine[grepl("^[-_.]?(raw_?)?counts?([-_.]?matrix)?\\.(csv|tsv|txt)(\\.gz)?$", rest, ignore.case = TRUE)]
    img <- mine[grepl("^([-_.].*)?\\.(tiff?|png|jpe?g)(\\.gz)?$", rest, ignore.case = TRUE)]
    if (length(img) > 1) img <- img[which.max(file.size(img))]
    if (length(cnt) != 1 || length(img) != 1) return(NULL)
    data.frame(prefix = pre, coords = cf, counts = cnt, image = img)
  })
  out <- do.call(rbind, c(list(none), rows))
  if (!is.null(only)) {
    o <- normalizePath(only, winslash = "/")
    keep <- apply(out[, c("coords", "counts", "image"), drop = FALSE], 1, function(x) o %in% normalizePath(x, winslash = "/"))
    out <- out[keep, , drop = FALSE]
  }
  out
}

# A gzipped file unpacked to a temporary file (others returned as they are)
gunzip_tmp <- function(f) {
  if (!grepl("\\.gz$", f, ignore.case = TRUE)) return(f)
  out <- tempfile(fileext = paste0(".", tools::file_ext(sub("\\.gz$", "", f, ignore.case = TRUE))))
  i <- gzfile(f, "rb"); o <- file(out, "wb"); on.exit({ close(i); close(o) })
  repeat { b <- readBin(i, "raw", 1e7); if (!length(b)) break; writeBin(b, o) }
  out
}

# Dense counts -> genes x spots dgCMatrix. Orientation: whichever side has the barcodes.
read_geo_counts <- function(f) {
  dt <- data.table::fread(f, header = TRUE, data.table = TRUE, showProgress = FALSE)
  ids <- as.character(dt[[1]]); cols <- names(dt)[-1]
  bc <- function(x) mean(grepl("^[ACGTN]{10,}(-[0-9]+)?$", head(x, 200)))
  spots_in_rows <- bc(ids) >= bc(cols)
  blocks <- split(seq_along(cols), ceiling(seq_along(cols) / 1000))   # column blocks keep the dense copy small
  m <- do.call(cbind, lapply(blocks, function(ix) {
    x <- as.matrix(dt[, ix + 1L, with = FALSE]); storage.mode(x) <- "double"
    methods::as(Matrix::Matrix(x, sparse = TRUE), "CsparseMatrix")
  }))
  rm(dt)
  dimnames(m) <- list(ids, cols)
  if (spots_in_rows) m <- Matrix::t(m)
  rownames(m) <- make.unique(rownames(m))
  methods::as(m, "CsparseMatrix")
}

# tissue_positions with standard column names (v2 header, or v1 without one)
read_geo_positions <- function(f) {
  first <- readLines(if (grepl("\\.gz$", f)) gzfile(f) else f, n = 1)
  hdr <- grepl("barcode", first, ignore.case = TRUE)
  pos <- read.csv(f, header = hdr, stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(pos) < 6) stop(basename(f), " has ", ncol(pos), " columns; expected barcode, in_tissue, array_row, array_col, pxl_row, pxl_col")
  pos <- pos[, 1:6]
  names(pos) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
  pos
}

# Image as a [0, 1] height x width x 3 array
read_geo_image <- function(f) {
  src <- gunzip_tmp(f)
  ext <- tolower(tools::file_ext(src))
  need <- switch(ext, tif = , tiff = "tiff", png = "png", jpg = , jpeg = "jpeg", stop("unsupported image type: ", basename(f)))
  if (!requireNamespace(need, quietly = TRUE)) stop("reading ", basename(f), " needs the ", need, " package: install.packages(\"", need, "\")")
  img <- switch(need, tiff = tiff::readTIFF(src), png = png::readPNG(src), jpeg = jpeg::readJPEG(src))
  if (is.list(img)) img <- img[[1]]
  if (length(dim(img)) == 2) img <- array(img, c(dim(img), 3))
  if (dim(img)[3] == 4) img <- img[, , 1:3]
  if (dim(img)[3] == 2) img <- array(img[, , 1], c(dim(img)[1:2], 3))
  if (max(img) > 1) img <- img / 255
  img
}

# Every k-th pixel (k >= 1) so the longer side is at most `max_px`
shrink <- function(img, max_px) {
  k <- max(1L, as.integer(ceiling(max(dim(img)[1:2]) / max_px)))
  list(img = img[seq(1, dim(img)[1], by = k), seq(1, dim(img)[2], by = k), , drop = FALSE], k = k)
}

# Scale s and offset (dr, dc) mapping full-resolution spot pixels onto `img`
# (image row = s * pxl_row + dr), chosen so in-tissue spots land on tissue and the rest on glass.
# Tissue = saturated, not dark pixels (H&E pink/purple; glass, fiducials and background are grey).
fit_image_scale <- function(img, pos) {
  mx <- pmax(img[, , 1], img[, , 2], img[, , 3]); mn <- pmin(img[, , 1], img[, , 2], img[, , 3])
  tissue <- (mx - mn) / pmax(mx, 1e-6) > 0.15 & mx > 0.3
  H <- nrow(tissue); W <- ncol(tissue)
  pr <- pos$pxl_row_in_fullres; pc <- pos$pxl_col_in_fullres; inside <- pos$in_tissue == 1
  score <- function(s, dr, dc) {
    rr <- round(pr * s + dr) + 1; cc <- round(pc * s + dc) + 1
    ok <- rr >= 1 & rr <= H & cc >= 1 & cc <= W
    if (mean(ok) < 0.98) return(c(-1, NA, NA))
    v <- tissue[cbind(rr[ok], cc[ok])]; t <- inside[ok]
    on <- mean(v[t]); off <- if (any(!t)) mean(v[!t]) else 0
    c(on - off, on, off)
  }
  s_max <- min(H / max(pr), W / max(pc))
  grid <- function(ss, drs, dcs) {
    g <- expand.grid(s = ss, dr = drs, dc = dcs)
    sc <- t(mapply(score, g$s, g$dr, g$dc))
    cbind(g, score = sc[, 1], on = sc[, 2], off = sc[, 3])[which.max(sc[, 1]), ]
  }
  b <- grid(seq(0.5, 1.05, by = 0.005) * s_max, 0, 0)                                   # scale alone
  step <- round(0.01 * max(H, W))
  b <- grid(b$s * seq(0.97, 1.03, by = 0.005), seq(-7, 7) * step, seq(-7, 7) * step)    # + offset
  fs <- max(1, round(step / 5))
  b <- grid(b$s * seq(0.996, 1.004, by = 0.001), b$dr + seq(-5, 5) * fs, b$dc + seq(-5, 5) * fs)
  list(s = b$s, dr = b$dr, dc = b$dc, on = b$on, off = b$off, score = b$score, tissue = tissue)
}

# Picture of the fit: in-tissue spots green, others blue, on the image
alignment_png <- function(img, pos, fit, file) {
  sm <- shrink(img, 1200); k <- sm$k
  png(file, width = dim(sm$img)[2], height = dim(sm$img)[1])
  on.exit(dev.off())
  par(mar = c(0, 0, 0, 0))
  plot.new(); plot.window(c(0, dim(sm$img)[2]), c(dim(sm$img)[1], 0), xaxs = "i", yaxs = "i")
  rasterImage(sm$img, 0, dim(sm$img)[1], dim(sm$img)[2], 0)
  x <- (pos$pxl_col_in_fullres * fit$s + fit$dc) / k; y <- (pos$pxl_row_in_fullres * fit$s + fit$dr) / k
  points(x, y, pch = 1, cex = 0.35, col = ifelse(pos$in_tissue == 1, "#1a9e1a", "#3355ff"))
}

# Convert one sample (a row of geo_visium_groups()) into <its folder>/<prefix>/.
# Returns list(sample_id, path, reused, log).
convert_geo_visium <- function(g, progress = function(detail) NULL) {
  dir <- dirname(normalizePath(g$coords, winslash = "/"))
  id <- sanitize_id(g$prefix); sec <- file.path(dir, id)
  mdir <- file.path(sec, "filtered_feature_bc_matrix"); sdir <- file.path(sec, "spatial")
  done <- file.path(sdir, "scalefactors_json.json")
  if (file.exists(done) && file.mtime(done) >= max(file.mtime(c(g$coords, g$counts, g$image))))
    return(list(sample_id = id, path = id, reused = TRUE, log = sprintf("[%s] already converted (%s); reusing it", id, sec)))
  for (pk in c("data.table", "DropletUtils", "png"))
    if (!requireNamespace(pk, quietly = TRUE)) stop("converting ", id, " needs the ", pk, " package")

  progress("reading the counts")
  m <- read_geo_counts(g$counts)
  xs <- head(m@x, 2e6)
  if (length(xs) && (any(xs != round(xs)) || min(xs) < 0)) stop(basename(g$counts), " doesn't hold raw counts (non-whole numbers); QC needs raw counts")
  pos <- read_geo_positions(g$coords)
  miss <- setdiff(colnames(m), pos$barcode)
  if (length(miss)) stop(length(miss), " of ", ncol(m), " count-matrix barcodes are not in ", basename(g$coords), " (e.g. ", miss[1], ")")

  progress("reading the image and fitting it to the spots")
  img <- read_geo_image(g$image)
  hi <- shrink(img, 4000); rm(img)
  fit <- fit_image_scale(hi$img, pos)
  lo <- shrink(hi$img, 600)

  progress("writing the Space Ranger folder")
  dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
  DropletUtils::write10xCounts(mdir, m, barcodes = colnames(m), gene.id = rownames(m), gene.symbol = rownames(m),
                               gene.type = "Gene Expression", version = "3", overwrite = TRUE)
  # Space Ranger has no image offset: fold it into the full-resolution coordinates
  out <- pos
  out$pxl_row_in_fullres <- as.integer(round(pos$pxl_row_in_fullres + fit$dr / fit$s))
  out$pxl_col_in_fullres <- as.integer(round(pos$pxl_col_in_fullres + fit$dc / fit$s))
  write.csv(out, file.path(sdir, "tissue_positions.csv"), row.names = FALSE, quote = FALSE)
  # spot size from the spot spacing: centres 100 um apart (2 array columns), spots 55 um across
  key <- paste(pos$array_row, pos$array_col)
  nb <- match(paste(pos$array_row, pos$array_col + 2), key)
  spacing <- stats::median(sqrt((pos$pxl_row_in_fullres[nb] - pos$pxl_row_in_fullres)^2 +
                                (pos$pxl_col_in_fullres[nb] - pos$pxl_col_in_fullres)^2), na.rm = TRUE)
  sf <- c(spot_diameter_fullres = 0.55 * spacing, fiducial_diameter_fullres = 0.55 * spacing * 1.6,
          tissue_hires_scalef = fit$s, tissue_lowres_scalef = fit$s / lo$k)
  png::writePNG(hi$img, file.path(sdir, "tissue_hires_image.png"))
  png::writePNG(lo$img, file.path(sdir, "tissue_lowres_image.png"))
  alignment_png(hi$img, out, list(s = fit$s, dr = 0, dc = 0), file.path(sdir, "alignment_check_app.png"))
  writeLines(c("{", paste(sprintf('  "%s": %.10g', names(sf), sf), collapse = ",\n"), "}"), done)   # written last: marks the folder complete

  log <- sprintf("[%s] converted from %s, %s, %s -> %s: %d genes x %d spots; image fitted to the spots: scale %.5f, offset (%+.0f, %+.0f) px; %.0f%% of in-tissue spots on tissue, %.0f%% of the others",
                 id, basename(g$counts), basename(g$coords), basename(g$image), sec, nrow(m), ncol(m),
                 fit$s, fit$dr, fit$dc, 100 * fit$on, 100 * fit$off)
  warn <- fit$score < 0.3 || fit$on < 0.5
  if (warn) log <- c(log, sprintf("[%s] WARNING: the image fit is weak; check %s before trusting the spatial plots",
                                  id, file.path(sdir, "alignment_check_app.png")))
  list(sample_id = id, path = id, reused = FALSE, warn = warn, log = log)
}

# Convert every loose-file Visium sample in a folder (or the one `only` belongs to) and list them
# in <folder>/visium_converted.csv. Returns list(sheet, samples, log), or NULL if there were none.
convert_geo_visium_folder <- function(dir, only = NULL, progress = function(k, n, detail) NULL) {
  groups <- geo_visium_groups(dir, only)
  if (!nrow(groups)) return(NULL)
  dir <- normalizePath(dir, winslash = "/")
  dataset <- sanitize_id(sub("_RAW$", "", dataset_name_from_dir(dir), ignore.case = TRUE))
  log <- character(0); rows <- list()
  for (i in seq_len(nrow(groups))) {
    r <- convert_geo_visium(groups[i, ], progress = function(d) progress(i, nrow(groups), sprintf("%s: %s", groups$prefix[i], d)))
    log <- c(log, r$log)
    rows[[i]] <- data.frame(sample_id = r$sample_id, dataset = dataset, path = r$path,
                            notes = sprintf("converted from GEO files %s*", groups$prefix[i]), stringsAsFactors = FALSE)
  }
  samples <- do.call(rbind, rows)
  sheet <- file.path(dir, "visium_converted.csv")
  if (file.exists(sheet)) {   # keep samples converted earlier from this folder
    old <- tryCatch(read.csv(sheet, colClasses = "character"), error = function(e) NULL)
    if (!is.null(old)) samples <- rbind(old[!old$sample_id %in% samples$sample_id, names(samples), drop = FALSE], samples)
  }
  write.csv(samples, sheet, row.names = FALSE)
  list(sheet = sheet, samples = samples[samples$sample_id %in% sanitize_id(groups$prefix), , drop = FALSE], dataset = dataset, log = log)
}
