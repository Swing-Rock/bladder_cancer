# h5ad.R — reading AnnData (.h5ad, from Scanpy / cellxgene / spatialdata-io) in the app.
#
# The QC scripts don't read .h5ad, so the app reads it here and hands the counts to the
# script's own functions, the same way as any other input. Uses hdf5r (already needed for
# 10x .h5 files); no Python needed.
#
# Counts come from the first of: layers/counts, raw/X, X. QC needs raw counts, so if the
# chosen matrix has non-integer values the log says so (it is probably normalized).

h5_attr <- function(obj, name) if (name %in% hdf5r::h5attr_names(obj)) hdf5r::h5attr(obj, name) else NULL

# A matrix stored in an .h5ad (sparse group or dense dataset), returned as genes x cells dgCMatrix
h5ad_matrix <- function(node) {
  if (inherits(node, "H5Group")) {
    enc <- h5_attr(node, "encoding-type") %||% if (!is.null(h5_attr(node, "h5sparse_format"))) paste0(h5_attr(node, "h5sparse_format"), "_matrix") else ""
    shape <- as.integer(h5_attr(node, "shape") %||% h5_attr(node, "h5sparse_shape"))   # n_obs, n_var
    x <- as.numeric(node[["data"]][]); i <- as.integer(node[["indices"]][]); p <- as.integer(node[["indptr"]][])
    if (enc == "csr_matrix") {
      # cells x genes in CSR == genes x cells in CSC
      return(Matrix::sparseMatrix(i = i + 1L, p = p, x = x, dims = c(shape[2], shape[1]), repr = "C"))
    }
    if (enc == "csc_matrix") {
      return(Matrix::t(Matrix::sparseMatrix(i = i + 1L, p = p, x = x, dims = c(shape[1], shape[2]), repr = "C")))
    }
    stop("unsupported sparse matrix encoding '", enc, "'")
  }
  m <- node[, ]                 # hdf5r reverses the dimensions: a dense (n_obs, n_var) X arrives as n_var x n_obs
  methods::as(Matrix::Matrix(m, sparse = TRUE), "CsparseMatrix")
}

# A column or index of a dataframe group (obs / var): plain arrays, strings, or categoricals
h5ad_column <- function(grp, name) {
  node <- grp[[name]]
  if (inherits(node, "H5Group")) {
    enc <- h5_attr(node, "encoding-type") %||% ""
    if (enc == "categorical" || node$exists("codes")) {          # codes + categories
      codes <- as.integer(node[["codes"]][]); cats <- node[["categories"]][]
      out <- rep(NA_character_, length(codes)); ok <- codes >= 0
      out[ok] <- as.character(cats[codes[ok] + 1L])
      return(out)
    }
    if (node$exists("values")) {                                 # nullable-string/integer/boolean: values + mask
      v <- node[["values"]][]
      if (node$exists("mask")) v[as.logical(node[["mask"]][])] <- NA
      return(v)
    }
    stop("can't read the column '", name, "' (encoding '", enc, "')")
  }
  v <- node[]
  cats <- h5_attr(node, "categories")                        # old-style categoricals
  if (!is.null(cats)) {
    lev <- tryCatch(grp$file[["/"]][[sub("^/", "", cats)]][], error = function(e) NULL)   # reference into uns/__categories
    if (!is.null(lev)) v <- as.character(lev[as.integer(v) + 1L])
  }
  v
}
h5ad_index <- function(grp) {
  idx <- h5_attr(grp, "_index") %||% "_index"
  if (!grp$exists(idx)) idx <- intersect(c("_index", "index", "obs_names", "var_names"), names(grp))[1]
  as.character(h5ad_column(grp, idx))
}
h5ad_columns <- function(grp) setdiff(names(grp), c(h5_attr(grp, "_index") %||% "_index", "__categories"))

# Gene names: gene_column = 2 wants symbols, 1 wants Ensembl IDs (as for 10x features.tsv)
h5ad_gene_names <- function(var, gene_column) {
  idx <- h5ad_index(var)
  cols <- h5ad_columns(var)
  ens <- mean(grepl("^ENS[A-Z]*G[0-9]", head(idx, 200))) > 0.5
  want_ens <- identical(as.character(gene_column), "1")
  if (want_ens && !ens) {
    alt <- intersect(c("gene_ids", "gene_id", "ensembl_id", "ensembl_ids", "feature_id"), cols)[1]
    if (!is.na(alt)) idx <- as.character(h5ad_column(var, alt))
  } else if (!want_ens && ens) {
    alt <- intersect(c("gene_symbols", "gene_symbol", "feature_name", "gene_name", "gene_names", "symbol", "features"), cols)[1]
    if (!is.na(alt)) idx <- as.character(h5ad_column(var, alt))
  }
  make.unique(idx)
}

# Read an .h5ad. Returns list(counts = genes x cells dgCMatrix, obs = data.frame, obs_cols,
# spatial = n x 2 matrix or NULL, source = where the counts came from, integer = TRUE/FALSE)
read_h5ad <- function(path, gene_column = 2, obs_wanted = character(0)) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) stop("reading .h5ad needs the hdf5r package: install.packages(\"hdf5r\")")
  f <- hdf5r::H5File$new(path, mode = "r")
  on.exit(f$close_all(), add = TRUE)
  if (f$exists("layers") && f[["layers"]]$exists("counts")) {
    m <- h5ad_matrix(f[["layers/counts"]]); src <- "layers/counts"; var <- f[["var"]]
  } else if (f$exists("raw") && f[["raw"]]$exists("X")) {
    m <- h5ad_matrix(f[["raw/X"]]); src <- "raw/X"; var <- f[["raw/var"]]
  } else {
    m <- h5ad_matrix(f[["X"]]); src <- "X"; var <- f[["var"]]
  }
  if (inherits(m, "dgCMatrix") == FALSE) m <- methods::as(m, "CsparseMatrix")
  obs <- f[["obs"]]
  cells <- h5ad_index(obs)
  if (ncol(m) != length(cells) && nrow(m) == length(cells)) m <- Matrix::t(m)
  rownames(m) <- h5ad_gene_names(var, gene_column)
  colnames(m) <- make.unique(cells)
  obs_cols <- h5ad_columns(obs)
  md <- data.frame(row.names = colnames(m))
  for (cl in intersect(obs_wanted, obs_cols)) md[[cl]] <- h5ad_column(obs, cl)
  spatial <- if (f$exists("obsm") && f[["obsm"]]$exists("spatial")) {
    s <- f[["obsm/spatial"]][, ]
    if (nrow(s) != ncol(m) && ncol(s) == ncol(m)) s <- t(s)
    s[, 1:2, drop = FALSE]
  }
  xs <- head(m@x, 1e5)
  list(counts = m, obs = md, obs_cols = obs_cols, spatial = spatial, source = src,
       integer = !length(xs) || all(xs == round(xs)))
}

# Log line about where the counts came from
h5ad_note <- function(h, id) {
  msg <- sprintf("[%s] .h5ad: counts from %s (%d genes x %d cells)", id, h$source, nrow(h$counts), ncol(h$counts))
  if (!h$integer) msg <- paste0(msg, "; WARNING: these are not whole numbers, so they look normalized. QC expects raw counts ",
                                "(Scanpy keeps them in layers['counts'] or .raw when saved that way)")
  msg
}

## ---- per data type: the same shapes the scripts' load_sample() returns ----------------------
# scRNA: a count matrix, cells named <sample>_<barcode>
h5ad_scrna <- function(p) {
  h <- read_h5ad(p$path, p$gene_column, obs_wanted = c(ANNOT_MARKER_COLS, QC_MARKER_COLS))
  message(h5ad_note(h, p$sample_id))
  m <- h$counts
  colnames(m) <- paste0(p$sample_id, "_", colnames(m))
  list(counts = m, obs_cols = h$obs_cols)
}

# Visium: a Seurat object with a Spatial assay; coordinates from obsm['spatial'] (no tissue image)
h5ad_visium <- function(p) {
  h <- read_h5ad(p$path, p$gene_column)
  message(h5ad_note(h, p$sample_id))
  if (is.null(h$spatial)) stop("the .h5ad has no obsm['spatial'] coordinates, so it can't be read as Visium")
  m <- h$counts
  colnames(m) <- paste0(p$sample_id, "_", colnames(m))
  obj <- SeuratObject::CreateSeuratObject(m, assay = "Spatial", project = p$sample_id)
  obj$pxl_col <- h$spatial[, 1]; obj$pxl_row <- h$spatial[, 2]   # what the script's spot_xy() reads first
  list(obj = obj, obs_cols = h$obs_cols)
}

# Xenium: the list the script's load_sample() returns (counts + per-cell info), from
# spatialdata-io-style obs columns
h5ad_xenium <- function(p) {
  h <- read_h5ad(p$path, p$gene_column, obs_wanted = c("cell_area", "nucleus_area", "x_centroid", "y_centroid",
    "control_probe_counts", "control_codeword_counts", "unassigned_codeword_counts", "deprecated_codeword_counts",
    "genomic_control_counts", "blank_codeword_counts"))
  message(h5ad_note(h, p$sample_id))
  m <- h$counts; o <- h$obs
  ctrl <- grep("control|codeword", names(o), value = TRUE)
  num <- function(cl) if (cl %in% names(o)) as.numeric(o[[cl]]) else rep(NA_real_, ncol(m))
  xy <- if (!is.null(h$spatial)) h$spatial else cbind(num("x_centroid"), num("y_centroid"))
  ids <- paste0(p$sample_id, "_", colnames(m))
  colnames(m) <- ids
  info <- data.frame(x = xy[, 1], y = xy[, 2],
                     control_counts = if (length(ctrl)) rowSums(sapply(o[ctrl], as.numeric)) else 0,
                     cell_area = num("cell_area"), nucleus_area = num("nucleus_area"), row.names = ids)
  if (!length(ctrl)) message(sprintf("[%s] no negative-control counts in the .h5ad; %% control will be 0", p$sample_id))
  list(counts = m, info = info, obj = NULL, obs_cols = h$obs_cols)
}

## ---- combined Visium .h5ad: several sections in one file (e.g. GSE171351) --------------------
# Some GEO series ship every Visium section in one AnnData file. The app splits such a file into
# one Space-Ranger-style folder per section, written next to the file, then loads the sections as
# a dataset. Each section keeps its tissue image and scale factors, and visium_qc.R can read the
# folders too. Same output as R_scripts/visium/h5ad_to_spaceranger.py, without needing Python.
# The list of sections goes in <file>_sections.csv (a sample sheet); while it exists, dataset
# folders skip the combined file itself (see find_samples()).

# Does an hdf5r file/group hold this path (e.g. "uns/spatial")?
h5_has <- function(node, p) {
  for (k in strsplit(p, "/", fixed = TRUE)[[1]]) { if (!node$exists(k)) return(FALSE); node <- node[[k]] }
  TRUE
}

# The libraries (sections) in uns['spatial'] of an .h5ad, or character(0)
h5ad_visium_libraries <- function(path) {
  if (!length(path) || !grepl("\\.h5ad$", path, ignore.case = TRUE) || !file.exists(path) ||
      !requireNamespace("hdf5r", quietly = TRUE)) return(character(0))
  tryCatch({
    f <- hdf5r::H5File$new(path, mode = "r"); on.exit(f$close_all(), add = TRUE)
    if (h5_has(f, "uns/spatial")) names(f[["uns/spatial"]]) else character(0)
  }, error = function(e) character(0))
}
is_combined_visium_h5ad <- function(path) length(h5ad_visium_libraries(path)) > 1

# An image stored in uns['spatial'][lib]['images'] as a [0, 1] array (height x width x 3) for png
h5ad_image <- function(ds) {
  a <- ds$read()
  d <- dim(a)
  if (length(d) == 3 && d[1] %in% 3:4 && !d[3] %in% 3:4) a <- aperm(a, c(3, 2, 1))   # hdf5r reverses the dimensions
  else if (length(d) == 2) a <- t(a)
  a <- array(as.numeric(a), dim = dim(a))
  if (max(a, na.rm = TRUE) > 1) a <- a / 255
  if (length(dim(a)) == 3 && dim(a)[3] == 4) a <- a[, , 1:3]
  a[is.na(a)] <- 0
  pmin(pmax(a, 0), 1)
}

# Split a combined Visium .h5ad. Returns list(sheet = the <file>_sections.csv sample sheet,
# samples = its rows, dataset, reused = TRUE if an earlier split was reused, log = lines for the log).
# progress(k, n, name) is called before each section.
split_visium_h5ad <- function(path, progress = function(k, n, name) NULL) {
  path <- normalizePath(path, winslash = "/")
  dir <- dirname(path); sheet <- split_manifest(path)
  dataset <- sanitize_id(sub("(_combined.*)?\\.h5ad$", "", basename(path), ignore.case = TRUE))

  # Split before (and the file hasn't changed since): reuse the folders
  if (file.exists(sheet) && file.mtime(sheet) >= file.mtime(path)) {
    old <- tryCatch(read.csv(sheet, colClasses = "character"), error = function(e) NULL)
    if (!is.null(old) && nrow(old) && all(file.exists(file.path(dir, old$path, "spatial", "tissue_positions_list.csv"))))
      return(list(sheet = sheet, samples = old, dataset = old$dataset[1], reused = TRUE,
                  log = sprintf("[%s] already split into %d sections (%s); reusing them", basename(path), nrow(old), sheet)))
  }
  for (pk in c("hdf5r", "DropletUtils", "png"))
    if (!requireNamespace(pk, quietly = TRUE)) stop("splitting a combined .h5ad needs the ", pk, " package")

  f <- hdf5r::H5File$new(path, mode = "r"); on.exit(f$close_all(), add = TRUE)
  libs <- if (h5_has(f, "uns/spatial")) names(f[["uns/spatial"]]) else character(0)
  if (length(libs) < 2) stop(basename(path), " holds ", length(libs), " section(s) in uns['spatial']; nothing to split")
  if (!h5_has(f, "obsm/spatial")) stop(basename(path), " has no obsm['spatial'] (spot coordinates)")

  # Raw counts: the first of layers/counts, raw/X, X that holds whole numbers
  m <- NULL
  for (cd in list(c("layers/counts", "var"), c("raw/X", "raw/var"), c("X", "var"))) {
    if (!h5_has(f, cd[1]) || !h5_has(f, cd[2])) next
    x <- h5ad_matrix(f[[cd[1]]])
    xs <- head(x@x, 2e6)
    if (!length(xs) || (all(xs == round(xs)) && min(xs) >= 0)) { m <- x; src <- cd[1]; var <- f[[cd[2]]]; break }
  }
  if (is.null(m)) stop("no raw counts in ", basename(path), ": layers['counts'], .raw and X all hold non-whole numbers (normalized data). QC needs raw counts.")
  obs <- f[["obs"]]
  cells <- h5ad_index(obs)
  if (ncol(m) != length(cells) && nrow(m) == length(cells)) m <- Matrix::t(m)
  if (ncol(m) != length(cells)) stop("the count matrix has ", ncol(m), " columns but obs has ", length(cells), " spots")
  gene_ids <- h5ad_gene_names(var, 1)   # Ensembl IDs if the file has them, else the index
  symbols  <- h5ad_gene_names(var, 2)   # gene symbols if the file has them, else the index

  # Which obs column says which section each spot belongs to: values = the uns['spatial'] keys
  ocols <- h5ad_columns(obs)
  read_col <- function(cl) tryCatch(h5ad_column(obs, cl), error = function(e) NULL)
  lib_col <- NA_character_; lib_of <- NULL
  for (cl in ocols) {
    v <- read_col(cl)
    if (is.null(v) || length(v) != length(cells)) next
    v <- as.character(v)
    if (!anyNA(v) && all(unique(v) %in% libs)) { lib_col <- cl; lib_of <- v; break }
  }
  if (is.null(lib_of)) stop("can't tell which obs column assigns spots to the sections (", paste(libs, collapse = ", "),
                            "). Split the file with R_scripts/visium/h5ad_to_spaceranger.py --sample-col <column> instead.")

  # Name each section after a patient/sample column when it has one value per section (as --name-col)
  name_col <- NA_character_; name_of <- NULL
  for (cl in intersect(c("Patient", "patient", "patient_id", "Patient_ID", "donor", "donor_id", "sample", "sample_id", "Sample"), ocols)) {
    v <- read_col(cl)
    if (is.null(v) || length(v) != length(cells)) next
    v <- as.character(v)
    if (all(tapply(v, lib_of, function(z) length(unique(z))) == 1)) { name_col <- cl; name_of <- v; break }
  }

  to_int <- function(cl, default) {
    v <- if (cl %in% ocols) read_col(cl)
    if (is.null(v) || length(v) != length(cells)) return(rep(as.integer(default), length(cells)))
    if (is.logical(v)) return(as.integer(v))
    x <- suppressWarnings(as.integer(round(as.numeric(v))))
    tf <- toupper(as.character(v)) %in% c("TRUE", "FALSE")
    x[tf] <- as.integer(toupper(as.character(v[tf])) == "TRUE")
    x[is.na(x)] <- as.integer(default)
    x
  }
  in_tissue <- to_int("in_tissue", 1); arow <- to_int("array_row", 0); acol <- to_int("array_col", 0)
  xy <- f[["obsm/spatial"]][, ]
  if (nrow(xy) != length(cells) && ncol(xy) == length(cells)) xy <- t(xy)   # column 1 = x (image column), 2 = y (image row)

  log <- sprintf("[%s] combined Visium .h5ad: %d sections (obs['%s']), counts from %s, %d genes x %d spots%s",
                 basename(path), length(libs), lib_col, src, nrow(m), ncol(m),
                 if (!is.na(name_col)) sprintf("; sections named <%s>_<library>", name_col) else "")
  rows <- list(); used <- character(0)
  for (k in seq_along(libs)) {
    lib <- libs[k]; sel <- which(lib_of == lib)
    if (!length(sel)) { log <- c(log, sprintf("  %s: no spots, skipped", lib)); next }
    name <- sanitize_id(if (!is.na(name_col)) paste0(name_of[sel[1]], "_", lib) else lib)
    name <- make.unique(c(used, name))[length(used) + 1]; used <- c(used, name)
    progress(k, length(libs), name)
    sec <- file.path(dir, name); mdir <- file.path(sec, "filtered_feature_bc_matrix"); sdir <- file.path(sec, "spatial")
    # barcodes: drop the suffixes added when the sections were concatenated (AAAC...-1-0-0 -> AAAC...-1)
    bcs <- sub("^([ACGTN]+-[0-9]+)(-[0-9]+)+$", "\\1", cells[sel])
    if (anyDuplicated(bcs)) bcs <- cells[sel]
    dir.create(sec, recursive = TRUE, showWarnings = FALSE)   # write10xCounts writes via a temporary folder inside it
    DropletUtils::write10xCounts(mdir, m[, sel, drop = FALSE], barcodes = bcs, gene.id = gene_ids, gene.symbol = symbols,
                                 gene.type = "Gene Expression", version = "3", overwrite = TRUE)
    dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
    # Space Ranger v1 positions (no header): barcode, in_tissue, array_row, array_col, pxl_row, pxl_col
    pos <- data.frame(bcs, in_tissue[sel], arow[sel], acol[sel], as.integer(round(xy[sel, 2])), as.integer(round(xy[sel, 1])))
    write.table(pos, file.path(sdir, "tissue_positions_list.csv"), sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
    lg <- f[["uns/spatial"]][[lib]]
    sf <- if (lg$exists("scalefactors")) { g <- lg[["scalefactors"]]; vapply(names(g), function(nm) as.numeric(g[[nm]]$read())[1], 0) } else numeric(0)
    if (!any(c("tissue_lowres_scalef", "tissue_hires_scalef") %in% names(sf)))
      stop("section ", lib, " has no tissue_lowres_scalef / tissue_hires_scalef in uns['spatial']['", lib, "']['scalefactors']")
    writeLines(c("{", paste(sprintf('  "%s": %.15g', names(sf), sf), collapse = ",\n"), "}"), file.path(sdir, "scalefactors_json.json"))
    kinds <- if (lg$exists("images")) intersect(c("lowres", "hires"), names(lg[["images"]])) else character(0)
    if (!length(kinds)) stop("section ", lib, " has no images in uns['spatial']['", lib, "']['images']")
    for (kind in kinds) png::writePNG(h5ad_image(lg[["images"]][[kind]]), file.path(sdir, sprintf("tissue_%s_image.png", kind)))
    rows[[length(rows) + 1]] <- data.frame(sample_id = name, dataset = dataset, path = name,
                                           notes = sprintf("from %s, library %s", basename(path), lib), stringsAsFactors = FALSE)
    log <- c(log, sprintf("  %s -> %s: %d spots, images: %s", lib, sec, length(sel), paste(kinds, collapse = ", ")))
  }
  samples <- do.call(rbind, rows)
  write.csv(samples, sheet, row.names = FALSE)
  log <- c(log, sprintf("  sample sheet: %s", sheet))
  list(sheet = sheet, samples = samples, dataset = dataset, reused = FALSE, log = log)
}
