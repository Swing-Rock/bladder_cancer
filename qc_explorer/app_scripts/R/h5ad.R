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
