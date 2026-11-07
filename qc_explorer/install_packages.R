# Packages for the QC app. Run once from the qc_explorer folder:  source("install_packages.R")
# (The QC scripts' own packages are listed in R_scripts/*/README.md;
#  this installs those too, so a fresh machine can run both.)
cran <- c("shiny", "bslib", "shinyFiles", "Seurat", "Matrix", "ggplot2", "patchwork", "data.table", "hdf5r",
          "pheatmap", "viridis", "scales", "remotes", "callr",
          "png", "tiff", "jpeg")   # last row: tissue images (GEO Visium .tif files)
bioc <- c("SingleCellExperiment", "scDblFinder", "celda", "DropletUtils",
          "SingleR", "celldex", "scran", "BiocParallel", "SummarizedExperiment")   # last row: annotation
optional <- c("arrow")   # only for Xenium folders that have cells.parquet but no cells.csv.gz

need <- cran[!vapply(cran, requireNamespace, TRUE, quietly = TRUE)]
if (length(need)) install.packages(need)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
need <- bioc[!vapply(bioc, requireNamespace, TRUE, quietly = TRUE)]
if (length(need)) BiocManager::install(need, update = FALSE, ask = FALSE)
# Malignant-cell calling (annotation, scRNA-seq)
if (!requireNamespace("copykat", quietly = TRUE)) remotes::install_github("navinlabcode/copykat")
if (!requireNamespace("infercnv", quietly = TRUE))   # needs JAGS installed first: https://mcmc-jags.sourceforge.io
  tryCatch(BiocManager::install("infercnv", update = FALSE, ask = FALSE),
           error = function(e) message("inferCNV NOT installed (", conditionMessage(e), "). It is the default malignant-cell method: ",
                                     "install JAGS (see README > Requirements) and run this file again"))
message("Done. Optional: install.packages(\"arrow\") for Xenium cells.parquet.")
