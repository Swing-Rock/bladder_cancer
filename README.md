# Bladder cancer project: analysis code

This repository holds the shared analysis code for the bladder cancer project, which uses single-cell RNA-seq (scRNA-seq), Xenium and Visium spatial transcriptomics data.

> **Work in progress.** Code is added as the project moves forward, so expect folders, scripts and conventions to change. Each folder has its own README describing what's there and how to run it.

## QC Explorer

[`qc_explorer/`](qc_explorer/) is the first step for every dataset: quality control and cell-type annotation, done the same way for everyone so results can be compared and merged across datasets. It has two parts that share the same code, so they always give the same results:

- **The QC app.** A point-and-click Shiny app that runs in your browser. Load a sample, a whole dataset folder or a sample sheet, look at the QC plots, set the cutoffs with sliders, and save. It records every cutoff you change from the defaults.
- **The QC scripts.** The same steps as R scripts driven by a sample sheet, for running without the app, for example on Expanse.

> **The app holds only a few samples at a time.** Everything loaded sits in your computer's memory, so the app loads up to 8 samples at once, and fewer if memory runs short. On a 16 GB laptop that's usually 4–8 scRNA-seq samples, about 8 Visium sections or 1–3 Xenium regions. For a larger dataset, either:
> - **work through it in batches in the app.** Save the samples on screen, then click *Load next batch*. *Auto-QC the rest* can finish the remaining samples with the default cutoffs once you trust them; or
> - **run the QC script** on the dataset's sample sheet, locally or on Expanse. It handles one sample at a time with no limit on how many. Copy any cutoffs you set by hand in the app into the sheet's `cutoff_*` columns.

| Step | scRNA-seq | Visium | Xenium |
| --- | --- | --- | --- |
| **QC**: per-sample filters on counts, genes and % mitochondrial reads | ✓ | ✓ | ✓ |
| Doublet removal (scDblFinder) and ambient-RNA estimation (DecontX) | ✓ | | |
| **Cell-type annotation**: SingleR and/or marker genes, per cell and per cluster | ✓ | ✓ (dominant type per spot) | ✓ |
| **Malignant-cell calls** from inferred copy-number changes (copyKAT / inferCNV) | ✓ | | |
| Cell-type proportions per spot (RCTD deconvolution) | | on hold | |

Every run writes the same output layout: one QC'd Seurat object per sample, a per-cell table, QC and annotation plots, and a `samplesheet.csv` that records every setting used. Give that sheet back to the app or the scripts to reproduce the run exactly.

**Get started:** read [`qc_explorer/README.md`](qc_explorer/README.md). It covers installing the packages, where to put your data and how to run the app or the scripts.

**Other docs in `qc_explorer/`:**

| Read | For |
| --- | --- |
| [`unconventional_datasets.md`](qc_explorer/DATASETS.md) | GEO datasets packaged unusually (GSE149652, GSE171351, GSE246011) and how to load them |
| [`app_scripts/README.md`](qc_explorer/app_scripts/README.md) | Using the app step by step |
| [`R_scripts/annotation/README.md`](qc_explorer/R_scripts/annotation/README.md) | Annotation settings, malignant-cell calls, and running copyKAT on Expanse |

never push data into the github repository. The folder name `datasets/` is included git ignores.

## Contributing

If you want editing access, message Shirley your GitHub username so she can add you as a contributor.