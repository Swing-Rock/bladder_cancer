# QC Explorer

Quality control for **scRNA-seq**, **Visium** and **Xenium** data, written for the Ongkeko Lab bladder cancer project. It has two parts:

- **The QC scripts** (`R_scripts/`): sample-sheet-driven QC, one script per data type, plus cell-type annotation. They run on their own, e.g. on Expanse.
- **The QC app** (`app_scripts/`): a Shiny app for setting cutoffs by eye. It loads a sample or a whole dataset folder, shows the QC plots, lets you move the cutoffs with sliders, and saves the results. The app runs the scripts' own functions and has no separate copy of the QC logic, so both always give the same results.

> **Work in progress.** Folders, scripts and defaults may still change. Each folder has its own README with the details.

## What's here

```
qc_explorer/
├── run_app.R                 starts the QC app
├── install_packages.R        installs the app's and the scripts' packages (run once)
├── app_scripts/              the QC app (Shiny)
└── R_scripts/
    ├── scRNA/                scrna_qc.R: MAD filters, cell calling, scDblFinder, DecontX
    ├── visium/               visium_qc.R: spot QC
    │                         visium_deconvolution.R: cell-type proportions per spot (RCTD)
    │                         h5ad_to_spaceranger.py: split combined .h5ad files into per-section folders
    ├── xenium/               xenium_qc.R: cell QC for Xenium
    └── annotation/           annotate_cells.R: SingleR labels + malignant calls (copyKAT / inferCNV)
```

## Status

| Step | Data | What it does | Status | Location |
| --- | --- | --- | --- | --- |
| QC | scRNA-seq | Per-sample MAD filters, emptyDrops cell calling for raw matrices, doublet removal (scDblFinder), ambient-RNA estimation (DecontX) | Available | [`R_scripts/scRNA/`](R_scripts/scRNA/) |
| QC | Visium | Per-sample MAD filters on spots, maps of removed spots | Available | [`R_scripts/visium/`](R_scripts/visium/) |
| QC | Xenium | Per-sample MAD filters, cell-area window, negative-control cap | New; default filters need review | [`R_scripts/xenium/`](R_scripts/xenium/) |
| Annotation | all three | SingleR labels per cell and per cluster; malignant vs normal from inferred CNVs (scRNA-seq only) | New; the CNV step has not yet been run end to end on real data | [`R_scripts/annotation/`](R_scripts/annotation/) |
| Deconvolution | Visium | Cell-type proportions per spot (RCTD) | On hold, waiting on an annotated bladder scRNA-seq reference | [`R_scripts/visium/`](R_scripts/visium/) |
| QC app | all three | Sliders for cutoffs, dataset folders in batches, annotation preview, records every change from the defaults | Available | [`app_scripts/`](app_scripts/) |

## Getting started

1. **Clone the repo** and open it in RStudio, with the working directory set to the repo root. The scripts' default paths assume this.
2. **Install the packages** once:

   ```r
   source("install_packages.R")   # the app's packages and the scripts'
   ```

3. **Put your data** under `datasets/<dataset ID>/`, one folder per sample. `datasets/` is in `.gitignore`, so data never ends up in the repo.
4. **Run QC**, either in the app or with the scripts.

### With the app

Open `run_app.R` (in the repo root) in RStudio and click **Source**. Then:

- in the Load data window that opens, choose the data type and point to one sample, or to a dataset folder with one sub-folder per sample (loading filters nothing, so already-QC'd data can be looked at as it is)
- click **Apply QC**, then adjust the cutoffs with the sliders
- choose an output folder and save

See [`app_scripts/README.md`](app_scripts/README.md) for the full walkthrough: dataset batches, patient IDs, annotation and dataset-level runs.

### With the scripts

Write a sample sheet (start from the template in each script's folder), then run from the repo root:

```bash
Rscript R_scripts/scRNA/scrna_qc.R    datasets/GSE130001/samplesheet_GSE130001.csv
Rscript R_scripts/visium/visium_qc.R  datasets/<ID>/samplesheet_visium.csv
Rscript R_scripts/xenium/xenium_qc.R  datasets/<ID>/samplesheet_xenium.csv
Rscript R_scripts/annotation/annotate_cells.R datasets/<ID>/qc_output
```

Paths in a sample sheet are relative to the sheet's own folder, so the same sheet works on everyone's machine. Each script writes to `<sample sheet folder>/qc_output/` unless you give an output folder as the second argument.

## Recording QC decisions

Any change to the QC parameters should be reported. Both the scripts and the app write the parameters they used:

- **The scripts** write `samplesheet_resolved.csv`, the final value of every parameter for every sample.
- **The app** also writes:
  - `qc_parameters_used.csv` and `annotation/annotation_parameters_used.csv`, which list every deviation from the defaults in plain words, including cutoffs set by hand;
  - `app_samplesheet_rows.csv`, with sample-sheet rows that make the scripts reproduce the app's exact cutoffs.

## Documentation

| Read | For |
| --- | --- |
| [`app_scripts/README.md`](app_scripts/README.md) | Using the app, what it saves, its limits |
| [`R_scripts/scRNA/README.md`](R_scripts/scRNA/README.md) | scRNA-seq QC: method, sample-sheet columns, reading the plots |
| [`R_scripts/visium/README.md`](R_scripts/visium/README.md) | Visium QC and deconvolution |
| [`R_scripts/xenium/README.md`](R_scripts/xenium/README.md) | Xenium QC and its default filters |
| [`R_scripts/annotation/README.md`](R_scripts/annotation/README.md) | Annotation references, settings and the malignant-cell calls |

## Requirements

R ≥ 4.3 with Seurat 5. The scripts were developed with R 4.3 and the app is used with R 4.5. Bioconductor packages (scDblFinder, celda, SingleR, celldex) and copyKAT / inferCNV are installed by `install_packages.R`.

Software to install outside R, before running `install_packages.R`:

| Software | Needed for | Notes |
| --- | --- | --- |
| [JAGS 4](https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/) | inferCNV (`malignant_method = "infercnv"` or `"copykat+infercnv"`) | inferCNV won't install or load without it, even though the annotation script doesn't run the step that uses it. On Windows, run the `JAGS-4.x.y.exe` installer with admin rights and the default location, then restart R. Check with `library(rjags)`, which should print "Linked to JAGS 4.x.y". |
| Python 3 | installing inferCNV | Any Python 3 (inferCNV's `argparse` dependency calls it). R must find it: on Windows, tick "Add python.exe to PATH" when installing, or point R to an existing one (e.g. Anaconda) with `PYTHON="C:/path/to/python.exe"` in `~/.Renviron`. The `WindowsApps` python stub doesn't work. |

Neither is needed if you don't use inferCNV: copyKAT alone works without them.
