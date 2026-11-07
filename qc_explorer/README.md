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

## Project defaults

For this project, QC everyone's data with the defaults below, so samples can be compared and merged across datasets. They are the defaults in the scripts and the app, so you get them unless you change something. If you have to change one for a sample, it's recorded in `samplesheet.csv` and `log.txt`; mention it when you share results.

| Step | Default | Setting |
| --- | --- | --- |
| Cell-type annotation | **SingleR** | `label_methods = SingleR` |
| Annotation reference | **HPCA** only (Human Primary Cell Atlas, broad cell types). It fills `cell_label`, the label everything downstream uses | `references = HPCA`, `label_level = main` |
| Malignant cells (scRNA-seq only) | **inferCNV** | `malignant_method = infercnv` |

### When inferCNV can't run: use copyKAT

inferCNV compares each cell's inferred copy-number profile with **normal reference cells from the same sample**: the immune and stromal cells, picked by their `cell_label`. It needs at least **50** of them (`cnv_min_reference`). With fewer, it skips the sample, which happens with tumour-enriched or epithelial sorts and some small samples. You'll see:

- `log.txt`: `inferCNV skipped: only <n> reference (immune/stromal) cells`
- `summary.csv`: the `malignant` column says `skipped` instead of giving counts, and the sample's cells have no `cnv_status`.

**For those samples, run copyKAT instead.** It can find the normal (diploid) cells itself when there's no reference. Set `annot_malignant_method = copykat` on the sample's row in `samplesheet.csv` (or pick copyKAT in the app's Annotation panel) and annotate again. The labels are reused, so only copyKAT runs. copyKAT is slow: about an hour or more per sample, and much longer on Windows. For more than a sample or two, run it on Expanse (`R_scripts/annotation/README.md`, *Running copyKAT on Expanse*; list just those samples in `SAMPLES`).

Run copyKAT the same way if inferCNV fails with an error (`log.txt`: `WARNING: inferCNV failed`).

Other cases:

- **Samples with no tumour cells** (e.g. sorted immune cells such as GSE149652): skip the malignant-cell step altogether (`annot_malignant_method = none`). Neither method has anything to find. See [`Unconventional_datasets.md`](Unconventional_datasets.md).
- **Visium and Xenium**: there is no malignant-cell step. Spots mix several cells, and a targeted panel can't show copy-number changes.
- **inferCNV not installed** (it needs JAGS; see [Requirements](#requirements)): annotation stops before it starts and lists what to install. Install it rather than switching the whole dataset to copyKAT.

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

- in the Load data window that opens, choose the data type and point to one sample, to a dataset folder with one sub-folder per sample, or to a sample sheet (every row is loaded, each with its own settings). Loading filters nothing, so already-QC'd data can be looked at as it is
- click **Apply QC**, then adjust the cutoffs with the sliders
- optionally run scDblFinder + DecontX and the annotation preview
- choose an output folder and save. Saving writes what you have: it doesn't run any step itself

See [`app_scripts/README.md`](app_scripts/README.md) for the full walkthrough: dataset batches, patient IDs, annotation and dataset-level runs.

**If the whole page turns grey** and stops responding, the browser has lost its connection to R. The app runs inside R, in RStudio, and the page is only a window onto it. Common causes are the computer going to sleep, the browser putting the tab to sleep, R stopping after an error, or R running out of memory with too many large samples loaded. What to do:

1. **Look at the RStudio console.**
   - **The red stop sign is still showing** (R is still running the app): just reload the browser page.
   - **The console is back at the `>` prompt** (the app stopped): read the last lines for the error, then click **Source** on `run_app.R` again.
2. **Load your data again.** What was loaded lives in the page and is gone. Everything you **saved** is safe on disk. With *Skip samples already saved in the output folder* ticked, a dataset picks up where you left off.
3. **Annotation runs aren't lost** as long as R kept running. They finish in the background, and annotating the same sample again reuses their results instead of starting over.

To avoid it during long runs, keep the computer awake and stop the browser from putting the tab to sleep (Edge: *Settings → System and performance → Never put these sites to sleep*, add `127.0.0.1`; Chrome: *Settings → Performance → Always keep these sites active*). If it happens right after loading a big batch, R probably ran out of memory: lower *Max samples in memory at once* (in the *Advanced settings* panel, under *Dataset mode*; default 8) and load again.

Not the same thing: single **samples greyed out with a spinner** in the sample list are just queued for a background step. They become clickable when their turn is done.

### With the scripts

Write a sample sheet (start from the template in each script's folder), then run from the repo root:

```bash
Rscript R_scripts/scRNA/scrna_qc.R    datasets/GSE130001/samplesheet_GSE130001.csv
Rscript R_scripts/visium/visium_qc.R  datasets/<ID>/samplesheet_visium.csv
Rscript R_scripts/xenium/xenium_qc.R  datasets/<ID>/samplesheet_xenium.csv
Rscript R_scripts/annotation/annotate_cells.R datasets/<ID>/qc_output
```

Paths in a sample sheet are relative to the sheet's own folder, so the same sheet works on everyone's machine. Each QC script writes to `<sample sheet folder>/qc_output/` unless you give an output folder as the second argument. The annotation script writes into the QC output folder it is given.

## The output folder

The scripts and the app write the same layout:

```
qc_output/
├── samplesheet.csv        one row per sample: every setting used (QC, cutoffs set by hand, annotation)
├── summary.csv            one row per sample: cells in, removed per reason, kept, cutoffs, labels
├── log.txt                every run and every save adds to it
└── <dataset>/<sample>/
    ├── <sample>.rds          the QC'd Seurat object; the annotation adds its labels to it
    ├── <sample>_cells.csv    one row per cell, kept or not (Visium: <sample>_spots.csv)
    ├── qc/                   QC plots
    ├── doublets_ambient/     scDblFinder and DecontX plots (scRNA-seq, when they ran)
    └── annotation/           annotation plots (when the sample was annotated)
```

`samplesheet.csv` is also a sample sheet. Give it to the QC script, or load it in the app as a dataset, and every sample is run again with exactly the same settings, including the cutoffs set by hand in the app (its `cutoff_*` columns) and the annotation settings (its `annot_*` columns). Run on an output folder's own `samplesheet.csv`, a script writes its results back into that folder.

Any change to the QC parameters should be reported: `samplesheet.csv` holds every value used, and `log.txt` lists each sample's changes from the defaults in plain words.

Running QC on a sample again (script or app) replaces its folder, and its earlier annotation goes with it, because that annotation was made on the old cells. `working_files/` (copyKAT and inferCNV) only speeds up reruns and can be deleted.

## Documentation

| Read | For |
| --- | --- |
| [`app_scripts/README.md`](app_scripts/README.md) | Using the app, what it saves, its limits |
| [`R_scripts/scRNA/README.md`](R_scripts/scRNA/README.md) | scRNA-seq QC: method, sample-sheet columns, reading the plots |
| [`R_scripts/visium/README.md`](R_scripts/visium/README.md) | Visium QC and deconvolution |
| [`R_scripts/xenium/README.md`](R_scripts/xenium/README.md) | Xenium QC and its default filters |
| [`R_scripts/annotation/README.md`](R_scripts/annotation/README.md) | Annotation references, settings and the malignant-cell calls |
| [`Unconventional_datasets.md`](Unconventional_datasets.md) | GEO datasets packaged unusually (GSE149652, GSE171351, GSE246011) and how to load them |

## Requirements

R ≥ 4.3 with Seurat 5. The scripts were developed with R 4.3 and the app is used with R 4.5. Bioconductor packages (scDblFinder, celda, SingleR, celldex) and copyKAT / inferCNV are installed by `install_packages.R`.

Software to install outside R, before running `install_packages.R`:

| Software | Needed for | Notes |
| --- | --- | --- |
| [JAGS 4](https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/) | inferCNV (`malignant_method = "infercnv"`, the default, or `"copykat+infercnv"`) | inferCNV won't install or load without it, even though the annotation script doesn't run the step that uses it. On Windows, run the `JAGS-4.x.y.exe` installer with admin rights and the default location, then restart R. Check with `library(rjags)`, which should print "Linked to JAGS 4.x.y". |
| Python 3 | installing inferCNV | Any Python 3 (inferCNV's `argparse` dependency calls it). R must find it: on Windows, tick "Add python.exe to PATH" when installing, or point R to an existing one (e.g. Anaconda) with `PYTHON="C:/path/to/python.exe"` in `~/.Renviron`. The `WindowsApps` python stub doesn't work. |

Both are needed for the default malignant-cell method, inferCNV. Without them, choose `copykat` or `none` as `malignant_method`.
