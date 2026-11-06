# QC explorer: interactive QC for scRNA-seq, Visium and Xenium

A small Shiny app for tuning QC by eye. Load one sample, look at its QC plots, apply QC, drag sliders to move the cutoffs, and save the QC'd Seurat object to a folder you choose.

**The app does not have its own copy of the QC logic.** It loads the lab's scripts and calls their functions directly:

| Data type | Script used | Functions called |
| --- | --- | --- |
| scRNA-seq | `../R_scripts/scRNA/scrna_qc.R` | `resolve_sample`, `load_sample`, `call_cells`, `compute_qc_metrics`, `compute_thresholds`, `apply_cell_filters`, `detect_doublets`, `correct_ambient`, `save_sample`, `summarise_sample`, `make_plots` |
| Visium | `../R_scripts/visium/visium_qc.R` | same, with `apply_spot_filters` and `spot_xy` |
| Xenium | `../R_scripts/xenium/xenium_qc.R` (new, see below) | same structure |
| Annotation (all three) | `../R_scripts/annotation/annotate_cells.R` | `prepare_references`, `annotate_one_sample` (SingleR + copyKAT / inferCNV), `make_sample_plots` |

Each script is loaded into its own environment. Only the final "Run" block (`commandArgs()` / `main()` / `quit()`) is skipped. If you edit a script, the app reloads it the next time you click **Load**. The defaults, the help text under each setting, and the cutoff rules all come from the script.

## Running it

Both files are in the repo root (`qc_explorer/`), one folder up from here.

```r
# once, from the qc_explorer folder
source("install_packages.R")

# every time: open qc_explorer/run_app.R in RStudio and click Source
# (or open app_scripts/app.R and click "Run App")
```

It expects the repo layout below. Scripts in another location can be set with the `QCAPP_SCRIPTS_DIR` (the `R_scripts` folder) / `QCAPP_SCRNA_SCRIPT` / `QCAPP_VISIUM_SCRIPT` / `QCAPP_XENIUM_SCRIPT` / `QCAPP_ANNOTATION_SCRIPT` environment variables.

```
qc_explorer/
├── R_scripts/          the QC and annotation scripts (the app only reads these)
│   ├── scRNA/   visium/   xenium/   annotation/
└── app_scripts/        this app
```

## Using it

The app opens with the **Load data** window (reopen it with *Load data…* at the top of the sidebar).

1. **Choose the data type**: scRNA-seq, Visium or Xenium.
2. **Choose the input**:
   - **One sample**: paste a path or use *Folder… / File…*. You can also use **Fill in from a sample sheet**: pick an existing sample sheet and a row, and that row's path, sample ID, dataset and parameters are filled in. The script's `resolve_sample()` parses the row, so relative paths and blanks behave exactly as they do in the script.
   - **What can be read:**

     | Data type | Inputs |
     | --- | --- |
     | scRNA-seq | 10x folder or `.h5`, dense `.csv`/`.tsv`, `.rds` (Seurat, SingleCellExperiment or a count matrix), `.h5ad` |
     | Visium | Space Ranger `outs/` folder, Seurat `.rds` with a `Spatial` assay, `.h5ad` with `obsm['spatial']` |
     | Xenium | Xenium output folder, `LoadXenium()` `.rds`, `.h5ad` (e.g. from spatialdata-io: `cell_area`, `nucleus_area` and control counts in `obs`) |

     For `.h5ad`, the counts come from `layers['counts']` if it exists, else `.raw`, else `X`; the log warns if they aren't whole numbers (normalized data). `gene_column` picks symbols or Ensembl IDs from `var` as it does for 10x files.
   - **Dataset folder**: a folder with one sub-folder per sample (see *Dataset folders* below). The dataset name defaults to the folder name, or to its parent's name if the folder is called `data`, `raw`, `samples`, …
3. **Load.** This only runs the script's reading and metric steps: nothing is filtered, and for scRNA-seq no cells are called yet. The plots show the data as loaded, and a banner says what it looks like: a raw droplet matrix (mostly empty barcodes), cells or spots only, or an object that has already been QC'd (a `_qc.rds` / `_annotated.rds` file, or one with doublet / DecontX / SingleR columns). Already-QC'd data can be looked at without being filtered again.
4. **Apply QC** (in the *QC cutoffs* panel or the banner; *Apply QC to all samples in memory* in a dataset). For scRNA-seq this runs the script's cell calling (emptyDrops on a raw matrix; `raw_path` as DecontX background) and measures the called cells; then each sample gets its MAD cutoffs and the sliders appear. Saving, scDblFinder + DecontX and annotation need QC applied first.
   - **Automatic pipeline.** Turn on *Run the next steps automatically after Apply QC* (in the *Advanced settings* panel, under *Pipeline*) and Apply QC is followed by scDblFinder + DecontX (scRNA-seq) and then the annotation preview, in that order, on the cells that pass the MAD cutoffs. The steps run in the background; the banner at the top shows progress and has two links: *Stop after this one* finishes the current sample and cancels the rest, and *Terminate now* (after a confirmation) kills the running job straight away, cancels the rest and deletes the job's temporary files. Finished samples keep their results either way; a terminated sample keeps what it had before that step. Both links also appear above the sample list during dataset-level runs, and work for any background run, not only the automatic pipeline. With dataset-level processing, each step runs on every sample in memory before the next one starts, and a sample that fails a step is left out of the later ones. The annotation uses the *Annotation* panel's settings: pick `none` as the malignant-cell method for a quick pass, since copyKAT can take an hour per sample. If you move a slider afterwards, the results are out of date: run the steps again, or just save (saving redoes them).
5. **Adjust the cutoffs.**
   - Each slider starts at the cutoff the script would use (labelled *MAD*). If you leave a slider there, the exact MAD value is used.
   - Dragging a slider makes that cutoff *manual*. It's shown in orange with the MAD value beside it, and the plots update as you drag. Dragging a "max" handle all the way right means no maximum.
   - **Or type the value** in the *min* / *max* box under the slider and press Enter (or click elsewhere). The slider moves there, and the exact number you typed is used (the slider itself can only stop on its own steps, so its handle may show a nearby value). An empty *max* box means no maximum. The boxes always show the cutoffs in use. A min above the max is refused. Dragging the slider afterwards replaces the typed value; *Reset to MAD cutoffs* clears it.
   - *Reset to MAD cutoffs* puts every slider back.
   - **Advanced settings** under the sliders in *QC cutoffs* has the QC filtering parameters: cutoff rules (`nmads_*`, floors, caps, `mito_floor`/`mito_ceiling`, …), doublets / ambient RNA, the gene filter, and gene names / mito pattern / emptyDrops cell calling (these last ones apply when you load again and apply QC). The separate **Advanced settings** panel (below *Annotation* in the sidebar) has the general ones: the automatic pipeline switch, `format`, `is_filtered`, `raw_path`, the seed, and *Max samples in memory at once*. Changing a *cutoff rule* (`nmads_*`, floors, caps, `mito_floor`/`mito_ceiling`, …) recomputes the MAD cutoffs live, and any slider still on MAD follows. *Reading the data* settings (format, `gene_column`, `mito_pattern`, `is_filtered`, …) take effect when you load again; a yellow banner reminds you.
6. **Look at the plots.**
   - *Distributions*: one histogram per metric, removed zones shaded.
   - *Scatter*: UMIs vs genes and UMIs vs % mito, coloured by why each cell was removed.
   - *Spatial* (Visium/Xenium): where removed spots/cells sit. Visium spots are drawn on the H&E.
   - *Doublets & ambient* (scRNA): runs scDblFinder + DecontX on the cells that currently pass.
   - *Summary*: counts per removal reason, the cutoffs (MAD value vs value used), and medians.
   - *Dataset* (dataset mode): every sample's status, cells in/kept and manual cutoffs.
7. **Patient ID.** Under *Patient & sample info*, type the on-screen sample's patient ID. In dataset mode you can instead import a `.csv` with `sample_id` and `patient_id` columns; any other columns (tissue, stage, treatment, …) are imported too. These go into every cell's metadata (`patient`, …), into `qc_parameters_used.csv` and `all_samples_summary.csv`, and into the annotated object. A sample sheet loaded in one-sample mode can carry a `patient_id` column too.
8. **Annotation.** The *Annotation* panel holds the annotation script's settings, in sections: *Labelling methods* (SingleR and/or marker genes, and which one wins), *SingleR*, *Marker genes*, *Clustering*, *Malignant cells* and *Other settings*. Each section shows only its main settings; the rarely changed ones (custom reference, custom marker list and cutoffs, PCs / variable genes, the normal-reference regex, inferCNV's settings, cores, seed, …) are under its own **Advanced** link, which shows "(n changed)" when any of them differ from the script's defaults. Sections and settings for methods that aren't selected are hidden: choose SingleR only and the *Marker genes* section disappears; inferCNV's settings only appear when inferCNV is chosen; `label_priority` only when both labelling methods are. Labels are kept when only the malignant-cell settings change, so switching the method doesn't rerun SingleR. The *Annotation* tab previews the labels, clusters and malignant calls for the cells that currently pass. With *Annotate each sample when saving* ticked (the default), every save also writes `annotation/<dataset>/<sample>_annotated.rds`. If the cells and settings haven't changed since a preview, the preview's result is reused. Annotating runs in a **background process** (with the `callr` package), both the preview and the annotation done when saving, so the app stays responsive during a long copyKAT run, and the banner at the top shows it with *Stop after this one* / *Terminate now*. When saving, the annotation starts once the QC results are written, and the save dialog says it's running; its files appear in `annotation/` when it finishes. If an annotation package is missing, the app says which one and offers to save the QC results without annotating. See [`../R_scripts/annotation/README.md`](../R_scripts/annotation/README.md) for what each setting does.
9. **Save.** Choose an output folder, then click *Save this sample*, or *Save all samples in memory* in dataset mode. For scRNA, scDblFinder and DecontX run first if they haven't already for the current cutoffs. You're asked before anything is overwritten, and a summary of each sample's changes from the defaults is shown.

### Dataset folders

```
datasets/GSE.../data/          <- point the app here
├── patient1/                  each sub-folder = one sample; its name becomes the sample ID
│   └── filtered_feature_bc_matrix/   (or outs/filtered_feature_bc_matrix.h5, GEO-prefixed files, ...)
├── patient2/
├── GSM0000003_counts.h5       scRNA only: a count file per sample also works
└── patient4.h5ad              any data type: one .rds or .h5ad file per sample also works
```

- **Finding each sample's input.** For scRNA, the app picks the 10x folder or `.h5` in each sub-folder, preferring *filtered* over *raw*, then a single count file (including `.rds` and `.h5ad`). For Visium and Xenium, the sub-folder itself is the input, since those scripts search inside it; a sub-folder holding just one `.rds` or `.h5ad` uses that file. `.rds` and `.h5ad` files sitting directly in the dataset folder are one sample each (a `_qc` / `_annotated` ending is dropped from the sample ID). Sub-folders with nothing readable are skipped and listed in a notice.
- **The sample list.** Samples appear as a vertical list to the left of the plots. Click a loaded sample to view it. **Each sample keeps its own slider positions** and its own scDblFinder/DecontX results, while the *Advanced settings* apply to all of them.
- **How many samples are held in memory (x).** Samples are loaded one at a time until either:
  - *Max samples in memory at once* is reached (default 8; under *Advanced settings → Dataset mode*), or
  - loading another would leave less than 2 GB of free RAM. The next sample's size is estimated as 1.5× the largest sample loaded so far.

  On a 16 GB laptop this is usually 4–8 scRNA samples, around 8 Visium sections, or 1–3 Xenium regions. The note under the list says how many are loaded and why it stopped.
- **Dataset-level processing.** Switch it on under *Advanced settings → Dataset mode*. *Apply QC*, *Run scDblFinder + DecontX* and *Annotate (preview)* then run on every sample in memory instead of just the one on screen (*Apply QC* on the ones that don't have it yet; the other two on the ones that do):
  - Samples run one after another, each in a **separate background R process** (the `callr` package), so the app stays usable while they run. You can click any sample, including queued ones, to look at its plots and earlier results.
  - Waiting samples are greyed out with a spinner in the sample list. The banner at the top shows progress, and *Stop after this one* ends the run early.
  - Each sample uses its own cutoffs as they were when the run started.
  - Each background process loads R packages first, which adds roughly 10–20 seconds per sample.
  - Without `callr` installed, samples run inside the app instead, and clicks wait until the current sample is done.
  - Applying QC to Visium and Xenium is only the cutoffs, which is instant, so it runs inside the app.
- **Next batches.** *Load next batch* unloads the current samples and loads the next ones. If any aren't saved, you can save them first or skip them; skipped samples come back after the rest.
- **Auto-QC the rest.** This processes the samples that haven't been loaded, one at a time: load, apply QC with the MAD cutoffs and the current Advanced settings, save, unload. It's for when you've reviewed enough samples to trust the defaults.
- **Resuming.** With *Skip samples already saved in the output folder* ticked, samples that already have a `<sample>_qc.rds` are marked done and not loaded again, so a large dataset can be finished over several sessions.

## What gets saved

The layout is the scripts' layout, so `visium_deconvolution.R` and anything else reading `<dataset>/<sample>_qc.rds` works unchanged:

```
<output folder>/
├── all_samples_summary.csv     one row per sample (each save adds or replaces that sample's row)
├── qc_parameters_used.csv      every parameter + final cutoff per sample, and every deviation from the defaults
├── samplesheet_resolved.csv    every parameter value used, per sample (same columns as the scripts write)
├── app_samplesheet_rows.csv    sample-sheet rows that reproduce the app's cutoffs with the script
├── qc_log.txt                  the script's log messages, appended (never wiped)
├── annotation/                 if annotating: <dataset>/<sample>_annotated.rds, labels, plots,
│                               annotation_parameters_used.csv (every setting + deviations), annotation_log.txt
└── <dataset>/
    ├── <sample_id>_qc.rds      saved by the script's save_sample()
    ├── qc_thresholds.csv       cutoffs used + mad_<cutoff> columns + manual_cutoffs
    ├── <sample_id>_cell_qc_metadata.csv   (spot_… for Visium) every cell, kept or not, with qc_reason
    └── plots/<sample_id>/      the script's make_plots() figures for this sample
```

If the script can't draw a sample's plots (for example a Visium `.h5ad`, which has no tissue image), the app saves its own versions instead (`01_distributions_app.png`, `02_counts_vs_genes_app.png`, `03_spatial_qc_app.png`) and says so in `qc_log.txt`.

There are two differences from running the script, both because the app saves one sample at a time. Metadata tables and plot folders are written per sample instead of one per dataset. Summary tables are updated in place instead of rewritten.

### Documenting manual changes

The scRNA README asks that any change to QC parameters be reported. The app records them in four places:

- **`qc_parameters_used.csv`** has one row per sample. It holds every sample-sheet parameter, each final cutoff with its source (`cutoff_min_counts`, `cutoff_min_counts_source` = `MAD rule` or `manual`), the QC script used, and two summary columns:
  - `n_deviations_from_default`
  - `deviations_from_default`, in plain words, e.g. `nmads_mito: default 3 -> 2; min_counts cutoff: MAD rule 569 -> set by hand 1995`. Samples with no changes say `none (all script defaults)`.

  The defaults it compares against are read from the script's `DEFAULTS`. `format` is not counted, since `auto` just becomes the detected format.

- `qc_thresholds.csv` lists every cutoff used, next to the MAD value it replaced (`mad_min_counts`, …), with `manual_cutoffs` naming the ones you moved.
- `qc_log.txt` gets a line like `cutoffs min_counts=1585 (manual; MAD 548), …`.
- `app_samplesheet_rows.csv` gets a row you can paste into the dataset's sample sheet so the script reproduces the result:
  - With no manual changes, the row is just the parameters.
  - With manual changes, every cutoff is pinned to its final value (`nmads_* = 99` so the MAD rule never binds, plus hard floors/caps, and `mito_floor = mito_ceiling` for mito). The app then re-runs the script's `compute_thresholds()` on that row and writes `reproduces_app_cutoffs = TRUE` only if the cutoffs match exactly. The extra `notes` and `reproduces_app_cutoffs` columns are ignored by the script.

## Xenium: new code, defaults need review

There was no Xenium QC script, so `R_scripts/xenium/xenium_qc.R` was written to match the scRNA and Visium scripts. Its default filters are first choices that still need review. See [`../R_scripts/xenium/README.md`](../R_scripts/xenium/README.md) for what it reads, the metrics, the defaults and how to run it on its own.

## Files

| File | What it is |
| --- | --- |
| `app.R` | the interface: sliders, settings, plots, save |
| `R/engine.R` | loads the scripts and maps each data type onto their functions |
| `R/export.R` | writes outputs, updates summary tables, records deviations, builds the reproducible sample-sheet rows |
| `R/annotation.R` | runs the annotation script from the app, caches references, writes annotation outputs |
| `R/background.R` | runs dataset-level steps in background R processes (callr) |
| `R/dataset.R` | dataset-folder scanning and the memory rule for batches |
| `R/h5ad.R` | reads `.h5ad` (AnnData) files with `hdf5r`, for all three data types |
| `R/sliders.R` | slider ranges and the rule that turns slider positions into cutoffs |
| `R/plots_live.R` | the fast interactive plots (saved figures come from the scripts' `make_plots()`) |
| `../install_packages.R`, `../run_app.R` | setup and launcher (in the repo root) |

## Limits

- Advanced settings are shared by every sample in a dataset. To give one sample different *rules*, load it on its own in one-sample mode (or use the scripts with a sample sheet; `app_samplesheet_rows.csv` feeds that). Its *cutoffs* can always be set per sample with the sliders.
- For very large samples, the interactive scatter and spatial plots draw a random subset (60k and 150k points); the counts and cutoffs always use every cell. The saved plots use every cell, except Xenium spatial maps, which are capped at 300k.
- **If the page disconnects** (laptop asleep, a browser tab put to sleep, a closed tab): what you had loaded is gone, since it lives in the page's session, but an annotation that was running **keeps going to the end**. A save's annotation still writes its outputs; a preview's result is kept. The next page that opens says which annotations are still running or have finished. Load the sample again with the same cutoffs and settings, and annotating or saving it reuses that run's labels and malignant-cell results instead of computing them again (as long as R itself keeps running). Other background steps (QC, scDblFinder + DecontX) stop when the page disconnects.
- To avoid disconnects during long runs, keep the computer awake, and stop the browser from putting the tab to sleep: in Edge, *Settings → System and performance → Never put these sites to sleep*, add `127.0.0.1`; in Chrome, *Settings → Performance → Always keep these sites active*.
- Loading, and *Apply QC* on one sample, block the app until they finish. Large raw scRNA matrices (emptyDrops, at Apply QC) and DecontX can take minutes; with *Dataset-level processing* on, they run in the background instead.
- `.h5ad` files are read by the app, not the QC scripts (which don't read `.h5ad`). Results are the same, but to rerun a sample with a script, save its counts as `.rds` or 10x first; `app_samplesheet_rows.csv` notes this. Visium `.h5ad` files load without the tissue image (the spot coordinates come from `obsm['spatial']`).
