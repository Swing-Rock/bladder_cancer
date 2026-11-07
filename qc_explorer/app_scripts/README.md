# QC explorer: interactive QC for scRNA-seq, Visium and Xenium

A small Shiny app for tuning QC by eye. Load one sample, look at its QC plots, apply QC, drag sliders to move the cutoffs, and save the QC'd Seurat object to a folder you choose.

**The app does not have its own copy of the QC logic.** It loads the lab's scripts and calls their functions directly:

| Data type | Script used | Functions called |
| --- | --- | --- |
| scRNA-seq | `../R_scripts/scRNA/scrna_qc.R` | `resolve_sample`, `load_sample`, `call_cells`, `compute_qc_metrics`, `compute_thresholds`, `apply_cell_filters`, `detect_doublets`, `correct_ambient`, `save_sample`, `summarise_sample`, `make_plots`, `write_sample_outputs` (saving) |
| Visium | `../R_scripts/visium/visium_qc.R` | same, with `apply_spot_filters` and `spot_xy` |
| Xenium | `../R_scripts/xenium/xenium_qc.R` (new, see below) | same structure |
| Annotation (all three) | `../R_scripts/annotation/annotate_cells.R` | `prepare_references`, `annotate_one_sample` (SingleR + copyKAT / inferCNV), `row_params` (a sample sheet's `annot_*` columns), `write_annotation_outputs` (saving) |

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

     **Combined Visium `.h5ad`** (several sections in one file, e.g. GSE171351's `GSE171351_combined_visium.h5ad`): the app spots it (more than one library in `uns['spatial']`) and splits it into one Space-Ranger-style folder per section, written next to the file, with each section's tissue image and scale factors. Then it loads the sections as a dataset. This works whether you pick the file itself (in either mode) or its folder as a dataset. The sections are listed in `<file>_sections.csv`, a sample sheet that `visium_qc.R` can also read. Sections are named `<patient>_<library>` when `obs` has a patient/sample column with one value per section, otherwise by library. Loading the file again reuses the folders unless the `.h5ad` has changed. While `_sections.csv` exists, a dataset folder scan skips the combined file itself. Delete it to split again.

     **GEO Visium as loose files** (e.g. GSE246011: `<GSM>_<name>-count.csv.gz`, `<GSM>_<name>_spot_coordinates.csv.gz`, `<GSM>_<name>.tif.gz` per sample, no Space Ranger folder): pick any of a sample's files, or the folder holding them (e.g. the unpacked `GSE246011_RAW/`). The app converts each sample into a Space-Ranger-style folder `<GSM>_<name>/` next to the files, lists them in `visium_converted.csv`, and loads them. The counts may have spots or genes in rows (whichever side has the barcodes). These series don't include `scalefactors_json.json`, so the app works out the image scale, and any crop offset, by fitting the spots to the tissue: in-tissue spots should land on stained tissue, the others on glass. **Check `spatial/alignment_check_app.png`** in each sample folder (green = in-tissue spots, blue = the rest). The log gives the fit, and a weak fit is flagged. Reading `.tif` images needs the `tiff` package (`install_packages.R` installs it). Converting takes about a minute per sample, and later loads reuse the folders.
   - **Dataset**: either a folder with one sub-folder per sample (see *Dataset folders* below; the dataset name defaults to the folder name, or to its parent's name if the folder is called `data`, `raw`, `samples`, …), or a **sample sheet** `.csv`, whose rows are all loaded, each with its own settings (see *A dataset from a sample sheet* below). *Fill in from a sample sheet* is for one sample only, so it is hidden in dataset mode.
3. **Load.** This only runs the script's reading and metric steps: nothing is filtered, and for scRNA-seq no cells are called yet. The plots show the data as loaded, and a banner says what it looks like: a raw droplet matrix (mostly empty barcodes), cells or spots only, or an object that has already been QC'd (a `<sample>.rds` saved by the scripts or the app, an older `_qc.rds` / `_annotated.rds` file, or one with doublet / DecontX / SingleR columns). Already-QC'd data can be looked at without being filtered again.
4. **Apply QC** (in the *QC cutoffs* panel or the banner; *Apply QC to all samples in memory* in a dataset). For scRNA-seq this runs the script's cell calling (emptyDrops on a raw matrix; `raw_path` as DecontX background) and measures the called cells; then each sample gets its MAD cutoffs and the sliders appear. Saving, scDblFinder + DecontX and annotation need QC applied first.
   - **Automatic pipeline.** Turn on *Run the next steps automatically after Apply QC* (in the *Advanced settings* panel, under *Pipeline*) and Apply QC is followed by scDblFinder + DecontX (scRNA-seq) and then the annotation preview, in that order, on the cells that pass the MAD cutoffs. The steps run in the background; the banner at the top shows progress and has two links: *Stop after this one* finishes the current sample and cancels the rest, and *Terminate now* (after a confirmation) kills the running job straight away, cancels the rest and deletes the job's temporary files. Finished samples keep their results either way; a terminated sample keeps what it had before that step. Both links also appear above the sample list during dataset-level runs, and work for any background run, not only the automatic pipeline. With dataset-level processing, each step runs on every sample in memory before the next one starts, and a sample that fails a step is left out of the later ones. The annotation uses the *Annotation* panel's settings: pick `none` as the malignant-cell method for a quick pass, since copyKAT can take an hour per sample. If you move a slider afterwards, the results are out of date: run the steps again before saving (saving doesn't run them).
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
7. **Patient ID.** Under *Patient & sample info*, type the on-screen sample's patient ID. In dataset mode you can instead import a `.csv` with `sample_id` and `patient_id` columns; any other columns (tissue, stage, treatment, …) are imported too. These go into every cell's metadata (`patient`, …), into the per-cell table, and into `samplesheet.csv` and `summary.csv`. A sample sheet can carry a `patient` (or `patient_id`) column too.
8. **Annotation.** The *Annotation* panel holds the annotation script's settings, in sections: *Labelling methods* (SingleR and/or marker genes, and which one wins), *SingleR*, *Marker genes*, *Clustering*, *Malignant cells* and *Other settings*. Each section shows only its main settings; the rarely changed ones (custom reference, custom marker list and cutoffs, PCs / variable genes, the normal-reference regex, inferCNV's settings, cores, seed, …) are under its own **Advanced** link, which shows "(n changed)" when any of them differ from the script's defaults. Sections and settings for methods that aren't selected are hidden: choose SingleR only and the *Marker genes* section disappears; inferCNV's settings only appear when inferCNV is chosen; `label_priority` only when both labelling methods are. Labels are kept when only the malignant-cell settings change, so switching the method doesn't rerun SingleR. The *Annotation* tab annotates the cells that currently pass and shows the labels, clusters and malignant calls. Annotating runs in a **background process** (with the `callr` package), so the app stays responsive during a long copyKAT run, and the banner at the top shows it with *Stop after this one* / *Terminate now*. Saving includes the annotation if it was made for the cells being saved; if the cutoffs (or the doublet removal) changed since, the tab says so, and saving leaves it out until you annotate again. A sample loaded from a sample sheet uses its row's `annot_*` columns on top of the panel's settings (the tab lists them). See [`../R_scripts/annotation/README.md`](../R_scripts/annotation/README.md) for what each setting does.
9. **Save.** Choose an output folder, then click *Save this sample*, or *Save all samples in memory* in dataset mode. Saving writes what is there and runs nothing: the QC result, plus scDblFinder + DecontX and the annotation when they were run on the cells being saved. The save dialog lists what was left out and why (not run, or run before the cutoffs changed). You're asked before anything is overwritten, and a summary of each sample's changes from the defaults is shown.

### Dataset folders

```
datasets/GSE.../data/          <- point the app here
├── patient1/                  each sub-folder = one sample; its name becomes the sample ID
│   └── filtered_feature_bc_matrix/   (or outs/filtered_feature_bc_matrix.h5, GEO-prefixed files, ...)
├── patient2/
├── GSM0000003_counts.h5       scRNA only: a count file per sample also works
└── patient4.h5ad              any data type: one .rds or .h5ad file per sample also works
```

- **Finding each sample's input.** A sub-folder saved by the scripts or the app (`<sample>/<sample>.rds`) uses that object, so a dataset folder of an output folder can be loaded again. For scRNA, the app picks the 10x folder or `.h5` in each sub-folder, preferring *filtered* over *raw*, then a single count file (including `.rds` and `.h5ad`). For Visium and Xenium, the sub-folder itself is the input, since those scripts search inside it; a sub-folder holding just one `.rds` or `.h5ad` uses that file. `.rds` and `.h5ad` files sitting directly in the dataset folder are one sample each (a `_qc` / `_annotated` ending is dropped from the sample ID). Sub-folders with nothing readable are skipped and listed in a notice.
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
- **Resuming.** With *Skip samples already saved in the output folder* ticked, samples that already have a `<dataset>/<sample>/<sample>.rds` there are marked done and not loaded again, so a large dataset can be finished over several sessions.

### A dataset from a sample sheet

Choose *Dataset* and give a sample sheet `.csv` as the path (the same sheet the QC scripts read: `sample_id`, `path`, and any of the script's columns). Every row is loaded, in batches as for a folder, and **each row sets its own rules**:

- **QC settings** (`nmads_low`, `mito_ceiling`, `run_doublets`, …) filled in on a row stay with that sample: the Advanced settings don't change them for it. Blank cells take the Advanced settings. The heading above the sliders lists what the row sets.
- **`cutoff_<name>`** columns (e.g. `cutoff_max_mito = 12`, `Inf` = no maximum) put that sample's slider there, as if typed by hand. *Reset to MAD cutoffs* still goes back to the MAD rule.
- **`annot_<name>`** columns (e.g. `annot_malignant_method = none`, `annot_references = HPCA;Monaco`) override the Annotation panel's settings for that sample.
- `dataset` and `patient` (or `patient_id`) columns are used per row; the *Dataset* box fills in rows without one. Any other column (notes, …) is kept and written back to `samplesheet.csv`.

Relative paths are read from the sheet's folder, as in the scripts. An output folder's own `samplesheet.csv` works as a sample sheet: loading it reproduces that run, with every setting, the cutoffs set by hand and the annotation settings, and the output folder defaults to that folder.

## What gets saved

The scripts' layout, written by the scripts' own functions (`write_sample_outputs()` in the QC script, `write_annotation_outputs()` in the annotation script):

```
<output folder>/
├── samplesheet.csv        one row per sample: every setting used. Reusable as a sample sheet (app or script)
├── summary.csv            one row per sample: cells in, removed for each reason, kept, cutoffs used, labels
├── log.txt                every save adds to it: cutoffs (MAD value and value used), changes from the defaults
└── <dataset>/<sample>/
    ├── <sample>.rds           the QC'd object; with the annotation's labels when the annotation was saved too
    ├── <sample>_cells.csv     every cell, kept or not, with qc_reason (+ the label columns). Visium: _spots.csv
    ├── qc/                    the script's make_plots() figures
    ├── doublets_ambient/      scDblFinder and DecontX plots (scRNA-seq, if they ran)
    └── annotation/            annotation plots (if the annotation was saved)
```

Saving a sample replaces its folder and its rows. An annotation saved there earlier goes too (it was made on the old cells), unless it is saved again with it.

If the script can't draw a sample's plots (for example a Visium `.h5ad`, which has no tissue image), the app saves its own versions in `qc/` instead (`01_distributions_app.png`, `02_counts_vs_genes_app.png`, `03_spatial_qc_app.png`) and says so in `log.txt`.

### Documenting manual changes

The scRNA README asks that any change to QC parameters be reported:

- **`samplesheet.csv`** holds every parameter used for every sample. Cutoffs you set by hand are its `cutoff_*` columns, with the exact values, so giving the sheet to the QC script (or loading it in the app) reproduces the app's cutoffs exactly. The annotation settings are its `annot_*` columns.
- **`summary.csv`** has the cutoffs used (`cutoff_min_counts`, …) and `fixed_cutoffs`, naming the ones set by hand.
- **`log.txt`** gets, for each save, a line like `cutoffs min_counts=1585 (set by hand; MAD rule 548), …` and the sample's changes from the script's defaults in plain words, e.g. `nmads_mito: default 3 -> 2`. The save dialog shows the same.

## Xenium: new code, defaults need review

There was no Xenium QC script, so `R_scripts/xenium/xenium_qc.R` was written to match the scRNA and Visium scripts. Its default filters are first choices that still need review. See [`../R_scripts/xenium/README.md`](../R_scripts/xenium/README.md) for what it reads, the metrics, the defaults and how to run it on its own.

## Files

| File | What it is |
| --- | --- |
| `app.R` | the interface: sliders, settings, plots, save |
| `R/engine.R` | loads the scripts and maps each data type onto their functions |
| `R/export.R` | saves a sample with the scripts' writers, records the deviations from the defaults |
| `R/annotation.R` | runs the annotation script from the app, caches references, draws the Annotation tab |
| `R/background.R` | runs dataset-level steps in background R processes (callr) |
| `R/dataset.R` | dataset-folder scanning, reading a dataset's sample sheet, and the memory rule for batches |
| `R/h5ad.R` | reads `.h5ad` (AnnData) files with `hdf5r`, for all three data types; splits combined Visium `.h5ad` files into section folders |
| `R/geo_visium.R` | converts GEO Visium loose files (counts .csv + spot coordinates + image) into Space Ranger folders |
| `R/sliders.R` | slider ranges and the rule that turns slider positions into cutoffs |
| `R/plots_live.R` | the fast interactive plots (saved figures come from the scripts' `make_plots()`) |
| `../install_packages.R`, `../run_app.R` | setup and launcher (in the repo root) |

## Limits

- Advanced settings are shared by every sample in a dataset loaded from a folder. To give samples different *rules*, load the dataset from a sample sheet with those settings on their rows. *Cutoffs* can always be set per sample with the sliders.
- For very large samples, the interactive scatter and spatial plots draw a random subset (60k and 150k points); the counts and cutoffs always use every cell. The saved plots use every cell, except Xenium spatial maps, which are capped at 300k.
- **If the page disconnects** (laptop asleep, a browser tab put to sleep, a closed tab): what you had loaded is gone, since it lives in the page's session, but an annotation that was running **keeps going to the end**. Its result is cached. The next page that opens says which annotations are still running or have finished. Load the sample again with the same cutoffs and settings, and annotating it reuses that run's labels and malignant-cell results instead of computing them again (as long as R itself keeps running). Other background steps (QC, scDblFinder + DecontX) stop when the page disconnects.
- To avoid disconnects during long runs, keep the computer awake, and stop the browser from putting the tab to sleep: in Edge, *Settings → System and performance → Never put these sites to sleep*, add `127.0.0.1`; in Chrome, *Settings → Performance → Always keep these sites active*.
- Loading, and *Apply QC* on one sample, block the app until they finish. Large raw scRNA matrices (emptyDrops, at Apply QC) and DecontX can take minutes; with *Dataset-level processing* on, they run in the background instead.
- `.h5ad` files are read by the app, not the QC scripts (which don't read `.h5ad`). Results are the same, but to rerun a sample with a script, save its counts as `.rds` or 10x first (and change `path` and `format` in `samplesheet.csv`). Visium `.h5ad` files load without the tissue image (the spot coordinates come from `obsm['spatial']`).
