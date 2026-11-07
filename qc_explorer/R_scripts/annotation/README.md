# annotate_cells — cell-type annotation (SingleR and/or marker genes) and malignant-cell calling (copyKAT / inferCNV)

`annotate_cells.R` labels the cells of QC'd samples. It runs after `scRNA/scrna_qc.R`, `visium/visium_qc.R`, `xenium/xenium_qc.R`, or the QC app itself, and before merging and Harmony. The QC app runs the same functions from its **Annotation** panel.

## Contents

- [Quick start](#quick-start)
- [Installation](#installation)
- [Input](#input)
- [What it does, by data type](#what-it-does-by-data-type)
- [How it works](#how-it-works)
- [Outputs](#outputs)
- [Interrupted runs and reruns](#interrupted-runs-and-reruns)
- [Running copyKAT on Expanse](#running-copykat-on-expanse)
- [Settings](#settings)
- [Adding a malignant-cell method](#adding-a-malignant-cell-method)
- [Caveats](#caveats)

---

## Quick start

```bash
# from the repo root (where "bladder cancer.Rproj" is)
Rscript R_scripts/annotation/annotate_cells.R datasets/GSE130001/qc_output
# -> the labels go into each qc_output/<dataset>/<sample>/<sample>.rds and <sample>_cells.csv,
#    plots into <sample>/annotation/, the settings into samplesheet.csv, the labels into summary.csv
```

In RStudio, set `INPUT` (and `PARAMS` if needed) at the top of the script and click **Source**.

The input is a QC output folder (from the QC scripts or the app), or its `samplesheet.csv`. Every sample in that sheet is annotated; without a `samplesheet.csv`, every `<dataset>/<sample>/<sample>.rds` in the folder. The sheet's `patient` column is added to every cell.

**Per-sample settings.** `PARAMS` applies to every sample. To change a setting for one sample, give its row in `samplesheet.csv` an `annot_<setting>` column, e.g. `annot_malignant_method = none` or `annot_references = HPCA;Monaco` (several values separated by `;`, `none` for an empty list). Blank cells use `PARAMS`. After a run, the sheet's `annot_*` columns hold every setting each sample was annotated with, so rerunning the script on the same folder reproduces it.

---

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("SingleR", "celldex", "scran", "BiocParallel", "SummarizedExperiment"))
install.packages(c("pheatmap", "viridis", "scales"))

# malignant-cell calling (scRNA-seq)
remotes::install_github("navinlabcode/copykat")   # malignant_method = "copykat" (the default)
BiocManager::install("infercnv")                  # malignant_method = "infercnv" or "copykat+infercnv"; also needs JAGS: https://mcmc-jags.sourceforge.io
```

The built-in references download once from ExperimentHub and are cached. On an HPC node without internet, run the script once on a login node first, or point `custom_reference` at a local object.

The script checks for every package the chosen settings need before it starts, and lists what to install.

---

## Input

The script takes Seurat objects with **raw counts**, one per sample: the `<dataset>/<sample>/<sample>.rds` files of a QC output folder (from the QC scripts or the app). It finds the right assay itself (`RNA`, `Spatial` or `Xenium`). Set `assay = "decontX"` to annotate scRNA-seq on the ambient-corrected counts instead.

---

## What it does, by data type

| | Cell labels (SingleR / markers) | Cluster labels | Malignant calls (CNV) |
| --- | --- | --- | --- |
| **scRNA-seq** | per cell | yes | yes |
| **Xenium** | per cell, but only on the panel's genes | yes | no (a targeted panel can't show copy-number changes) |
| **Visium** | each spot's *dominant* type only | yes | no |

For Visium, a 55 µm spot holds several cells, so one label per spot hides the mixture. Use `visium/visium_deconvolution.R` (RCTD) for cell-type proportions; SingleR on spots is only a rough check. For Xenium, a few hundred panel genes give coarse labels. A custom single-cell reference restricted to the same genes works best.

---

## How it works

### 1. SingleR (per cell)
Each cell's log-expression is correlated with every cell type in the reference, using the reference's marker genes, and the best match wins.

- `singler_<ref>_label` is the label SingleR assigned.
- `singler_<ref>_pruned` is the same label, set to `NA` where the call is low-confidence (the winner barely beats the runner-up).
- `singler_label` is the pruned label from the *primary* reference: the custom reference if you gave one, otherwise the first built-in one.

Each cell is labelled independently of every other cell. So annotating a sample on its own gives the same labels as annotating a merged object, and no integration is needed first.

### 1b. Marker genes (per cell), `label_methods` includes `markers`
A second, independent way to label cells: score each cell for each cell type's marker genes.

- **Marker list:** `MARKERS` in the script has canonical markers for the major cell types of bladder tumours: Epithelial (EPCAM, KRT8/18/19, uroplakins, …), T cells, NK cells, B cells, plasma cells, myeloid cells, mast cells, fibroblasts, endothelial cells, and smooth muscle / pericytes. Add your own types, or replace a built-in one, with `custom_markers`: a `.csv` with columns `cell_type` and `gene`, one row per gene. With `marker_list_mode = "replace"` only your list is used. A type needs at least 2 of its markers in the data to be scored; the log lists the types it skipped.
- **Score:** each marker gene's log-expression is divided by its 99th percentile in the sample (capped at 1), so every marker counts the same whatever its expression level. A type's score is the mean over its markers, from 0 (none expressed) to 1.
- `marker_label`: the type with the highest score. `marker_score`: that score. `marker_delta`: how far it beats the runner-up.
- `marker_pruned`: `marker_label`, set to `NA` when the score is below `marker_min_score` or the lead is below `marker_min_delta`.

With clustering on, `marker_cluster_label` labels each cluster by its average marker expression, and `marker_cluster_majority` / `_frac` give the most common per-cell marker label. Cluster labels are steadier than per-cell ones when the data are sparse.

Marker labels are easy to check and to explain, and they work without a reference. SingleR uses the whole transcriptome and can name finer types. Running both and looking where they disagree is a quick way to find clusters that need a manual look. `plots/<sample>/06_marker_dotplot.png` shows every marker by final label.

### 1c. The final label: `cell_label`
`label_methods` picks the method(s): `SingleR`, `markers`, or both. When both run, `label_priority` says which wins:

- `cell_label`: the priority method's confident label (`singler_label` or `marker_pruned`). Where it has none (pruned to `NA`), the other method's label fills in.
- `cell_label_source`: which method each cell's label came from.
- `cell_lineage`: `cell_label` as a broad lineage (`T cells`, `Myeloid`, `Epithelial`, …; see `LINEAGES` in the script). The two methods and the references name types differently (`T_cells`, `CD8+ T-cells`, `T cells`), so use this column to compare or count across methods and samples.
- `label_agreement` (both methods only): `agree` / `disagree` (same broad lineage or not), `one method only` (just one made a confident call), or `NA` (neither did).

`cell_label` is the label the rest of the pipeline uses: the normal reference for the malignant-cell step, the plots, and the app.

### 2. Cluster-level labels (`cluster_level = TRUE`)
Each sample is clustered on its own (Seurat: normalize, variable genes, PCA, neighbours, Louvain, UMAP). Two checks on the per-cell calls are added:
- `singler_<ref>_cluster_label`: SingleR run on each cluster's pooled expression.
- `singler_<ref>_cluster_majority` and `_frac`: the most common per-cell label in the cluster, and its share of the cluster's cells.

If the per-cell and cluster labels disagree, that cluster needs a manual look.

### 3. Malignant cells (scRNA-seq)
This step is its own function, `call_malignant()`, run after the SingleR labels. `malignant_method` picks the method: `none` (skip the step), `copykat`, `infercnv` or `copykat+infercnv`. The choices come from `MALIGNANT_METHODS` in the script, and the app's dropdown lists the same ones, so a new method shows up in both once it's added there (see [Adding a malignant-cell method](#adding-a-malignant-cell-method)).

SingleR can't tell tumour from normal epithelium: both are "Epithelial cells". Tumour cells, though, carry large copy-number changes (in MIBC: 9p/9q loss, 1q gain, 17p loss, …). copyKAT and inferCNV estimate these from expression. They average expression over windows of neighbouring genes along each chromosome and compare each cell with known-normal cells. A run of genes that is consistently higher or lower than normal means a gain or a loss.

- **Normal reference:** cells whose final label (`cell_label`) matches `cnv_reference_types` (immune and stromal types, by default), taken from the same sample. Every other cell is tested.
- **copyKAT** calls each cell `aneuploid` or `diploid`. With fewer than `cnv_min_reference` reference cells, copyKAT finds diploid cells itself.
- **inferCNV** gives a denoised CNV profile per cell. Each cell gets two numbers:
  - `cnv_score`: the mean squared deviation from normal.
  - `cnv_cor`: the correlation with the average profile of the sample's top-scoring cells.

  A cell is called malignant when its score is above the `cnv_score_quantile` of the reference cells' scores **and** `cnv_cor > cnv_cor_cutoff`. This is the approach of Tirosh et al. 2016 and Puram et al. 2017.
- **`cnv_status`** combines the calls of the chosen method(s): `malignant`, `normal`, `normal (reference)`, `uncertain` (the methods disagree), `unresolved` (no call) or `not tested`. Every method's columns are always written (blank for methods not run), so the output has the same columns whatever is chosen.

Both tools run per sample, because each patient's tumour has its own CNV profile and the best normal baseline is the same sample's own normal cells.

---

## Outputs

Everything goes into the QC output folder, next to the QC results:

```
<qc output folder>/
├── samplesheet.csv              annot_* columns: every setting each sample was annotated with
├── summary.csv                  annotated_at, cell_labels (the commonest final labels), malignant (cnv_status counts)
├── log.txt                      every run adds to it: settings that differ from PARAMS, references, labels
├── working_files/<dataset>/<sample>/   copyKAT / inferCNV working files and cached results (safe to delete)
└── <dataset>/<sample>/
    ├── <sample>.rds             the QC'd object, now with the labels (+ cnv_status, patient); raw counts untouched
    ├── <sample>_cells.csv       the label columns added: final label, SingleR, markers, scores, CNV calls
    └── annotation/
        ├── 01_score_heatmap_<ref>.png   each cell vs every label; a clean call = one bright row
        ├── 02_delta_<ref>.png           how far each winning label is ahead; low = pruned
        ├── 03_umap_labels.png           clusters, final label, lineage, agreement, each method's per-cell and cluster labels
        ├── 04_cnv_status.png            malignant calls on the UMAP + inferCNV score vs correlation
        ├── 05_*                         copyKAT's and inferCNV's own heatmaps
        └── 06_marker_dotplot.png        marker genes by final label (when markers run)
```

`<sample>.rds` also keeps the raw SingleR results, the settings and what the labels were computed from, in `obj@misc$annotation` (used to reuse the labels; see below). The scaled data made for clustering is not saved (it is large; `NormalizeData` / `ScaleData` again if you need it).

Running QC on a sample again (script or app) replaces its object, so its annotation goes with it: annotate it again afterwards. The app saves an annotation with the same function, `write_annotation_outputs()`, when you save a sample you annotated in it.

---

## Interrupted runs and reruns

The malignant-cell step (copyKAT especially) is the slow part, so the script protects the work before it:

1. **Labels are saved first.** Right after SingleR and clustering, `<sample>.rds` is written with the labels (and `obj@misc$annotation$stage = "labels saved"`). If the run is then interrupted (Ctrl+C, a closed window, a killed job), the labels are still there. The file is rewritten with the malignant calls when the step finishes (`stage = "complete"`).
2. **Rerunning picks up where it stopped** (`reuse_previous = TRUE`, the default). The script reuses the labels already in `<sample>.rds` when they were made from the same cells, counts and label settings. It reuses a method's cached result in `working_files/` when the cells, reference cells and that method's settings are the same. The log says what it reused. Change anything those depend on and that part is computed again.
3. **Switching methods is cheap.** Running with `malignant_method = "none"` first gives the labels in minutes; a later run with `copykat` reuses them and only runs copyKAT. Switching back to a method that already ran reuses its result.

Set `reuse_previous = FALSE` to compute everything again. In the app, the same reuse works within a session: previews keep their labels and method results in the app's temporary folder, so changing only the malignant-cell settings doesn't rerun SingleR.

---

## Running copyKAT on Expanse

copyKAT is slow and memory-hungry on a laptop. On Windows it also runs on a single core, so a sample of a few thousand cells can take an hour or more and tie up most of the RAM. If you don't want to run it locally, run it on Expanse with the provided sbatch script instead. There `n_cores` actually parallelizes it, and nothing on your machine can interrupt it.

> **To be added:** the sbatch script and instructions for using it.

A convenient split:

1. **Locally:** run with `malignant_method = "none"` (or pick *none* in the app) to get the SingleR labels and clusters in a few minutes.
2. **On Expanse:** run with `malignant_method = "copykat"` using the sbatch script.

Because labels are reused when nothing they depend on has changed (see [above](#interrupted-runs-and-reruns)), copying the QC output folder (with the labelled `<sample>.rds` files) lets the Expanse run skip SingleR and go straight to copyKAT. For one sample only, set `annot_malignant_method = copykat` on its row in `samplesheet.csv`.

---

## Settings

All settings are in `PARAMS` at the top of the script; the app shows the same ones. Any of them can be set for one sample with an `annot_<setting>` column in `samplesheet.csv` (see [Quick start](#quick-start)).

| Setting | Default | Meaning |
| --- | --- | --- |
| `label_methods` | `SingleR` | `SingleR`, `markers`, or both. In the app, a multi-select. |
| `label_priority` | `SingleR` | When both run: the method that fills `cell_label` (the other fills its gaps). |
| `references` | `HPCA`, `Blueprint` | Built-in references for SingleR (see `REFERENCES`). `Monaco` / `DICE` add fine immune subsets. |
| `label_level` | `main` | `main` = broad types; `fine` = subtypes. |
| `custom_reference`, `custom_label_col` | blank | An annotated single-cell Seurat object and its label column, e.g. a published bladder atlas. **Needed for states the atlases don't have** (myCAF / iCAF, urothelial subtypes). |
| `primary` | blank | Reference that fills `singler_label`. Blank = custom if given, else the first built-in. |
| `custom_markers` | blank | Your own marker list: a `.csv` with columns `cell_type`, `gene`. Relative paths are found the same way as `infercnv_gene_order`. |
| `marker_list_mode` | `add` | `add`: your types are added to the built-in list (a type with the same name replaces the built-in one). `replace`: only your list. |
| `marker_min_score`, `marker_min_delta` | `0.1`, `0.05` | A cell's top marker score must reach `marker_min_score` (0–1) and beat the runner-up by `marker_min_delta`, or `marker_pruned` is `NA`. Raise them for stricter calls. |
| `assay` | `RNA` | `decontX` = ambient-corrected counts from `scrna_qc.R`. |
| `cluster_level`, `cluster_resolution` | `TRUE`, `0.8` | Per-sample clustering and cluster labels. |
| `malignant_method` | `copykat` | `none`, `copykat`, `infercnv` or `copykat+infercnv` (formerly `cnv_method`; the old value `both` still works). |
| `cnv_reference_types` | immune + stromal regex | Regex on the final label (`cell_label`) picking the normal reference cells. It covers the built-in marker types (T_cells, Myeloid, Fibroblasts, …) and the common reference names; check it matches your custom reference or marker type names. |
| `cnv_score_quantile`, `cnv_cor_cutoff` | `0.95`, `0.4` | inferCNV malignant cutoffs (see above). |
| `infercnv_gene_order` | blank | Gene positions for inferCNV (`gene  chr  start  end`). Blank = copyKAT's built-in hg38 table. Or use `hg38_gencode_v27.txt` from https://data.broadinstitute.org/Trinity/CTAT/cnv/. A relative path (e.g. `references/hg38_gencode_v27.txt`) is looked up in the folder R runs in, then in `qc_explorer/`, then in the folder holding it (`bladder cancer`), so the same value works in the script and the app; `custom_reference` works the same way. |
| `cnv_max_cells` | `15000` | Larger samples run CNV on a random subset (copyKAT is slow); untested cells get `not tested`. |
| `reuse_previous` | `TRUE` | Reuse saved labels and cached malignant-cell results when the cells and settings are unchanged (see above). `FALSE` = compute everything again. |

---

## Adding a malignant-cell method

Methods are listed in `MALIGNANT_METHODS` near the top of `annotate_cells.R`. To add one (for example an expression-based classifier such as scCancer2):

1. Write `run_<method>(counts, ref_cells, sample_id, work_dir, params)`. It gets the raw counts of the tested cells, the normal reference cells and a working folder. It returns `NULL` if it skipped the sample, or a list with:
   - `call`: a named vector over the tested cells: `"malignant"`, `"normal"` or `NA`;
   - `columns`: a data.frame (row names = cells) with the method's own output, e.g. a probability;
   - `heatmap`: the path to a figure to copy into `plots/`, or `NA`.
2. Add an entry to `MALIGNANT_METHODS`: a display `label`, the `packages` it needs (with install commands), the `params` its result depends on, the metadata `columns` it adds (with their blank values), and `run`.
3. Optionally add combinations to `MALIGNANT_COMBOS`, e.g. `"copykat+mymethod"`.

The method then appears in `malignant_method` and in the app's dropdown, is included in `cnv_status`, and gets cached and reused like the others. A method written in Python can be called from `run_<method>()` with `system2()`: write the counts to a file, run the script, and read its per-cell output back.

---

## Caveats

- **Near-diploid tumours** (some low-grade, FGFR3-driven NMIBC) can be missed by CNV-based calling. This is less of a problem in MIBC.
- **Few normal cells** (e.g. epithelium-enriched sorts) give a weak inferCNV baseline, and inferCNV is skipped below `cnv_min_reference`. copyKAT can still run without one.
- **Cell-cycle and stress programs** can look like small CNVs. Trust arm-scale changes in the heatmaps.
- **The reference sets what can be found.** Labels not in the reference come out as the nearest type that is, or as `NA` after pruning.
- **Runtime:** copyKAT can take an hour or more for ~10k cells. Run large datasets with the script on Expanse rather than in the app; see [Running copyKAT on Expanse](#running-copykat-on-expanse).

## References

- Aran et al. 2019, *Nat Immunol* — SingleR.
- Gao et al. 2021, *Nat Biotechnol* — copyKAT.
- Tirosh et al. 2016, *Science*; Puram et al. 2017, *Cell* — inferCNV-based malignant scoring.
- inferCNV of the Trinity CTAT Project: https://github.com/broadinstitute/infercnv
