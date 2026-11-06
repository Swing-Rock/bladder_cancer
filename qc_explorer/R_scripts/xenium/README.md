# xenium_qc — sample-sheet-driven QC for 10x Xenium (human data)

`xenium_qc.R` filters low-quality cells from Xenium runs, one Seurat object per sample. It follows the same conventions as [`../scRNA/scrna_qc.R`](../scRNA/) and [`../visium/visium_qc.R`](../visium/): the same sample-sheet rules, the same per-sample MAD cutoffs with hard floors and caps, and the same output layout. The QC app (`app_scripts/`) uses its functions directly.

> **New script; defaults need review.** There was no Xenium QC in the lab pipeline before. The default filters below are first choices, not settled lab standards. Check the plots on your own data before relying on them.

## Contents

- [Quick start](#quick-start)
- [Installation](#installation)
- [Input](#input)
- [Sample sheet columns](#sample-sheet-columns)
- [How the QC works](#how-the-qc-works)
- [Outputs](#outputs)
- [Reading the QC plots](#reading-the-qc-plots)
- [Troubleshooting](#troubleshooting)

---

## Quick start

```bash
# from the repo root
Rscript R_scripts/xenium/xenium_qc.R datasets/<ID>/samplesheet_xenium.csv
# -> datasets/<ID>/qc_output/<dataset>/<sample>_qc.rds, tables, plots/
```

In RStudio, set `SAMPLE_SHEET_PATH` at the top of the script and click **Source**. Start the sample sheet from [`samplesheet_template_xenium.csv`](samplesheet_template_xenium.csv). Or use the QC app's **Xenium** mode to set cutoffs by eye.

---

## Installation

```r
install.packages(c("Seurat", "Matrix", "ggplot2", "patchwork", "data.table", "hdf5r"))
install.packages("arrow")   # only if a sample has cells.parquet but no cells.csv.gz
```

The script checks for the packages it needs before it starts.

---

## Input

One row per Xenium region or output folder. `path` can be:

- **A Xenium output folder** (e.g. `outs/`). The script searches inside it, so GEO downloads with prefixed file names work too. It needs:
  - the count matrix: `cell_feature_matrix.h5`, or the `cell_feature_matrix/` folder (`matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz`)
  - the per-cell table: `cells.csv.gz`, or `cells.parquet` (needs the arrow package)
- **A Seurat `.rds` made with `LoadXenium()`.** These objects usually lack cell areas, so the area filters are skipped.

Only `Gene Expression` features become the counts. Negative-control probes, negative-control codewords and unassigned/deprecated codewords are kept only as a per-cell control count.

---

## Sample sheet columns

Required: `sample_id`, `path`. A blank cell means "use the default". Relative paths are relative to the sheet's own folder.

| Column | Default | Meaning |
| --- | --- | --- |
| `dataset` | `dataset1` | Output folder name; samples are never merged. |
| `format` | `auto` | `xenium` (output folder) or `rds`. `auto` decides from the path. |
| `nmads_low` | `3` | MADs below the median (log scale) for transcripts and genes per cell. |
| `nmads_high` | `5` | MADs above the median (log scale) for transcripts and genes per cell. |
| `nmads_area` | `5` | MADs below and above the median (log scale) for cell area. |
| `min_counts` | `10` | Hard floor on transcripts per cell, always applied. |
| `min_features` | `5` | Hard floor on genes per cell, always applied. |
| `max_counts`, `max_features` | blank | Optional hard caps. |
| `min_cell_area`, `max_cell_area` | blank | Optional hard floor / cap on cell area (µm²). |
| `max_pct_control` | `5` | Maximum % of a cell's decoded counts that come from negative controls. |
| `min_cells_gene` | `0` | Keep genes detected in at least this many kept cells. `0` keeps the whole panel. |
| `seed` | `1234` | For reproducible random steps. |

Extra columns (e.g. `notes`) are ignored, and the log says so.

---

## How the QC works

1. **Read** the matrix and the cells table, and match them by cell ID. Cell names are prefixed with the sample ID.
2. **Measure every cell:**

   | Metric | What it is |
   | --- | --- |
   | `nCount_Xenium` | transcripts (gene-expression counts) |
   | `nFeature_Xenium` | genes detected |
   | `cell_area`, `nucleus_area` | segmentation areas (µm²) from the cells table; `nucleus_ratio` = nucleus / cell |
   | `pct_control` | negative-control counts as a % of all decoded counts in the cell: a per-cell false-positive rate |

3. **Set this sample's cutoffs**, the same way as the scRNA and Visium scripts. Transcripts and genes use median − `nmads_low` MADs and median + `nmads_high` MADs on the log scale, but never looser than the hard floors or caps. Cell area uses median ± `nmads_area` MADs on the log scale, within the optional hard limits. `pct_control` uses the fixed cap.
4. **Flag cells** outside any cutoff. The first failing check becomes `qc_reason`: `low_counts`, `high_counts`, `low_features`, `high_features`, `small_area`, `large_area` or `high_control`.
5. **Save** the kept cells with raw counts.

**Not done here, on purpose:** doublet detection and ambient-RNA correction. Those tools assume droplet-based data. Imaging-based segmentation errors (merged or split cells) show up instead as unusual cell areas and transcript counts, which the area and count cutoffs catch.

---

## Outputs

```
qc_output/
├── all_samples_summary.csv     every sample: cells in, removed by reason, kept, status
├── samplesheet_resolved.csv    the final value of every parameter, per sample
├── qc_log.txt
└── <dataset>/
    ├── <sample_id>_qc.rds      one Seurat object per sample
    ├── qc_thresholds.csv       the cutoffs actually used, per sample
    ├── cell_qc_metadata.csv    every cell (kept and removed): metrics + qc_reason
    └── plots/                  01–04 PNGs, one panel per sample
```

Each `<sample_id>_qc.rds` has:
- an `Xenium` assay with the raw counts of the kept cells
- an FOV named `fov` with the cell centroids, so `ImageDimPlot()` works
- metadata: `sample`, `dataset`, `control_counts`, `pct_control`, `cell_area`, `nucleus_area`, `nucleus_ratio`, `x`, `y`

For `.rds` input, the original object is subset instead, which keeps its own FOVs.

---

## Reading the QC plots

| Plot | Looks fine when | Change something when |
| --- | --- | --- |
| `01_violin_before_after.png`: transcripts, genes, cell area, % control, before (grey) vs after (blue), cutoffs dashed | Only the tails are removed. | A cutoff cuts into the main body → loosen that `nmads_*` or hard limit. |
| `02_counts_vs_genes_by_reason.png` | Removed cells sit at the low end or off the diagonal. | A large coloured group sits inside the grey cloud → that rule is too strict. |
| `03_spatial_removed_cells.png`: removed cells on the tissue | Removed cells are scattered, or sit at tissue edges and folds. | Removed cells form a whole region (e.g. one tissue compartment) → a real low-RNA cell type may be being lost; loosen `nmads_low` / `min_counts`. |
| `04_spatial_transcripts.png` | Transcript density follows tissue structure. | A patch of very low counts → a section or imaging problem worth noting. |

Spatial maps draw at most 300,000 cells per sample (a random subset). The counts and cutoffs always use every cell.

---

## Troubleshooting

| Message | Fix |
| --- | --- |
| `No cell_feature_matrix.h5 ...` / `No cells.csv.gz ...` | `path` must be the Xenium output folder, or a folder containing those files. |
| `only cells.parquet was found, which needs the arrow package` | `install.packages("arrow")`. |
| `No cell IDs shared between the matrix and the cells table` | The matrix and cells table come from different runs or regions. |
| `no negative-control counts found; % control will be 0` | Older Xenium output without control columns; the `pct_control` filter has no effect. |
| `More than one ... file under ...` | Put each sample in its own folder. |
