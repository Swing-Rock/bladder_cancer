# Visium preprocessing: QC and deconvolution

Two scripts for 10x Visium spatial transcriptomics, run in this order:

| Script | What it does | Input | Output |
| --- | --- | --- | --- |
| `h5ad_to_spaceranger.py` (only if needed) | Splits a combined Scanpy `.h5ad` into one Space-Ranger-style folder per section | A Visium `.h5ad` with images in `uns["spatial"]` | Per-section folders + a sample sheet for `visium_qc.R` |
| `visium_qc.R` | Removes low-quality spots with per-sample MAD cutoffs (same method as [`../scRNA/scrna_qc.R`](../scRNA/README.md)) and maps where the removed spots are | A sample sheet listing Space Ranger folders (any layout) or Seurat `.rds` files | One `<sample>_qc.rds` per sample (raw counts + tissue image), tables, plots |
| `visium_deconvolution.R` | Estimates the cell-type proportions in every spot with RCTD | The `*_qc.rds` files, plus an **annotated single-cell Seurat object** as the reference | `<sample>_deconv.rds`, per-spot proportion tables, spatial maps |

Human data only (mitochondrial genes are matched with `^MT-`).

---

## Quick start

Run from the repo root (open `bladder cancer.Rproj` in RStudio, or `cd` there first).

**1. Write a sample sheet** in your dataset folder, e.g. `datasets/<ID>/samplesheet_visium.csv`. Only `sample_id` and `path` are required, and paths are relative to the sheet's folder:

```csv
sample_id,dataset,path
slide1_A1,MY_VISIUM,data/slide1_A1/outs
slide1_B1,MY_VISIUM,data/slide1_B1/outs
```

**2. Run QC.**

```bash
Rscript 01_preprocessing+metadata/visium/visium_qc.R datasets/<ID>/samplesheet_visium.csv
# -> datasets/<ID>/qc_output/<dataset>/<sample>_qc.rds, tables, plots/
```

**3. Run deconvolution**, once you have an annotated scRNA-seq reference:

```bash
Rscript 01_preprocessing+metadata/visium/visium_deconvolution.R \
  datasets/<ID>/qc_output/<dataset> \
  datasets/<reference>/reference_annotated.rds \
  celltype
# -> datasets/<ID>/qc_output/<dataset>/deconvolution/
```

The three arguments are: the folder of `*_qc.rds` files, the reference `.rds`, and the name of the reference's cell-type metadata column. In RStudio, set the same three values at the top of the script and click **Source**.

**4. Check the results.** Outputs go to `datasets/<ID>/qc_output/`, which is git-ignored along with all data. After each run:

1. Open `all_samples_summary.csv`: spots kept per section, and why the others were removed.
2. Look at `plots/03_spatial_removed_spots.png` for every section, to confirm that no real tissue region was cut out (see [How the QC works](#how-the-qc-works)).
3. After deconvolution, read `deconvolution_log.txt` (reference cell types, genes shared, skipped spots) and check that the proportion maps make anatomical sense.
4. Note anything unusual about the dataset in the root README's dataset table.

### Example: GSE171351

GSE171351 ships all four sections in one `.h5ad`, so it needs the conversion step first (see [Combined `.h5ad` files](#combined-h5ad-files-scanpyanndata)):

```bash
# 1. split the combined .h5ad into one folder per section (also writes the sample sheet)
python  01_preprocessing+metadata/visium/h5ad_to_spaceranger.py \
        datasets/GSE171351/data/GSE171351_combined_visium.h5ad datasets/GSE171351/data/sections \
        --name-col Patient --dataset GSE171351

# 2. QC
Rscript 01_preprocessing+metadata/visium/visium_qc.R datasets/GSE171351/samplesheet_visium.csv

# 3. deconvolution, once the annotated bladder reference exists
Rscript 01_preprocessing+metadata/visium/visium_deconvolution.R \
        datasets/GSE171351/qc_output/GSE171351 <reference.rds> celltype
```

---

## Installation

R ≥ 4.3.

```r
install.packages(c("Seurat", "Matrix", "ggplot2", "patchwork", "hdf5r"))

# RCTD (deconvolution only). Either version works; the script detects which one is installed.
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("spacexr", "SpatialExperiment"))   # Bioconductor version, R >= 4.5
# or the original GitHub version:
# remotes::install_github("dmcable/spacexr", build_vignettes = FALSE)
```

spacexr exists in two versions with different function names. The Bioconductor version (since Bioconductor 3.21) uses `createRctd()` / `runRctd()` on SpatialExperiment objects. The original GitHub version uses `Reference()` / `SpatialRNA()` / `create.RCTD()` / `run.RCTD()`. The script checks which is installed, uses the matching calls, and logs which one it used.

---

## Input layouts `visium_qc.R` can read

Set `format` to `spaceranger` for a folder or `rds` for a Seurat object, or leave it blank to guess from the path.

| Layout | What the script looks for |
| --- | --- |
| **Space Ranger `outs/`** | `filtered_feature_bc_matrix.h5` (or the `filtered_feature_bc_matrix/` folder) plus `spatial/tissue_positions.csv`, `scalefactors_json.json`, `tissue_lowres_image.png`. |
| **GEO download** | The same files flattened into one folder, with a prefix (`GSM123_tissue_positions_list.csv.gz`), images and JSON gzipped, Space Ranger v1 positions without a header, and the matrix as `.mtx` + `barcodes.tsv`. If there's no `features.tsv` in the sample folder, the script looks **one folder up**, since some GEO series share one features file. |
| **Seurat `.rds`** | A Seurat object with a `Spatial` assay (and ideally an image). |

### Combined `.h5ad` files (Scanpy/AnnData)

Some GEO series ship all sections in one AnnData file instead of Space Ranger folders, e.g. GSE171351's `GSE171351_combined_visium.h5ad`. Split it into one Space-Ranger-style folder per section with `h5ad_to_spaceranger.py`, then run `visium_qc.R` on the sample sheet it writes:

```bash
pip install anndata scipy pandas pillow
python 01_preprocessing+metadata/visium/h5ad_to_spaceranger.py \
  datasets/GSE171351/data/GSE171351_combined_visium.h5ad \
  datasets/GSE171351/data/sections \
  --name-col Patient --dataset GSE171351
# -> datasets/GSE171351/data/sections/<Patient>_<library>/ + datasets/GSE171351/samplesheet_visium.csv
Rscript 01_preprocessing+metadata/visium/visium_qc.R datasets/GSE171351/samplesheet_visium.csv
```

The converter:

- **Counts:** takes raw counts from `layers["counts"]`, `.raw` or `.X`, in that order, and stops if none of them are whole numbers (i.e. only normalized data was shared).
- **Sections:** finds which `obs` column assigns spots to sections (the one whose values match the library IDs in `uns["spatial"]`), or uses `--sample-col`.
- **Barcodes:** strips the suffixes Scanpy adds when concatenating (`AAAC...-1-0-0-0` → `AAAC...-1`).
- **Spatial files:** writes positions from `obsm["spatial"]` (Scanpy convention: column 0 = image x, column 1 = image y), the scale factors, and the low- and high-res images.

Check the result once by eye: open `plots/03_spatial_removed_spots.png` and a `SpatialFeaturePlot` of one section, and confirm the spots sit on the tissue.

Files are found by name pattern anywhere inside the sample folder, so extra files don't matter. If a folder has both `filtered_` and `raw_` matrices, the filtered one is used. Each sample needs its own folder.

**Only in-tissue spots are kept.** If the matrix includes every spot on the capture area, the script keeps only those marked `in_tissue = 1` in the positions file and logs how many that was.

---

## How the QC works

1. **Read** the counts, positions, scale factors and low-res image (hi-res if there's no low-res), and attach the image to the Seurat object. Spot names get the sample ID as a prefix (`slide1_A1_AAACAAGTATCTCCCA-1`). Each spot's array and pixel positions are saved in the metadata (`array_row`, `array_col`, `pxl_row`, `pxl_col`), so plots and deconvolution don't depend on Seurat's image format.
2. **Measure each spot:** UMIs (`nCount_Spatial`), genes (`nFeature_Spatial`) and % mitochondrial. % ribosomal and % hemoglobin are recorded but never used to remove spots.
3. **Remove outlier spots** with per-sample MAD cutoffs. This is the same method, defaults and safety rails as the scRNA-seq QC; the [scRNA README](../scRNA/README.md#4-decide-what-counts-as-unusual-mad-based-outlier-detection) explains the algorithm in detail:
   - UMIs and genes: more than `nmads_low` (3) MADs below or `nmads_high` (5) MADs above the median, on a log10 scale
   - % mito: more than `nmads_mito` (3) MADs above the median, clamped to `mito_floor`–`mito_ceiling` (10–25%)
   - hard floors: at least `min_counts` (250) UMIs and `min_features` (200) genes per spot
4. **Drop rare genes** (seen in fewer than `min_spots_gene` = 3 kept spots) and save `<sample>_qc.rds` with **raw counts**.

**Different from scRNA-seq QC, on purpose:**

- **No doublet detection or ambient-RNA correction.** A 55 µm Visium spot holds several cells by design, so "doublet" doesn't apply, and scRNA ambient-RNA tools assume one cell per droplet.
- **No normalization.** RCTD needs raw counts. Normalize (e.g. `SCTransform(obj, assay = "Spatial")`) in your downstream analysis.
- **The floors are lower** (250 UMIs, 200 genes) because a spot's counts depend on the tissue under it. These match the lower cutoffs of the earlier fixed-threshold Visium QC.

**Check where the removed spots are.** In tissue, low counts aren't always bad quality: stroma, necrotic or fibrotic areas, and the tissue edge capture less RNA. A real region could be flagged just for being a minority of the slide. `plots/03_spatial_removed_spots.png` shows this:

- **Scattered** removed spots, or ones at a torn edge, are what you expect.
- A **solid block** of removed spots matching a histological region means the cutoff is cutting out real tissue. Loosen `nmads_low` for that sample and rerun.

---

## How the deconvolution works

Each spot is a mixture of several cells. **RCTD** (Robust Cell Type Decomposition; Cable et al. 2022) does the following:

1. It learns an average expression profile for each cell type from the single-cell reference.
2. It models each spot's counts as a weighted mix of those profiles, accounting for the technical differences between the scRNA-seq and Visium platforms.
3. It estimates the weights: the proportion of the spot made up by each cell type.

The script runs RCTD in **`full` mode**, which allows any number of cell types per spot. That's appropriate for Visium's multi-cell spots. The other modes (`doublet`, `multi`) assume one or two cells per spot and are meant for higher-resolution platforms like Slide-seq.

### The reference

The reference can be any annotated single-cell Seurat object. It sets what RCTD can find, so it matters more than any parameter:

- **Same tissue.** Use bladder / urothelial carcinoma scRNA-seq, ideally your own QC'd and annotated data. A reference from another cancer will mislabel bladder-specific cell states.
- **Every cell type that's in the tissue.** RCTD can only assign spots to types it knows. A missing type gets split among the most similar ones that are present.
- **Raw counts.** The script stops if the counts aren't whole numbers.
- **At least 25 cells per type.** Smaller types are dropped, with a log message, because RCTD can't learn a profile from so few cells. Types with more than 10,000 cells are randomly downsampled to 10,000.
- **Labels at the resolution you care about.** For CAF states, label the fibroblast subtypes (e.g. myCAF, iCAF) in the reference, not just "Fibroblast". Very similar subtypes are harder to separate in 55 µm spots, so check that the proportion maps make anatomical sense.
- **Matching gene symbols.** The log reports how many genes the reference and the spots share, and warns if it's under 5,000.

`/` in cell-type names is replaced with `-`, because RCTD rejects it. Names are also cleaned up to be valid column names, e.g. `CD8 T` becomes `RCTD_CD8.T`.

---

## Outputs

```
datasets/<ID>/qc_output/
├── all_samples_summary.csv          every sample: spots in/kept, removed by reason, status
├── samplesheet_resolved.csv         final value of every parameter per sample
├── qc_log.txt
└── <dataset>/
    ├── <sample>_qc.rds              kept spots, raw counts, tissue image, positions in metadata
    ├── qc_thresholds.csv            cutoffs actually used, per sample
    ├── spot_qc_metadata.csv         every in-tissue spot: metrics, position, qc_reason
    ├── plots/
    │   ├── 01_violin_before_after.png
    │   ├── 02_counts_vs_genes_by_reason.png
    │   ├── 03_spatial_removed_spots.png    where removed spots are on the tissue
    │   ├── 04_spatial_umis.png             UMIs across the tissue
    │   └── 05_spatial_mito.png             % mito across the tissue
    └── deconvolution/                       written by visium_deconvolution.R
        ├── <sample>_deconv.rds              QC'd object + RCTD_<type>, RCTD_dominant, RCTD_max_proportion
        ├── <sample>_proportions.csv         one row per spot
        ├── <sample>_dominant_cell_type.png
        ├── <sample>_proportions.png         one map per cell type
        ├── reference_cell_types.csv         cell types and cell counts actually used
        ├── deconvolution_summary.csv        mean proportion of each type per sample
        └── deconvolution_log.txt
```

`RCTD_max_proportion` is how much the dominant type accounts for. Spots where it's low are true mixtures, so `RCTD_dominant` means little for them. RCTD also applies its own minimum-UMI cutoff per spot (100 UMIs by default in the GitHub version); any spot it skips gets `NA` proportions and is logged.

---

## Sample sheet columns (`visium_qc.R`)

| Column | Default | What it does |
| --- | --- | --- |
| `sample_id` | *required* | Unique name; output file name and spot-name prefix. |
| `path` | *required* | Sample folder or `.rds`, relative to the sheet's folder. |
| `dataset` | `dataset1` | Output folder name. Samples are never merged. |
| `format` | `auto` | `spaceranger` or `rds`. |
| `gene_column` | `2` | Features file column with gene names: 1 = Ensembl ID, 2 = symbol. |
| `mito_pattern` | `^MT-` | Regex for mitochondrial genes. |
| `nmads_low` | `3` | MADs below the median (log10) for UMIs and genes. |
| `nmads_high` | `5` | MADs above the median (log10) for UMIs and genes. |
| `nmads_mito` | `3` | MADs above the median for % mito. |
| `min_counts` | `250` | UMI floor per spot, always applied. |
| `min_features` | `200` | Gene floor per spot, always applied. |
| `max_counts`, `max_features` | blank | Optional hard caps. |
| `mito_floor`, `mito_ceiling` | `10`, `25` | Clamp for the MAD mito cutoff (%). |
| `min_spots_gene` | `3` | Keep genes seen in at least this many kept spots. |
| `seed` | `1234` | Random seed. |

`samplesheet_template_visium.csv` has example rows for each input layout.

## Deconvolution settings (`PARAMS` at the top of `visium_deconvolution.R`)

| Setting | Default | What it does |
| --- | --- | --- |
| `reference_assay` | `RNA` | Assay in the reference with raw counts. |
| `min_cells_per_type` | `25` | Drop reference cell types with fewer cells. |
| `max_cells_per_type` | `10000` | Downsample bigger types to this many cells. |
| `doublet_mode` | `full` | RCTD mode; keep `full` for Visium. |
| `n_cores` | `4` | Parallel cores. |
| `min_shared_genes` | `5000` | Warn below this many shared genes. |
| `seed` | `1234` | Random seed (downsampling). |

---

## Troubleshooting

| Message | Cause and fix |
| --- | --- |
| `No ... file found under ...` | A required Space Ranger file is missing. The QC needs a count matrix, `tissue_positions`, `scalefactors_json.json` and a tissue `.png`. |
| `More than one ... file under ...` | Two samples share a folder. Give each sample its own folder. |
| `No barcodes shared between the count matrix and the in-tissue spots` | The matrix and positions file come from different samples, or barcodes were renamed. |
| `WARNING: no genes match the mito pattern` | Genes are Ensembl IDs → set `gene_column = 2`. |
| `Column '...' not found in the reference metadata` | Wrong cell-type column name. The error lists the available columns. |
| `counts are not whole numbers` | The reference has normalized data. Use the raw-count assay (`reference_assay`). |
| `WARNING: only N shared genes` | Reference and spots use different gene identifiers or annotations. |
| `Package 'spacexr' (RCTD) is not installed` | See [Installation](#installation). |

---

## Changes from the earlier Visium scripts (`st_data_qc.R`, `st_deconvolution.R`)

| Before | Now | Why |
| --- | --- | --- |
| Fixed cutoffs (200–7,500 genes, 250–50,000 UMIs, mito < 15%, ribo < 40%) | Per-sample MAD cutoffs, with 250/200 floors and a 10–25% mito clamp | Same method as the scRNA-seq QC, and it adapts to each section's depth. The ribosomal filter was dropped because ribosomal content varies by cell type. |
| One sample per function call | Sample sheet, many samples, one failure doesn't stop the rest | Same workflow as the scRNA-seq QC. |
| `.mtx` loader renamed spots `spot_1, spot_2, ...` and attached positions by row order | Spots keep their barcodes and are matched to positions by barcode; out-of-tissue spots removed | Row-order matching misaligns spots whenever the positions file lists all spots or a different order. |
| `.mtx` loader didn't attach the image | Image attached for every layout | Spatial plots work for GEO downloads too. |
| Wrote a new `.h5` into the data folder as a side effect | Data folders are only read | Inputs stay untouched. |
| SCTransform inside QC | QC saves raw counts | Normalization belongs downstream, and RCTD needs raw counts. |
| Depended on `visium_scripts/st_helper.R` | Self-contained | Plot helpers rewritten inside the script. `.h5ad` input dropped. |
| Reference hard-coded to Puram et al. 2017 HNSCC (GSE103322, log-TPM converted to pseudo-counts) | Any annotated Seurat object, given at run time | Bladder tissue needs a bladder reference, and RCTD expects UMI counts. |
| GitHub spacexr only | Bioconductor or GitHub spacexr | The Bioconductor release renamed the functions. |

---

## References

- Cable DM, Murray E, Zou LS, Goeva A, Macosko EZ, Chen F, Irizarry RA (2022). Robust decomposition of cell type mixtures in spatial transcriptomics. *Nature Biotechnology* 40:517–526. <https://pmc.ncbi.nlm.nih.gov/articles/PMC8606190/>
- spacexr on Bioconductor: <https://bioconductor.org/packages/spacexr>
- MAD-based QC: see the references in [`../scRNA/README.md`](../scRNA/README.md#references) (OSCA; Germain et al. 2020; Heumos et al. 2023).
