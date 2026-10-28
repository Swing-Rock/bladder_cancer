# scrna_qc — sample-sheet-driven QC for human scRNA-seq

`scrna_qc.R` runs quality control on any number of single-cell RNA-seq samples. You list your samples in a CSV (one row per sample) and run one command. For each sample, it removes broken, dying and doublet cells, estimates ambient-RNA contamination, and saves a clean Seurat object, plus tables and plots that show exactly what was removed and why.

- **One row per sample.** 
- **Data-driven cutoffs.** Each sample is judged against its own distribution, not a fixed rule.
- **Accomodates many input formats.** Cell Ranger folders, `.h5`, dense CSV/TSV tables, and Seurat/SCE `.rds` files.
- **Fails safely.** The whole sample sheet is checked before anything runs, and one broken sample doesn't stop the others.

---

## Contents

1. [Quick start](#quick-start)
2. [Installation](#installation)
3. [The sample sheet](#the-sample-sheet)
4. [Input formats](#input-formats)
5. [How the QC works](#how-the-qc-works)
6. [Outputs](#outputs)
7. [Reading the QC plots](#reading-the-qc-plots)
8. [Troubleshooting](#troubleshooting)
9. [Code layout](#code-layout)
10. [Where the defaults come from](#where-the-defaults-come-from)

---

## Quick start

**1. Write a sample sheet.** Only `sample_id` and `path` are required, though the `format` column is highly recommended:

```csv
sample_id,dataset,path, format
sample1,GSE130001,GSM3729178, 10x
sample2,GSE130001,GSM3729179, 10x
```

Rest of the columns are quality control parameters. Leave them as default values unless there's a specific reason to change them. Then document and report any changes made to default parameters 

Relative paths are relative to the folder the sample sheet is in. Absolute paths (`C:/...`, `/home/...`) work too.


**2. Run it** from the repo root (the folder with `bladder cancer.Rproj`):

```bash
Rscript 01_preprocessing+metadata/scRNA/scrna_qc.R datasets/GSE130001/samplesheet_GSE130001.csv
# outputs go to datasets/GSE130001/qc_output/; add a second argument to choose another output folder
```

In RStudio,  set `SAMPLE_SHEET_PATH` at the top of the script, and click **Source**. After sourcing once, you can run another sheet with `main("datasets/<other>/samplesheet.csv")`.

**3. Check the results.** After each run:

1. Open `qc_output/all_samples_summary.csv`: how many cells each sample kept, and why the others were removed.
2. Look at the plots in each dataset folder (see [Reading the QC plots](#reading-the-qc-plots)).
3. Document and report anything unusual about the dataset in the root README's dataset table.

---

## Installation

You need R ≥ 4.3 (it was developed with R 4.3 and Seurat 5).

```r
install.packages(c("Seurat", "Matrix", "ggplot2", "patchwork", "data.table"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("SingleCellExperiment", "scDblFinder", "celda"))

# Only if you need them:
BiocManager::install("DropletUtils")   # for unfiltered (raw) matrices when is_filtered = FALSE
install.packages("hdf5r")              # 10x .h5 files
```

The script checks up front that every package your sample sheet needs is installed. If one is missing, it tells you exactly what to install.

| Package | Used for |
| --- | --- |
| Seurat, Matrix | reading data, sparse matrices, the output objects |
| ggplot2, patchwork | QC plots |
| data.table | fast reading of large dense tables (optional) |
| SingleCellExperiment, scDblFinder | doublet detection |
| celda | DecontX ambient-RNA correction |
| DropletUtils | `emptyDrops` cell calling on unfiltered data |
| hdf5r | reading 10x `.h5` files |

---

## The sample sheet

Each row is one sample. Every column except `sample_id` and `path` is optional, and a **blank cell means "use the default"**. Columns the script doesn't recognize (e.g. `notes`) are ignored. `samplesheet_template.csv` has example rows for each kind of input.

### Columns

| Column | Default | What it does |
| --- | --- | --- |
| **What the sample is** | | |
| `sample_id` | *required* | Unique name. Becomes the output file name (`<sample_id>_qc.rds`) and the cell-name prefix. |
| `path` | *required* | 10x folder, `.h5`, dense table or `.rds` (see [Input formats](#input-formats)). |
| `dataset` | `dataset1` | dataset name |
| `format` | `auto` | `10x`, `10x_h5`, `dense`, `dense_t`, `rds`. `auto` guesses from the path. |
| `gene_column` | `2` | Which column of a 10x features file names the genes: 1 = Ensembl ID, 2 = gene symbol. |
| `mito_pattern` | `^MT-` | Regex for mitochondrial genes. Change only for unusual naming (e.g. `^GRCh38-MT-`). |
| `is_filtered` | `auto` | Whether the matrix contain empty droplets. `auto` = detected from the matrix (raw if more than half the barcodes have ≤ `emptydrops_lower` UMIs). `TRUE` = cells only. `FALSE` = raw droplet matrix; cells are called with emptyDrops first. An explicit `TRUE`/`FALSE` that contradicts the matrix stops that sample with an error. || `raw_path` | blank | Raw (unfiltered) matrix for a filtered sample. DecontX uses its empty droplets as the ambient profile. |
| **Cell filters** | | |
| `nmads_low` | `3` | MADs below the median (log10 scale) for UMIs and genes. |
| `nmads_high` | `5` | MADs above the median (log10 scale) for UMIs and genes. |
| `nmads_mito` | `3` | MADs above the median for % mito. |
| `min_counts` | `500` | UMI floor, always applied. |
| `min_features` | `200` | Gene floor, always applied. |
| `max_counts` | blank | Optional hard UMI cap. |
| `max_features` | blank | Optional hard gene cap. |
| `mito_floor` | `10` | The MAD mito cutoff never goes below this. |
| `mito_ceiling` | `25` | The MAD mito cutoff never goes above this. |
| **Doublets and ambient RNA** | | |
| `run_doublets` | `TRUE` | Run scDblFinder. |
| `doublet_rate` | blank | Expected doublet fraction (e.g. `0.08`). Blank = scDblFinder's default of ~1% per 1,000 cells. |
| `run_decontx` | `TRUE` | Run DecontX. |
| `max_contamination` | blank | Remove cells above this ambient fraction (e.g. `0.5`). Blank = keep all. |
| **Cell calling** (only when `is_filtered = FALSE`) | | |
| `emptydrops_lower` | `100` | Barcodes with ≤ this many UMIs define the ambient profile. |
| `emptydrops_fdr` | `0.001` | FDR cutoff for calling a barcode a cell. |
| **Genes and reproducibility** | | |
| `min_cells_gene` | `3` | Keep genes detected in at least this many of the sample's kept cells. |
| `seed` | `1234` | Random seed for scDblFinder, DecontX and emptyDrops. |

TRUE/FALSE columns accept `TRUE/FALSE`, `T/F`, `yes/no` or `1/0`.

### Checks before anything runs

The script validates the whole sheet first and lists **every** problem at once. Nothing is run until the sheet is clean, so a typo in row 30 doesn't surface an hour into the job. It checks for:

- missing `sample_id` or `path` columns, and duplicate sample IDs
- paths that don't exist (including `raw_path`)
- unknown formats, non-numbers in number columns, non-TRUE/FALSE in logical columns
- impossible settings (e.g. `mito_floor` > `mito_ceiling`, `nmads_*` ≤ 0)

---

## Input formats

With `format` left blank, the script guesses from the path: a folder → `10x`, `.h5` → `10x_h5`, `.rds` → `rds`, `.csv`/`.tsv`/`.txt` (gzipped or not) → `dense`. Set `format` yourself only for `dense_t`, or when the guess is wrong.

| Format | `path` points to | Notes |
| --- | --- | --- |
| `10x` | a folder with `matrix.mtx`, `barcodes.tsv`, and `genes.tsv` or `features.tsv` | Gzipped or not, with or without a GEO prefix (`GSM123_...`). Works for Cell Ranger v2 and v3. Antibody/CRISPR features are dropped. **One sample per folder.** |
| `10x_h5` | a Cell Ranger `.h5` file | Needs `hdf5r`. |
| `dense` | text table, **genes in rows**, cells in columns; first column = gene names | The common GEO `counts.csv.gz` layout. |
| `dense_t` | text table, **cells in rows**, genes in columns | If `dense` errors with "non-numeric column", try this. |
| `rds` | a Seurat object (v4 or v5), a SingleCellExperiment, or a count matrix | Uses raw counts from the RNA assay. Seurat v5 split layers are joined first. |

**What's inside a 10x folder.** The three files together describe one genes × cells table of counts:

- `barcodes.tsv.gz`: one cell barcode per line (e.g. `AAACCTGAGATTACCC-1`). These are the column names.
- `genes.tsv.gz` / `features.tsv.gz`: one gene per line. The Ensembl ID and symbol are tab-separated (v3 adds a third column for feature type). These are the row names.
- `matrix.mtx.gz`: only the non-zero counts, one per line as `gene_row cell_column count`. For example, `32709 1 84` means gene #32709 (MT-CYB) has 84 UMIs in cell #1.

**Always give the script raw counts.** If the values aren't whole numbers, the log warns that the data looks normalized. QC on normalized data is meaningless.

---

## How the QC works

Each sample goes through the same steps, independently:

```
read counts ─► call cells* ─► QC metrics ─► cutoffs ─► filter cells ─► doublets ─► ambient RNA ─► save .rds
                (*only if unfiltered)
```

### 1. Read the counts

Any input format is turned into a genes × cells sparse count matrix. Cell names are prefixed with the sample ID (`sample1_AAACCTGAGATTACCC-1`), so cells from different samples never share a name if you combine them later. Gene names containing `_` are changed to `-` (Seurat forbids underscores), and duplicate symbols are made unique (`GENE`, `GENE.1`).

### 2. Call cells (unfiltered data only)

A raw 10x matrix contains every barcode, and most of them are empty droplets holding only floating "ambient" RNA. When `is_filtered = FALSE`, **emptyDrops** separates cells from empty droplets:

1. Barcodes with ≤ `emptydrops_lower` (100) UMIs are assumed to be empty. Pooled together, they define what ambient RNA looks like.
2. Every other barcode is tested: is its gene-expression profile significantly different from that ambient profile?
3. Barcodes with FDR ≤ `emptydrops_fdr` (0.001) are called cells. That means at most ~0.1% of called cells are expected to be empty droplets.

The empty droplets are kept aside as the **background** for DecontX (step 6). If your data is already filtered but you have the raw matrix, give it as `raw_path` to get the same benefit.

### 3. Measure every cell

| Metric | What it is | What a bad value means |
| --- | --- | --- |
| `nCount_RNA` | total UMIs (RNA molecules) | Very low: broken cell or poor capture. |
| `nFeature_RNA` | number of genes detected | Very low: same as above. |
| `percent.mt` | % of UMIs from `MT-` genes | High: dying cell. Its membrane leaks cytoplasmic RNA, but mitochondrial RNA stays trapped inside the mitochondria. |
| `percent.ribo` | % from ribosomal protein genes | Reported only. Varies by cell type. |
| `percent.hb` | % from hemoglobin genes | Reported only. Helps you spot red blood cells at annotation. |
| `log10GenesPerUMI` | complexity | Reported only. Low = reads piled onto a few genes. |

### 4. Decide what counts as "unusual": MAD-based outlier detection

The script doesn't use one fixed rule for every dataset (like "every cell needs 500 genes"). It asks, **separately for each sample**: *is this cell much further from the typical cell than cells usually are?* This is the adaptive-threshold approach described in the quality-control chapter of OSCA, the Bioconductor single-cell book ([OSCA: Quality Control](https://bioconductor.org/books/release/OSCA.basic/quality-control.html); Amezquita et al. 2020).

For each metric:

1. **Typical value:** the **median**, the middle cell when all cells are sorted.
2. **Normal spread:** the **MAD** (median absolute deviation). Measure every cell's distance from the median and take the median of those distances. It's multiplied by 1.4826 so that it matches a standard deviation for bell-shaped data.

   $$\text{MAD} = 1.4826 \times \operatorname{median}\big(\lvert x_i - \operatorname{median}(x)\rvert\big)$$

3. **Cutoff:** flag cells more than *n* MADs from the median, in the "bad" direction:

| Metric | Scale | Flag if | Default *n* |
| --- | --- | --- | --- |
| UMIs, genes (low side) | log10 | below median − *n* × MAD | `nmads_low = 3` |
| UMIs, genes (high side) | log10 | above median + *n* × MAD | `nmads_high = 5` |
| % mito | as is | above median + *n* × MAD | `nmads_mito = 3` |

**Why the median and MAD instead of the mean and SD?** Extreme cells barely move them. Take % mito for 9 cells: `4, 5, 5, 6, 6, 6, 7, 8, 40`.

- Median = 6. Distances from 6, sorted: 0, 0, 0, 1, 1, 1, 2, 2, 34, so MAD = 1 × 1.4826 = 1.48. Cutoff = 6 + 3 × 1.48 = **10.4%**, and the dying 40% cell is flagged.
- With mean and SD instead: mean = 9.7, SD = 11.4, so the cutoff is about **44%** and the dying cell is *kept*. It inflated the SD enough to hide itself.

**Why log10 for UMIs and genes?** Counts vary by multiples: 2,000 and 16,000 are each 8× away from 4,000. On a log scale those distances are equal, and the distribution becomes roughly symmetric, which the MAD assumes. OSCA log-transforms these two metrics for the same reasons: to "improve resolution at small values" and to avoid negative thresholds ([OSCA: Quality Control](https://bioconductor.org/books/release/OSCA.basic/quality-control.html)).

**Why 3 MADs?** This is OSCA's default. Its rationale is that the cutoff "will retain 99% of non-outlier values that follow a normal distribution", so only clear outliers are removed ([OSCA: Quality Control](https://bioconductor.org/books/release/OSCA.basic/quality-control.html)). OSCA calls this a rough rationale, not a tested optimum. pipeComp (below) is the benchmark evidence. The high side uses a looser 5 MADs, because unusually large cells are better handled by doublet detection (step 5). A tight upper cutoff would also remove real large cells, such as tumor epithelium.

**Why per sample?** Samples differ in sequencing depth and quality. A cell that's typical in a shallow sample could look weak next to a deep one. Each sample is judged only against itself.

**Safety rails** on top of the MAD cutoffs:

- **Hard floors:** at least `min_counts` (500) UMIs and `min_features` (200) genes, always. This protects a uniformly poor sample, whose own median could be very low.
- **Mito clamp:** the MAD mito cutoff is clamped between `mito_floor` (10%) and `mito_ceiling` (25%). A very clean sample can't produce an over-strict 5% cutoff, and a messy one can't produce a meaningless 40% cutoff.
- **Optional caps:** `max_counts` and `max_features`. parameters are avaliable but not set due to sample to sample variations 

Each removed cell gets one `qc_reason`: the first check it failed, in the order low UMIs → high UMIs → low genes → high genes → mito.

### 5. Remove doublets (scDblFinder)

A doublet is two cells captured in one droplet. It looks like one cell expressing both cells' genes, and it can show up as a fake "hybrid" cell type (e.g. a fibroblast + T cell expressing COL1A1 and CD3E). **scDblFinder** finds them in three steps:

1. It clusters the cells and creates artificial doublets by adding together the counts of random pairs of real cells, mostly from different clusters.
2. It trains a classifier (gradient-boosted trees) to tell real cells from these artificial doublets.
3. It scores every real cell from 0 to 1. The cutoff is set so that the number called matches the expected doublet rate: about 1% per 1,000 cells for 10x Chromium (~3% for 3,000 cells, ~10% for 10,000). Use `doublet_rate` to override this.

It runs per sample, because doublets can only form within one capture. It uses raw counts. Samples with fewer than 50 cells skip this step with a warning.

Each sample gets its own **doublet score cutoff**, because the cutoff depends on how many doublets that sample is expected to have. The script records it as `doublet_threshold` in `qc_thresholds.csv`, along with `doublet_pct` (the % of scored cells called doublets), and draws it as a dashed line on `04_doublet_scores.png`. It also prints it in the log, e.g. `scDblFinder: 106 doublets among 3306 cells (3.2%; score cutoff 0.967)`. The cutoff is read from scDblFinder's own output. If that's unavailable, the script uses the lowest score among the cells scDblFinder called doublets.

### 6. Estimate ambient RNA (DecontX)

During dissociation, some cells burst. Their RNA floats in the suspension and gets captured in every droplet, so a fibroblast can carry a little epithelial RNA that isn't its own. **DecontX** models each cell's counts as a mixture of:

- its own cell population's expression profile, and
- a contamination profile pooled from the other populations (or, when a raw matrix is available, taken from the real empty droplets).

DecontX reports `decontX_contamination`, the estimated fraction of each cell's counts that is ambient, and a corrected count matrix. **No cells are removed here** unless you set `max_contamination`. The raw counts stay in the `RNA` assay and the corrected counts go into a separate `decontX` assay, so you can compare the two. The corrected counts are not whole numbers, so `round()` them for methods that need integers.

DecontX's default mode works well for single-cell data but poorly for **single-nucleus** data (Janssen et al. 2023). Be careful using it on snRNA-seq.

### 7. Drop rare genes and save

Genes detected in fewer than `min_cells_gene` (3) of the sample's kept cells are dropped. Each sample is then saved as its own Seurat object.

---

## Outputs

```
qc_output/
├── all_samples_summary.csv     every sample: cells in, removed by reason, kept, status
├── samplesheet_resolved.csv    the final value of every parameter, per sample
├── qc_log.txt                  timestamped log, including warnings
└── <dataset>/
    ├── <sample_id>_qc.rds      one Seurat object per sample
    ├── qc_thresholds.csv       the cutoffs actually used, per sample (incl. doublet score cutoff)
    ├── cell_qc_metadata.csv    every cell (kept and removed): metrics + qc_reason
    └── plots/                  01–05 PNGs, one panel per sample
```

**Each `<sample_id>_qc.rds`** contains:

- `RNA` assay: raw counts of the kept cells
- `decontX` assay: ambient-corrected counts (if DecontX ran)
- metadata: `sample`, `dataset`, `percent.mt`, `percent.ribo`, `percent.hb`, `log10GenesPerUMI`, `scDblFinder.score`, `decontX_contamination`

---

## Reading the QC plots

A good cutoff trims a tail and leaves the main body of the distribution alone. If a dashed red line cuts through the middle of a violin, the cutoff is too strict for that sample.

| Plot | Looks fine when | Change something when |
| --- | --- | --- |
| `01_violin_before_after.png`: UMIs, genes, % mito, before (grey) vs after (blue), with cutoffs dashed | Only the tails are removed. | A cutoff cuts into the main body → loosen that `nmads_*`. |
| `02_counts_vs_genes_by_reason.png`: UMIs vs genes, coloured by why each cell was removed | Removed cells sit at the bottom-left or off the diagonal. | A big coloured group sits inside the grey cloud → that rule is too strict. |
| `03_counts_vs_mito.png`: UMIs vs % mito | High-mito cells are mostly low-UMI, i.e. damaged. | Many *high-UMI* cells are above the line → it may be a real high-mito cell type; raise `mito_ceiling`. |
| `04_doublet_scores.png`: doublet scores, with each sample's cutoff dashed and labelled | Most cells score near 0, with a small group near 1, and the cutoff sits in the gap between them. | Doublets are far more than ~5–10% of a sample, or the cutoff cuts through a big peak → check `doublet_rate`. |
| `05_decontX_contamination.png` | Most cells are below ~0.2. | A large group is above 0.5 → consider `max_contamination` or the `decontX` assay. |

**if you decide to make any changes to the quality control parameters, please make sure to report and document your changes**

---

## Troubleshooting

| Message | Likely cause and fix |
| --- | --- |
| `Sample sheet problems (nothing was run)` | Fix each listed row; the message names the row and column. |
| `Missing packages: ...` | Install them with the command printed, or set `run_doublets` / `run_decontx` to `FALSE`. |
| `WARNING: no genes match the mito pattern '^MT-'` | Gene names are Ensembl IDs, mitochondrial genes were removed, or they use a different prefix. Check head(rownames(...)); for 10x input make sure gene_column is 2, otherwise convert IDs to symbols or set mito_pattern |
| `WARNING: counts are not whole numbers` | The file holds normalized data. Find the raw counts. |
| `More than one matrix.mtx file in ...` | Several samples share one folder. Put each sample in its own folder. |
| `non-numeric column(s) after the first` | The dense table is the other way round. Set `format = dense_t`. |
| `[sample] FAILED: ...` | That sample couldn't be read. The others still finish, and the failure is listed in `all_samples_summary.csv`. |
| `WARNING: scDblFinder failed` / `DecontX failed` | That step was skipped for that sample only. See `qc_log.txt` for the error. |

When run with `Rscript`, the script exits with a non-zero status if any sample failed, so a SLURM job shows as failed rather than silently half-done.

---

## Code layout

The file is organized top to bottom as:

| Section | What's there |
| --- | --- |
| Settings | `SAMPLE_SHEET_PATH`, `OUT_DIR` (used when run without arguments) |
| Defaults | `DEFAULTS`: every sample-sheet column and its default value |
| Logging | `start_log()`, `log_msg()` |
| Sample sheet | `read_sample_sheet()`, `resolve_sample()`, `check_packages()` |
| Readers | `read_counts()` and one reader per format |
| QC steps | `load_sample()` → `call_cells()` → `compute_qc_metrics()` → `compute_thresholds()` → `apply_cell_filters()` → `detect_doublets()` → `correct_ambient()`, chained by `qc_one_sample()` |
| Outputs | `save_sample()`, `summarise_sample()`, `make_plots()`, `write_dataset_outputs()`, `write_overview()` |
| `main()` | runs everything in order |
| Last lines | call `main()` |

To change a default for everyone, edit `DEFAULTS`. To change it for one sample, use the sample sheet. To add a new numeric parameter, add it to `DEFAULTS` and use `p$<name>` in the step that needs it. Type checking picks it up automatically.

---

## Where the defaults come from

Most defaults follow published recommendations or benchmarks. A few are conventions

| Choice | Basis |
| --- | --- |
| Per-sample thresholds | Heumos et al. 2023: QC "is performed at the sample level as thresholds can vary substantially between samples". |
| MAD outliers, log scale, 3 MADs | OSCA default ([Quality Control chapter](https://bioconductor.org/books/release/OSCA.basic/quality-control.html); Amezquita et al. 2020). |
| 5 MADs on the high side | Lenient filtering favored by pipeComp (Germain et al. 2020), which found stringent global filters discard real subpopulations. |
| 3 MADs for % mito | OSCA and sc-best-practices-style pipelines. |
| `mito_floor = 10%` | Osorio & Cai: recommended 10% standard for human tissues (5% for mouse). |
| `mito_ceiling = 25%` | Yates et al. 2025 found malignant cells often have higher mito %, and that common 10–20% cutoffs can remove viable tumor cells |
| `min_counts = 500`, `min_features = 200`, `min_cells_gene = 3` | Convention (Seurat tutorial values). They're only floors. |
| scDblFinder | Top performer in its authors' 16-dataset benchmark (Germain et al. 2022), which also reports it ranking best in an updated version of Xi & Li's independent benchmark. In Xi & Li's published 2021 benchmark (which predates it), DoubletFinder was the most accurate. Recommended by Heumos et al. 2023 for accuracy and speed. |
| Doublet rate ~1% per 1,000 cells | 10x multiplet-rate table (~0.8% at 1,000 cells, ~7.6% at 10,000). |
| DecontX | Performed similarly to CellBender on scRNA-seq in Janssen et al. 2023 (not recommended for snRNA-seq defaults). |
| emptyDrops FDR 0.001 | OSCA default (Lun et al. 2019). |

**Known departures from published guidance.** You may want to revisit these:

- Doublets are detected *after* the cell filters. scDblFinder's authors advise against stringent filtering before it, though the MAD filters here are moderate.
- A cell is removed if it fails *any one* metric. pipeComp's milder default requires failing at least two.

### References

- Amezquita RA, Lun ATL, Hicks SC, Gottardo R et al. (2020). Orchestrating single-cell analysis with Bioconductor. *Nature Methods* 17:137–145. <https://doi.org/10.1038/s41592-019-0654-x>
- OSCA online book, *Basics of Single-Cell Analysis with Bioconductor*, Chapter 1: Quality Control. <https://bioconductor.org/books/release/OSCA.basic/quality-control.html>
- Germain P-L, Sonrel A, Robinson MD (2020). pipeComp, a general framework for the evaluation of computational pipelines, reveals performant single cell RNA-seq preprocessing tools. *Genome Biology* 21:227. <https://pmc.ncbi.nlm.nih.gov/articles/PMC7465801>
- Heumos L et al. (2023). Best practices for single-cell analysis across modalities. *Nature Reviews Genetics*. <https://pmc.ncbi.nlm.nih.gov/articles/PMC10066026>
- Osorio D, Cai JJ. Systematic determination of the mitochondrial proportion in human and mice tissues for single-cell RNA-sequencing data quality control. <https://www.biorxiv.org/content/10.1101/2020.02.20.958793>
- Yates J, Kraft A, Boeva V (2025). Filtering cells with high mitochondrial content depletes viable metabolically altered malignant cell populations in cancer single-cell studies. *Genome Biology* 26:91. <https://link.springer.com/article/10.1186/s13059-025-03559-w>
- Germain P-L et al. (2022). Doublet identification in single-cell sequencing data using scDblFinder. *F1000Research* 10:979. <https://f1000research.com/articles/10-979>
- Yang S et al. (2020). Decontamination of ambient RNA in single-cell RNA-seq with DecontX. *Genome Biology* 21:57. <https://doi.org/10.1186/s13059-020-1950-6>
- Xi NM, Li JJ (2021). Benchmarking computational doublet-detection methods for single-cell RNA sequencing data. *Cell Systems* 12(2):176–194. <https://doi.org/10.1016/j.cels.2020.11.008>
- Hao Y et al. (2024). Dictionary learning for integrative, multimodal and scalable single-cell analysis (Seurat v5). *Nature Biotechnology*. <https://pmc.ncbi.nlm.nih.gov/articles/PMC10928517>
- Janssen P et al. (2023). The effect of background noise and its removal on the analysis of single-cell expression data. *Genome Biology*. <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10278251/>
- Lun ATL et al. (2019). EmptyDrops: distinguishing cells from empty droplets in droplet-based single-cell RNA sequencing data. *Genome Biology* 20:63. <https://doi.org/10.1186/s13059-019-1662-y>
