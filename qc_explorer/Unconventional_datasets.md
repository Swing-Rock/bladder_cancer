# Unconventional datasets

Most GEO series load as they are: point the app (or a sample sheet) at a Cell Ranger / Space Ranger folder, an `.h5`, an `.rds` or an `.h5ad`. The series below are packaged differently. This page says how to load each one and what to watch for. If a new dataset doesn't fit any of these, add it here.

## Quick reference

| Dataset | Data type | What's unusual | What to do |
| --- | --- | --- | --- |
| [GSE149652](#gse149652) | scRNA-seq | Plain CSV count tables; sorted T cells only | Load the per-sample `*_count_matrice.csv.gz` files as `dense`; skip the malignant-cell step |
| [GSE171351](#gse171351) | Visium | Every section in one `.h5ad` | Pick the `.h5ad`; the app splits it into one folder per section |
| [GSE246011](#gse246011) | Visium | Loose files per sample (CSV counts + coordinates + TIFF), no Space Ranger folder; mixed cancer types | Put the bladder samples' files in their own folder and load it; the app converts them |

| [GSE345644](#gse345644) | Xenium | Cell table is `cells.parquet.gz` | Unzip it to `cells.parquet` before QC |
| [DUTRENEO](#dutreneo-zenodo-19693316) | Xenium | One 30 GB Zenodo zip of full Xenium folders; cell table is `cells.parquet`; response labels in a separate file | Extract only `cell_feature_matrix.h5` + `cells.parquet` per sample; put labels in the sheet's `notes` |
| [GSE319536](#gse319536) | Visium | Loose per-sample GEO files with GSM prefixes and `.gz`, no Space Ranger folder | Rename into `<sample>/filtered_feature_bc_matrix.h5` + `<sample>/spatial/`; `format = spaceranger` |
| [zenodo17999961 / zenodo20199229](#zenodo17999961--zenodo20199229-cervilla-et-al) | Visium CytAssist / Xenium | Multi-cancer platform-comparison records; the bladder Xenium is a duplicate of DUTRENEO DU39 | Download only the `*_Bladder.zip` files; patient ID `DU39` |
| [Bulk RNA-seq cohorts](#bulk-rna-seq-imvigor210-gse176307-uc-genome) | bulk RNA-seq | No bulk script in the repo; three different formats | Processed with stand-alone Python scripts (in each upload) |

---

## GSE149652

**scRNA-seq, sorted CD4⁺ and CD8⁺ T cells from bladder tumours (anti-PD-L1 cohort).**

GEO ships each sample as three CSV files:

```
GSM4506643_Anti-PD-L1_A_tumor_CD4_droplet_count_matrice.csv.gz      <- the counts: use this one
GSM4506643_Anti-PD-L1_A_tumor_CD4_droplet_cellinfo_matrice.csv.gz   per-cell metadata (not read by QC)
GSM4506643_Anti-PD-L1_A_tumor_CD4_TCR_paired_matrice.csv.gz         paired TCR data (not read by QC)
```

**How to load it**

- **Use the per-sample files** from `GSE149652_RAW.tar`, one sample per file. Don't use the series-level `GSE149652_CD4_TIL_*` / `GSE149652_CD8_TIL_*` files: they combine several samples in one matrix, and QC treats each file as one sample.
- **Input type `dense`.** The `*_count_matrice.csv.gz` files are plain count tables. `dense` expects genes in rows and cells in columns. If loading fails with *"non-numeric column(s) after the first"*, the table is the other way round: use `dense_t`.
- The cellinfo and TCR files aren't needed for QC. To use them later (e.g. clonotypes), join them to the QC'd cells by barcode.
```
```

**Annotation**

The dataset description says the cells may all be T cells. Check the SingleR labels first, if the samples are mostly T cells, **skip the malignant-cell step** (copyKAT / inferCNV): there are no tumour cells to find, so the calls would be meaningless. 
* if running the scripts: Set `annot_malignant_method = none` on these samples' rows in `samplesheet.csv`.
* if using the app: 

---

## GSE171351

**Visium, several bladder tumour sections stored together in one AnnData file** (`GSE171351_combined_visium.h5ad`).

**How to load it**

In the app, choose data type **Visium** and pick the `.h5ad` file itself (in either *One sample* or *Dataset* mode), or its folder as a dataset. The app:

1. sees that the file holds several sections (more than one library in `uns['spatial']`);
2. splits it into one Space-Ranger-style folder per section, **next to the `.h5ad`**, each with its counts, spot positions, scale factors and H&E images. Folders are named `<Patient>_<library>`;
3. writes `GSE171351_combined_visium_sections.csv` listing the sections. This is also a sample sheet for `visium_qc.R`;
4. loads the sections as dataset `GSE171351`, with output in `datasets/GSE171351/qc_output/`.

Loading the file again reuses the folders (unless the `.h5ad` changed). While `_sections.csv` exists, a folder scan skips the combined file, so the data isn't loaded twice. Delete the section folders and `_sections.csv` to split again.

If you are using the script instead of the app, `R_scripts/visium/h5ad_to_spaceranger.py` makes the same folders (see `R_scripts/visium/README.md`).

> **If anything goes wrong while splitting** (an error, a section missing, an image that looks wrong), please let me know.

---

## GSE246011

**Visium across several cancer types. Only 2 samples are bladder (`BLCA`); the rest are stomach adenocarcinoma (`STAD`).**

GEO ships each sample as loose files instead of a Space Ranger folder:

```
GSM7853987_BLCA-B1-count.csv.gz              counts: spots in rows, genes in columns (raw whole numbers)
GSM7853987_BLCA-B1_spot_coordinates.csv.gz   spot positions (Space Ranger tissue_positions layout)
GSM7853987_BLCA-B1.tif.gz                    H&E image (3000 x 3000 px)
GSM7853987_BLCA-B1-Metadata.csv.gz           spot metadata (not needed for QC)
```

**How to load it**

1. **Keep only the bladder samples.** `GSE246011_RAW.tar` contains every sample. Put the `BLCA` samples' files in their own folder (e.g. `datasets/GSE246011/data/`); otherwise loading the folder converts and loads the stomach samples too. Picking one of a sample's files loads just that sample.
2. In the app, choose **Visium** and pick that folder (or any one of a sample's files).
3. The app converts each sample into a Space-Ranger-style folder `<GSM>_<name>/` next to the files, lists them in `visium_converted.csv`, and loads them. Converting takes about a minute per sample; later loads reuse the folders.

**Check the image alignment**

GEO doesn't include Space Ranger's `scalefactors_json.json`, which says how spot coordinates map onto the image. The app works it out by fitting the spots to the tissue: in-tissue spots should land on stained tissue, the rest on bare glass. For `BLCA-B1` this works well (68% of in-tissue spots on tissue vs 3% of the others), but it's automatic, so:

- **open `spatial/alignment_check_app.png`** in each converted sample folder. Green circles = in-tissue spots, blue = the rest. Green should cover the tissue and blue should follow the gaps and edges;
- the log gives each sample's fit, and a weak fit is flagged in red. If a sample's spots are visibly shifted or scaled, its spatial plots are wrong; tell me and don't rely on them.

The alignment only affects where spots are drawn on the image. The counts and the QC metrics themselves don't depend on it.

Reading `.tif` images needs the `tiff` R package (`install_packages.R` installs it).

---

## GSE345644

**Xenium, one healthy bladder sample (480-gene panel).**

GEO ships the cell table as `GSM10012265_cells.parquet.gz`. `xenium_qc.R` looks for `cells.csv.gz` or `cells.parquet`, so it doesn't find a gzipped parquet.

**How to load it**

- Download `cell_feature_matrix.h5` and `cells.parquet.gz`, then `gunzip -k` the parquet. Reading `.parquet` needs the `arrow` R package.
- Also download `cells_barcodes_labels.csv.gz`: the authors' own cell labels (Basal, Intermediate, Umbrella, SMC, xFB-0…9, `low_quality_cell`). Useful for checking label transfer later.

**What to watch for**

- Very sparse (median 18 transcripts/cell); cutoffs end up at the hard floors (10 transcripts, 5 genes) and 66% of cells are kept.
- Every removed cell is one the authors labelled `low_quality_cell`, and all their real cell types are 100% kept. About 28k kept cells are still author-labelled low quality: consider dropping them downstream.

---

## DUTRENEO (Zenodo 19693316)

**Xenium (377-gene Multi-Tissue panel), 20 MIBC patients from the DUTRENEO trial, with treatment arm and pathological response.** GEO GSE328930 has the matched bulk RNA-seq (counts + metadata only, no Xenium).

Zenodo has one zip, `Dutreneo_xenium.zip` (~30 GB), holding a full Xenium output folder per patient (`DU15/` … `DU39/`, including large images QC doesn't need), plus `mdata_xenium.txt` with the labels.

**How to load it**

- Extract only what QC needs (a few GB instead of 30):
  ```bash
  unzip -n Dutreneo_xenium.zip '*/cell_feature_matrix.h5' '*/cells.parquet' '*/metrics_summary.csv' -d raw
  ```
  The cell table is `cells.parquet` (there is no `cells.csv.gz`), so `arrow` is needed.
- Sample sheet: `sample_id` = folder name (DU15…), `path` = `raw/DU15`. Put the arm, response, TIS group and biopsy ID from `mdata_xenium.txt` in `notes`. (`mdata_xenium.txt` has one more value per row than header names: the first value is an R row number.)
- `ID` in `mdata_xenium.txt` is the biopsy ID that links to the GEO bulk samples (`Sample_<ID>`); a few are written differently (`2709203-A`, `271510327-1881`).

**What to watch for**

- Depth varies a lot (median 6–81 transcripts/cell), and % kept follows it (36–95%). DU23, DU20, DU18 and DU17 keep < 50%.
- In those low-depth samples, immune cells are kept less often than epithelial cells (e.g. DU23: 47% vs 82%), because small, RNA-poor cells fall below the 10-transcript floor. Spatial immune features (density, CD8–tumour distance) should be depth-adjusted or checked without these samples.
- DU39 is the same Xenium run as zenodo20199229 `xenium_Bladder` (identical files; see below).

---

## GSE319536

**Visium (FFPE probe set, 18,085 genes), 22 primary MIBC.** Same study as GSE326226 (Xenium).

GEO ships each sample as loose, prefixed, gzipped files, not a Space Ranger folder:

```
GSM9519394_BC01_filtered_feature_bc_matrix.h5
GSM9519394_BC01_tissue_positions.csv.gz
GSM9519394_BC01_scalefactors_json.json.gz
GSM9519394_BC01_tissue_hires_image.png.gz
GSM9519394_BC01_tissue_lowres_image.png.gz
```

**How to load it**

- Rebuild the Space Ranger layout: `raw/BC01/filtered_feature_bc_matrix.h5` and `raw/BC01/spatial/{tissue_positions.csv, scalefactors_json.json, tissue_hires_image.png, tissue_lowres_image.png}` (strip the prefix, `gunzip` the rest).
- In the sample sheet, `format` must be `spaceranger` (`visium` is rejected).

**What to watch for**

- BC09: only 277 tissue spots and ~40% mito → not recommended.
- BC22: removed spots cover intact tissue including most tumour nests (median 36 UMIs/spot; overlay in `BC22/qc/qc_on_HE.png`) → not recommended.
- BC03 and BC20 keep < 50%, but the removed spots are necrotic/debris areas on the H&E, so the kept spots are usable.
- The mito cutoff hits its 25% cap in 8 samples (typical for FFPE).
- When matching QC'd spots back to the raw data, cell names in the `.rds` are prefixed with the sample (`BC22_AACACC…-1`).

---

## zenodo17999961 / zenodo20199229 (Cervilla et al.)

**Platform-comparison study across six cancer types; the same tumours on Visium v1/CytAssist/Visium HD (zenodo17999961), Xenium + CosMx (zenodo20199229) and Xenium MT/5K (zenodo18000256).**

**How to load it**

- Each record has one zip per cancer type and platform: download only the bladder ones (`cytassist_Bladder.zip`, 48 MB; `xenium_Bladder.zip`, 1.2 GB). `visiumhd_Bladder.zip` (3.4 GB) was skipped because Visium HD is on hold.
- `cytassist_Bladder.zip` is a normal Space Ranger `outs` folder (`format = spaceranger`). It is a Gene + Protein run, so the `.h5` has two modalities; QC uses the RNA only.
- `xenium_Bladder.zip`: extract `cell_feature_matrix.h5` + `cells.parquet` (as for DUTRENEO).

**What to watch for**

- **`xenium_Bladder` is DUTRENEO DU39.** Both files are byte-identical (same md5) and QC gives the same numbers (216,976 cells in, 188,933 kept). Don't count it as a separate sample.
- So the CytAssist (and Visium HD) bladder sample is the DU39 tumour on other platforms, a labelled sample (NAC arm, complete response). All of them use patient ID `DU39`.
- CytAssist QC: 2,340 of 2,428 spots kept (96.4%).

---

## Bulk RNA-seq: IMvigor210, GSE176307, UC-GENOME

**Immunotherapy-treated urothelial cancer cohorts with response labels.** The repo has no bulk script, so these were processed with stand-alone Python scripts (included in each upload with a README). No samples were removed; QC problems are flagged in the tables.

| Cohort | Source | Expression provided | Labelled (CR/PR vs SD/PD) |
| --- | --- | --- | --- |
| IMvigor210 | `IMvigor210CoreBiologies` `cds` object | raw counts (Entrez IDs; symbols + gene length in the object's rowData) | 68 vs 230 (298 of 348) |
| GSE176307 | GEO supplementary files | RSEM expected counts + Salmon TPM (symbols) | 16 vs 73 (89 of 90) |
| UC-GENOME | cBioPortal `blca_bcan_hcrn_2022` | **z-scores only** (no counts/TPM public) | 35 vs 59 (94 of 180) |

**What to watch for**

- GSE176307: expression columns are RNA-seq IDs (`RS-…`); match them to patients with `Sample_description` in the series matrix. Two TPM columns have no GEO sample. BACI165 has two samples. Five different PD-1/PD-L1 drugs and mixed lines of therapy.
- UC-GENOME: the cBioPortal datahub link is `https://datahub.assets.cbioportal.org/blca_bcan_hcrn_2022.tar.gz`. Because only z-scores exist, combine cohorts by z-scoring each gene within IMvigor210 and GSE176307 too (16k shared genes). 4 samples are recommended for exclusion (≥ 35% genes missing, PCA outliers). 11 have survival-from-immunotherapy longer than total survival. ODbL licence.
- IMvigor210: only 195 of 348 samples are from the bladder (rest upper tract/metastases); platinum timing varies (161 before, 105 after, 82 unknown).
