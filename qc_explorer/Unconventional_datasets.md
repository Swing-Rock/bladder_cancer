# Unconventional datasets

Most GEO series load as they are: point the app (or a sample sheet) at a Cell Ranger / Space Ranger folder, an `.h5`, an `.rds` or an `.h5ad`. The series below are packaged differently. This page says how to load each one and what to watch for. If a new dataset doesn't fit any of these, add it here.

## Quick reference

| Dataset | Data type | What's unusual | What to do |
| --- | --- | --- | --- |
| [GSE149652](#gse149652) | scRNA-seq | Plain CSV count tables; sorted T cells only | Load the per-sample `*_count_matrice.csv.gz` files as `dense`; skip the malignant-cell step |
| [GSE171351](#gse171351) | Visium | Every section in one `.h5ad` | Pick the `.h5ad`; the app splits it into one folder per section |
| [GSE246011](#gse246011) | Visium | Loose files per sample (CSV counts + coordinates + TIFF), no Space Ranger folder; mixed cancer types | Put the bladder samples' files in their own folder and load it; the app converts them |

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
