# Bladder cancer project: analysis code

This repository holds the shared analysis code for the Ongkeko Lab bladder cancer project, which uses single-cell RNA-seq (scRNA-seq) and Visium spatial transcriptomics data.

> **Work in progress.** Code is added as the project moves forward, so expect folders, scripts and conventions to change. Each stage folder has its own README describing what's there and how to run it.

## Status

| Stage | Data | What it does | Status | Location |
| --- | --- | --- | --- | --- |
| 01 Preprocessing | scRNA-seq | Quality control: per-sample MAD filters, doublet removal (scDblFinder), ambient-RNA estimation (DecontX) | Available | [`01_preprocessing+metadata/scRNA/`](01_preprocessing+metadata/scRNA/) |
| 01 Preprocessing | Visium | Quality control of spots: per-sample MAD filters, maps of removed spots | Available | [`01_preprocessing+metadata/visium/`](01_preprocessing+metadata/visium/) |
| 01 Preprocessing | Visium | Deconvolution: cell-type proportions per spot (RCTD) | Available, waiting on an annotated bladder scRNA-seq reference | [`01_preprocessing+metadata/visium/`](01_preprocessing+metadata/visium/) |

Later stages (integration, annotation, downstream analysis) will be added here as they're written.

## Repository layout

Scripts and data are kept separate. Scripts are tracked in git. Data lives in `datasets/`, which is **git-ignored**: it never goes to GitHub, and each person keeps their own copy.

```
bladder cancer/                          repo root
├── README.md                            this file
├── bladder cancer.Rproj                 open this in RStudio
├── 01_preprocessing+metadata/           stage 1: QC and preprocessing
│   ├── README.md
│   ├── scRNA/
│   │   ├── scrna_qc.R                   scRNA-seq QC pipeline
│   │   ├── samplesheet_template.csv     example sample sheet
│   │   └── README.md                    how the QC works and how to run it
│   └── visium/
│       ├── h5ad_to_spaceranger.py       splits a combined Visium .h5ad into per-section folders
│       ├── visium_qc.R                  Visium spot QC
│       ├── visium_deconvolution.R       RCTD cell-type deconvolution
│       ├── samplesheet_template_visium.csv
│       └── README.md
└── datasets/                            NOT in git: your local data
    └── <dataset ID>/                    e.g. GSE130001/
        ├── samplesheet_<ID>.csv
        ├── data/                        raw downloads (one folder per sample)
        └── qc_output/                   written by the QC script
```

New stages get their own numbered folder (`02_...`, `03_...`), in the order they run.

## Getting started

1. **Clone the repo** and open `bladder cancer.Rproj` in RStudio. This sets the working directory to the repo root, which the scripts' default paths assume.
2. **Set up your data.** Create `datasets/<dataset ID>/` and put each sample's files in its own folder under `data/`, keeping the layout above.
3. **Write a sample sheet** for the dataset (start from `01_preprocessing+metadata/scRNA/samplesheet_template.csv`). Paths in a sample sheet are relative to the sheet's own folder, so the same sheet works on everyone's machine.
4. **Run the stage you need.** For example, scRNA-seq QC, run from the repo root:

   ```bash
   Rscript 01_preprocessing+metadata/scRNA/scrna_qc.R datasets/GSE130001/samplesheet_GSE130001.csv
   ```

   See [`01_preprocessing+metadata/scRNA/README.md`](01_preprocessing+metadata/scRNA/README.md) for installation, the sample sheet columns, the QC method and how to read the outputs.

## Datasets

| Dataset | Type | Samples | Notes |
| --- | --- | --- | --- |
| GSE130001 | scRNA-seq (10x, Cell Ranger v2 filtered) | 2 muscle-invasive bladder cancer tumors | QC run with `scrna_qc.R` |
| GSE171351 | Visium (Space Ranger 1.1, one combined `.h5ad`) | 4 treatment-naive muscle-invasive bladder cancer sections (Gouin et al.) | Split with `h5ad_to_spaceranger.py`; spots were already filtered by the authors |

Add a row when you start working with a new dataset, so everyone knows what's in use.

## Conventions

- **Never commit data or outputs.** Keep them under `datasets/`. If you need to share results, share them outside git (or ask before adding an exception to `.gitignore`).
- **No machine-specific paths in scripts.** Use paths relative to the repo root or to the sample sheet, so scripts run on anyone's copy.
- **One README per stage folder.** It should say what the scripts do, what they need as input, and what they produce.
- **Record what you ran.** Note the dataset, script and parameters for results you report, e.g. by keeping the `samplesheet_resolved.csv` the QC script writes.

## Contributing

1. Pull the latest `main` before starting.
2. Work on a branch for anything larger than a small fix.
3. Keep commits focused, with messages that say what changed.
4. Update the relevant README, and the status table above, when you add or change a script.
