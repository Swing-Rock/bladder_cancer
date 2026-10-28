#!/usr/bin/env python3
"""
h5ad_to_spaceranger.py — split a (combined) Scanpy/AnnData Visium .h5ad into one
Space-Ranger-style folder per tissue section, so visium_qc.R can read it.

Usage
    python h5ad_to_spaceranger.py <input.h5ad> <out_dir> [--sample-col COL] [--name-col COL]
                                  [--dataset NAME] [--sheet PATH]

    <input.h5ad>    AnnData file with spots in rows, raw counts (in .X, layers["counts"] or .raw),
                    obsm["spatial"] pixel coordinates and uns["spatial"][<library>] images + scalefactors
    <out_dir>       where the per-section folders are written (e.g. datasets/GSE171351/data/sections)
    --sample-col    obs column saying which section each spot belongs to; its values must be the
                    keys of uns["spatial"]. Default: detected automatically.
    --name-col      obs column used to name each section's folder (e.g. a patient ID); the library ID
                    is appended so names stay unique. Default: just the library ID.
    --dataset       value for the sample sheet's dataset column. Default: name of the input file.
    --sheet         where to write a ready-to-use sample sheet for visium_qc.R.
                    Default: <out_dir>/../../samplesheet_visium.csv (the dataset folder).

What each section folder gets (the layout visium_qc.R expects):
    <section>/filtered_feature_bc_matrix/matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz
    <section>/spatial/tissue_positions_list.csv, scalefactors_json.json,
                      tissue_lowres_image.png, tissue_hires_image.png

Requires: anndata, numpy, pandas, scipy, pillow
"""
import argparse
import gzip
import json
import os
import re
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
from PIL import Image


def die(msg):
    sys.exit(f"ERROR: {msg}")


def pick_counts(a):
    """Return (matrix, source) holding raw integer counts, spots x genes."""
    candidates = []
    if "counts" in a.layers:
        candidates.append((a.layers["counts"], 'layers["counts"]', a.var))
    if a.raw is not None:
        candidates.append((a.raw.X, ".raw.X", a.raw.var))
    candidates.append((a.X, ".X", a.var))
    for m, src, var in candidates:
        vals = m.data if sp.issparse(m) else np.asarray(m).ravel()
        vals = vals[: 2_000_000]
        if vals.size and np.all(np.mod(vals, 1) == 0) and vals.min() >= 0:
            return sp.csr_matrix(m), src, var
    die("no raw integer counts found in layers['counts'], .raw or .X (the file seems to hold "
        "normalized data only). QC and RCTD need raw counts.")


def detect_sample_col(a, libraries):
    for col in a.obs.columns:
        vals = set(a.obs[col].astype(str).unique())
        if vals and vals <= libraries:
            return col
    return None


def strip_concat_suffix(barcodes):
    """'AAACCGTTCGTCCAGG-1-0-0-0' -> 'AAACCGTTCGTCCAGG-1' (suffixes added when sections were concatenated)."""
    return [re.sub(r"^([ACGTN]+-\d+)(-\d+)+$", r"\1", b) for b in barcodes]


def to_uint8(img):
    img = np.asarray(img)
    if img.dtype != np.uint8:
        img = np.clip(img * (255 if img.max() <= 1.0 else 1), 0, 255).astype(np.uint8)
    if img.ndim == 3 and img.shape[2] == 4:
        img = img[:, :, :3]
    return img


def write_gz_lines(path, lines):
    with gzip.open(path, "wt") as fh:
        fh.write("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("h5ad")
    ap.add_argument("out_dir")
    ap.add_argument("--sample-col")
    ap.add_argument("--name-col")
    ap.add_argument("--dataset")
    ap.add_argument("--sheet")
    args = ap.parse_args()

    a = ad.read_h5ad(args.h5ad)
    print(f"Read {args.h5ad}: {a.n_obs} spots x {a.n_vars} genes")

    if "spatial" not in a.obsm:
        die('obsm["spatial"] (spot pixel coordinates) is missing; this is not usable as spatial data')
    libs = a.uns.get("spatial", {})
    if not libs:
        die('uns["spatial"] (images and scale factors) is missing')
    libraries = set(libs.keys())

    counts, src, var = pick_counts(a)
    print(f"Raw counts taken from {src}")

    sample_col = args.sample_col
    if sample_col is None:
        if len(libraries) == 1:
            sample_col = None
        else:
            sample_col = detect_sample_col(a, libraries)
            if sample_col is None:
                die(f"can't tell which obs column assigns spots to the libraries {sorted(libraries)}; "
                    f"pass --sample-col")
    if sample_col is not None:
        print(f"Section column: obs['{sample_col}']")
        lib_of_spot = a.obs[sample_col].astype(str).values
    else:
        lib_of_spot = np.array([next(iter(libraries))] * a.n_obs)

    gene_ids = var["gene_ids"].astype(str).values if "gene_ids" in var else np.asarray(var.index, dtype=str)
    symbols = np.asarray(var.index, dtype=str)
    ftype = var["feature_types"].astype(str).values if "feature_types" in var else np.array(["Gene Expression"] * len(symbols))
    features = [f"{i}\t{s}\t{t}" for i, s, t in zip(gene_ids, symbols, ftype)]

    os.makedirs(args.out_dir, exist_ok=True)
    dataset = args.dataset or re.sub(r"(_combined.*)?\.h5ad$", "", os.path.basename(args.h5ad))
    sheet = args.sheet or os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(args.out_dir))), "samplesheet_visium.csv")
    sheet_dir = os.path.dirname(os.path.abspath(sheet))
    rows = []

    for lib in sorted(libraries):
        m = lib_of_spot == lib
        n = int(m.sum())
        if n == 0:
            print(f"  {lib}: no spots, skipped")
            continue
        name = lib
        if args.name_col:
            vals = a.obs.loc[m, args.name_col].astype(str).unique()
            if len(vals) != 1:
                die(f"--name-col {args.name_col} has {len(vals)} values within library {lib}")
            name = f"{vals[0]}_{lib}"
        sec = os.path.join(args.out_dir, name)
        mdir = os.path.join(sec, "filtered_feature_bc_matrix")
        sdir = os.path.join(sec, "spatial")
        os.makedirs(mdir, exist_ok=True)
        os.makedirs(sdir, exist_ok=True)

        bcs = strip_concat_suffix(a.obs_names[m])
        if len(set(bcs)) != len(bcs):
            print(f"  {lib}: barcodes not unique after removing concatenation suffixes; keeping them as they are")
            bcs = list(a.obs_names[m])

        # counts: genes x spots, integers
        mat = counts[m].T.tocoo()
        mat = sp.coo_matrix((mat.data.astype(np.int64), (mat.row, mat.col)), shape=mat.shape)
        mtx_path = os.path.join(mdir, "matrix.mtx")
        scipy.io.mmwrite(mtx_path, mat, field="integer")
        with open(mtx_path, "rb") as fi, gzip.open(mtx_path + ".gz", "wb") as fo:
            fo.write(fi.read())
        os.remove(mtx_path)
        write_gz_lines(os.path.join(mdir, "barcodes.tsv.gz"), bcs)
        write_gz_lines(os.path.join(mdir, "features.tsv.gz"), features)

        # positions: Space Ranger v1 layout (no header):
        # barcode, in_tissue, array_row, array_col, pxl_row_in_fullres, pxl_col_in_fullres
        # obsm["spatial"] follows the Scanpy convention: column 0 = x (image column), column 1 = y (image row)
        xy = np.asarray(a.obsm["spatial"])[m]
        obs = a.obs.loc[m]
        pos = pd.DataFrame({
            "barcode": bcs,
            "in_tissue": obs["in_tissue"].astype(int).values if "in_tissue" in obs else 1,
            "array_row": obs["array_row"].astype(int).values if "array_row" in obs else 0,
            "array_col": obs["array_col"].astype(int).values if "array_col" in obs else 0,
            "pxl_row": np.round(xy[:, 1]).astype(int),
            "pxl_col": np.round(xy[:, 0]).astype(int),
        })
        pos.to_csv(os.path.join(sdir, "tissue_positions_list.csv"), header=False, index=False)

        # images + scale factors
        lib_uns = libs[lib]
        sf = {k: float(np.asarray(v)) for k, v in lib_uns.get("scalefactors", {}).items()}
        if "tissue_lowres_scalef" not in sf and "tissue_hires_scalef" not in sf:
            die(f"library {lib} has no tissue_lowres_scalef/tissue_hires_scalef in uns['spatial'][{lib}]['scalefactors']")
        with open(os.path.join(sdir, "scalefactors_json.json"), "w") as fh:
            json.dump(sf, fh, indent=2)
        imgs = lib_uns.get("images", {})
        if not imgs:
            die(f"library {lib} has no images in uns['spatial'][{lib}]['images']")
        for kind in ("lowres", "hires"):
            if kind in imgs:
                Image.fromarray(to_uint8(imgs[kind])).save(os.path.join(sdir, f"tissue_{kind}_image.png"))

        rows.append({"sample_id": name, "dataset": dataset,
                     "path": os.path.relpath(sec, sheet_dir).replace(os.sep, "/"),
                     "notes": f"from {os.path.basename(args.h5ad)}, library {lib}"})
        print(f"  {lib} -> {sec}: {n} spots, images: {', '.join(k for k in ('lowres', 'hires') if k in imgs)}")

    pd.DataFrame(rows).to_csv(sheet, index=False)
    print(f"Sample sheet for visium_qc.R: {sheet}")


if __name__ == "__main__":
    main()
