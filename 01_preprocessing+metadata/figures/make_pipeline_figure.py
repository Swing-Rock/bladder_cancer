#!/usr/bin/env python3
"""
make_pipeline_figure.py — draws the stage-1 preprocessing overview figure.

Usage:  python make_pipeline_figure.py [out_dir]
Writes pipeline_overview.png (300 dpi), .pdf and .svg (text kept as editable text).
Requires: matplotlib. Edit the BOXES / ARROWS below to change the figure.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib import font_manager

out_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))

# Arial/Helvetica if available (journal standard), otherwise a clean sans-serif
available = {f.name for f in font_manager.fontManager.ttflist}
family = next((f for f in ["Arial", "Helvetica", "Inter", "DejaVu Sans"] if f in available), "sans-serif")
plt.rcParams.update({"font.family": family, "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})

INK, QUIET, LINE = "#222222", "#5F5F5F", "#9A9A9A"
SC = {"edge": "#2F6DB5", "fill": "#EAF1FA", "head": "#2F6DB5"}     # scRNA-seq: blue
VS = {"edge": "#C46A1B", "fill": "#FCF0E4", "head": "#C46A1B"}     # Visium: orange
OUT = {"edge": "#7A7A7A", "fill": "#F4F4F4"}                       # data products

W, H = 7.2, 5.6                                                     # inches (double-column width)
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 100); ax.set_ylim(-3, 78); ax.axis("off")

LX, RX, BW = 37, 78, 34          # lane centres and box width
LABEL_X = 1.5


def box(cx, cy, w, h, title, lines=(), style=SC, dashed=False, title_size=8.6, line_size=7.2, mono_title=True):
    ax.add_patch(FancyBboxPatch((cx - w / 2, cy - h / 2), w, h,
                                boxstyle="round,pad=0,rounding_size=1.2",
                                linewidth=1.0, edgecolor=style["edge"], facecolor=style["fill"],
                                linestyle=(0, (3, 2)) if dashed else "-", zorder=2))
    n = len(lines)
    top = cy + (n * 2.35) / 2 + 0.2
    ax.text(cx, top, title, ha="center", va="center", fontsize=title_size, color=INK,
            fontweight="bold", family="monospace" if mono_title else None, zorder=3)
    for i, ln in enumerate(lines):
        ax.text(cx, top - 2.6 - i * 2.35, ln, ha="center", va="center", fontsize=line_size, color=INK, zorder=3)


def pill(cx, cy, w, h, text, size=7.2):
    ax.add_patch(FancyBboxPatch((cx - w / 2, cy - h / 2), w, h, boxstyle="round,pad=0,rounding_size=%.2f" % (h / 2),
                                linewidth=0.8, edgecolor=OUT["edge"], facecolor=OUT["fill"], zorder=2))
    ax.text(cx, cy, text, ha="center", va="center", fontsize=size, color=INK, zorder=3)


def arrow(x1, y1, x2, y2, color=LINE, dashed=False, lw=1.1):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=8, linewidth=lw,
                                 color=color, linestyle=(0, (3, 2)) if dashed else "-", shrinkA=0, shrinkB=0, zorder=1))


def row_label(y, step, name):
    n = name.count("\n") + 1
    top = y + 1.6 + 1.1 * (n - 1)           # keeps the label block centred on the row
    ax.text(LABEL_X, top, step, ha="left", va="center", fontsize=7.2, color=QUIET, fontweight="bold")
    ax.text(LABEL_X, top - 1.9, name, ha="left", va="top", fontsize=7.6, color=INK, linespacing=1.2)


# ---- title and lane headers ----------------------------------------------------------
ax.text(LABEL_X, 75.5, "Stage 1 preprocessing: bladder cancer scRNA-seq and Visium", ha="left", va="center",
        fontsize=10, color=INK, fontweight="bold")
for cx, label, st in [(LX, "scRNA-seq", SC), (RX, "Visium spatial transcriptomics", VS)]:
    ax.plot([cx - BW / 2, cx + BW / 2], [70.2, 70.2], color=st["head"], linewidth=2.2, solid_capstyle="butt")
    ax.text(cx, 72.0, label, ha="center", va="center", fontsize=8.8, color=st["head"], fontweight="bold")

# ---- rows ------------------------------------------------------------------------------
Y_IN, Y_FMT, Y_QC, Y_OBJ, Y_STEP4, Y_RES = 64.0, 54.0, 40.0, 27.5, 16.5, 6.0
row_label(Y_IN, "INPUT", "Downloaded data")
row_label(Y_FMT, "STEP 1", "Format")
row_label(Y_QC, "STEP 2", "Quality control")
row_label(Y_OBJ, "OUTPUT", "QC'd samples")
row_label(Y_STEP4, "STEPS 3–4", "Reference &\ndeconvolution")
row_label(Y_RES, "RESULT", "Spot composition")

# inputs
box(LX, Y_IN, BW, 6.2, "10x folder · .h5 · dense table · .rds", style=SC, mono_title=False, title_size=7.6)
box(RX, Y_IN, BW, 6.2, "Space Ranger · GEO files · combined .h5ad", style=VS, mono_title=False, title_size=7.6)

# format
ax.text(LX + 1.2, Y_FMT, "read directly", ha="left", va="center", fontsize=7.0, color=QUIET, style="italic")
box(RX, Y_FMT, BW, 7.4, "h5ad_to_spaceranger.py", ["only for a combined .h5ad: one folder per section"],
    style=VS, dashed=True)

# QC
box(LX, Y_QC, BW, 12.6, "scrna_qc.R",
    ["per-sample MAD filters: UMIs, genes, % mito",
     "doublet removal: scDblFinder",
     "ambient RNA: DecontX (separate assay)"], style=SC)
box(RX, Y_QC, BW, 12.6, "visium_qc.R",
    ["per-sample MAD filters: UMIs, genes, % mito",
     "in-tissue spots only",
     "maps of removed spots on the tissue"], style=VS)

# per-sample objects
pill(LX, Y_OBJ, BW, 4.6, "<sample>_qc.rds per sample (raw counts)")
pill(RX, Y_OBJ, BW, 4.6, "<sample>_qc.rds per section (raw counts)")

# reference (planned) and deconvolution
box(LX, Y_STEP4, BW, 9.0, "Reference (planned)", ["merge · integrate · annotate cell types",
                                                  "one shared cell-type label"], style=SC, dashed=True, mono_title=False)
box(RX, Y_STEP4, BW, 9.0, "visium_deconvolution.R", ["RCTD, full mode", "cell-type proportions per spot"], style=VS)

# result
pill(RX, Y_RES, BW, 4.6, "cell-type proportions + dominant type per spot")

# ---- arrows ------------------------------------------------------------------------------
gap = 0.5
arrow(LX, Y_IN - 3.1 - gap, LX, Y_QC + 6.3 + gap)                            # scRNA input -> QC
arrow(RX, Y_IN - 3.1 - gap, RX, Y_FMT + 3.7 + gap, dashed=True)              # Visium input -> format
arrow(RX, Y_FMT - 3.7 - gap, RX, Y_QC + 6.3 + gap)                           # format -> QC
ax.annotate("", xy=(RX + BW / 2 + 1.6, Y_QC + 6.3 + gap), xytext=(RX + BW / 2 + 1.6, Y_IN - 3.1 - gap),
            arrowprops=dict(arrowstyle="-|>", color=LINE, lw=1.1, connectionstyle="arc3,rad=0", shrinkA=0, shrinkB=0,
                            mutation_scale=8))                               # Space Ranger / GEO bypass
ax.text(RX + BW / 2 + 2.4, Y_FMT, "Space Ranger /\nGEO files", ha="left", va="center", fontsize=6.4, color=QUIET,
        style="italic", linespacing=1.1)
arrow(LX, Y_QC - 6.3 - gap, LX, Y_OBJ + 2.3 + gap)
arrow(RX, Y_QC - 6.3 - gap, RX, Y_OBJ + 2.3 + gap)
arrow(LX, Y_OBJ - 2.3 - gap, LX, Y_STEP4 + 4.5 + gap, dashed=True)
arrow(RX, Y_OBJ - 2.3 - gap, RX, Y_STEP4 + 4.5 + gap)
arrow(LX + BW / 2 + gap, Y_STEP4, RX - BW / 2 - gap, Y_STEP4, color=SC["edge"], dashed=True, lw=1.3)
ax.text((LX + RX) / 2, Y_STEP4 + 1.6, "annotated", ha="center", va="bottom", fontsize=6.6, color=SC["edge"])
ax.text((LX + RX) / 2, Y_STEP4 - 1.6, "reference", ha="center", va="top", fontsize=6.6, color=SC["edge"])
arrow(RX, Y_STEP4 - 4.5 - gap, RX, Y_RES + 2.3 + gap)

# ---- legend (bottom, left lane, below the reference box) --------------------------------------
lx, ly = LX - BW / 2, -1.2
items = [("solid", "script"), ("dashed", "optional or planned step"), ("pill", "data product")]
for (kind, label), x in zip(items, [lx, lx + 13, lx + 41]):
    if kind == "pill":
        ax.add_patch(FancyBboxPatch((x, ly - 1.0), 4.5, 2.0, boxstyle="round,pad=0,rounding_size=1.0",
                                    linewidth=0.8, edgecolor=OUT["edge"], facecolor=OUT["fill"]))
    else:
        ax.add_patch(FancyBboxPatch((x, ly - 1.0), 4.5, 2.0, boxstyle="round,pad=0,rounding_size=0.5",
                                    linewidth=0.9, edgecolor=LINE, facecolor="white",
                                    linestyle=(0, (3, 2)) if kind == "dashed" else "-"))
    ax.text(x + 5.5, ly, label, ha="left", va="center", fontsize=6.6, color=QUIET)

os.makedirs(out_dir, exist_ok=True)
for ext in ("png", "pdf", "svg"):
    fig.savefig(os.path.join(out_dir, f"pipeline_overview.{ext}"), dpi=300 if ext == "png" else None,
                facecolor="white", bbox_inches="tight", pad_inches=0.08)
print("wrote pipeline_overview.png/.pdf/.svg to", out_dir)
