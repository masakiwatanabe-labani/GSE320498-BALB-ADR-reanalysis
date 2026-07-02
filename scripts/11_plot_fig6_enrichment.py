#!/usr/bin/env python3
"""
Fig.6-style GSEA enrichment plots, built directly from the canonical fgsea run
(ADR_B_vs_A, DESeq2 Wald stat ranking, A1-excluded main analysis).

Data source: outputs/tables/enrichment_curves/ (written by
outputs/scripts/10_export_enrichment_curves.R). The running-ES curve is the
exact numeric output of fgsea's own algorithm (reimplemented line-for-line from
fgsea::plotEnrichment / fgsea::calcGseaStat), not an approximation -- NES/FDR
annotated on each panel come from the identical fgsea() call that produced the
curve, so curve and statistics cannot drift apart.

Layout per panel (top to bottom), matching manuscript Fig.6B-D style:
  1. Green running enrichment score (ES) curve, dashed zero line.
  2. Black barcode of gene-set member hit positions in the ranked list.
  3. Red (higher in ByJcl) -> blue (higher in AJcl) gradient strip showing the
     ranking metric (DESeq2 Wald stat) across the full ranked gene list.
NES and joint FDR (q) are annotated on panel 1. x-axis is labelled with
"ByJcl (higher)" / "AJcl (higher)" at the two ends, per the project's fixed
sign convention (positive = higher in ByJcl) recorded in README_analysis_log.md.
"""

import os
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.colors import Normalize
import matplotlib.cm as cm
import numpy as np

BASE = "/usr/local/jupyter/ADR_BALB/outputs"
CURVE_DIR = os.path.join(BASE, "tables/enrichment_curves")
FIG_DIR = os.path.join(BASE, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

PANELS = {
    "ECM_ORGANIZATION": {
        "display_name": "Reactome: Extracellular Matrix Organization",
        "fig_label": "Fig. 6C (proposed replacement panel)",
        "outfile": "Fig6_enrichment_ECM_ORGANIZATION",
    },
    "INTEGRIN_SIGNALING": {
        "display_name": "Reactome: Integrin Signaling",
        "fig_label": "Fig. 6D",
        "outfile": "Fig6_enrichment_INTEGRIN_SIGNALING",
    },
    "KARAISKOS_TOP50": {
        "display_name": "Karaiskos et al. 2018 podocyte markers (TOP50)",
        "fig_label": "Fig. 6B (proposed replacement for podocyte-ageing set)",
        "outfile": "Fig6_enrichment_KARAISKOS_TOP50",
    },
}


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def load_metadata():
    rows = read_tsv(os.path.join(CURVE_DIR, "panel_metadata.tsv"))
    return {r["panel"]: r for r in rows}


def load_ranking():
    rows = read_tsv(os.path.join(CURVE_DIR, "ranking_vector_ADR_B_vs_A_stat.tsv"))
    stat = np.array([float(r["stat"]) for r in rows])
    return stat


def load_curve(panel):
    rows = read_tsv(os.path.join(CURVE_DIR, f"curve_{panel}.tsv"))
    x = np.array([float(r["x"]) for r in rows])
    y = np.array([float(r["y"]) for r in rows])
    return x, y


def load_hits(panel):
    rows = read_tsv(os.path.join(CURVE_DIR, f"hits_{panel}.tsv"))
    return np.array([int(r["hit_position"]) for r in rows])


def fmt_p(x):
    x = float(x)
    return f"{x:.2e}" if x < 0.001 else f"{x:.4f}"


def plot_panel(panel_key, meta, stat, fig_dir):
    info = PANELS[panel_key]
    m = meta[panel_key]
    x, y = load_curve(panel_key)
    hits = load_hits(panel_key)
    n = len(stat)

    fig = plt.figure(figsize=(8, 6))
    gs = gridspec.GridSpec(3, 1, height_ratios=[5, 0.6, 1.1], hspace=0.05)

    # --- Panel 1: running ES curve ---
    ax0 = fig.add_subplot(gs[0])
    ax0.plot(x, y, color="green", linewidth=2)
    ax0.axhline(0, color="black", linewidth=1)
    ax0.set_xlim(0, n)
    ax0.set_ylabel("Enrichment score (ES)", fontsize=13)
    ax0.tick_params(labelbottom=False, labelsize=11)
    ax0.set_title(f"{info['display_name']}\n{info['fig_label']}", fontsize=14, fontweight="bold")

    nes = float(m["NES"])
    pval = m["pval"]
    fdr = m["FDR_joint"]
    size = m["size"]
    annot = (f"NES = {nes:+.2f}\n"
             f"nominal P = {fmt_p(pval)}\n"
             f"joint FDR (q) = {fmt_p(fdr)}\n"
             f"gene set size = {size}")
    y_top = max(y.max(), 0)
    y_bot = min(y.min(), 0)
    if nes >= 0:
        ax0.text(0.98, 0.95, annot, transform=ax0.transAxes, ha="right", va="top",
                 fontsize=11, family="monospace",
                 bbox=dict(boxstyle="round", facecolor="white", alpha=0.85, edgecolor="gray"))
    else:
        ax0.text(0.02, 0.05, annot, transform=ax0.transAxes, ha="left", va="bottom",
                 fontsize=11, family="monospace",
                 bbox=dict(boxstyle="round", facecolor="white", alpha=0.85, edgecolor="gray"))
    pad = (y_top - y_bot) * 0.15 if (y_top - y_bot) > 0 else 0.05
    ax0.set_ylim(y_bot - pad, y_top + pad)
    for spine in ("top", "right"):
        ax0.spines[spine].set_visible(False)

    # --- Panel 2: black hit barcode ---
    ax1 = fig.add_subplot(gs[1], sharex=ax0)
    ax1.vlines(hits, 0, 1, color="black", linewidth=0.8)
    ax1.set_xlim(0, n)
    ax1.set_ylim(0, 1)
    ax1.set_yticks([])
    ax1.tick_params(labelbottom=False)
    for spine in ax1.spines.values():
        spine.set_visible(False)

    # --- Panel 3: red-blue ranking-metric gradient ---
    ax2 = fig.add_subplot(gs[2], sharex=ax0)
    # clip color scale at the 99th percentile of |stat| -- a few extreme-stat
    # genes would otherwise wash out the gradient across the other ~98% of the
    # ranked list (min/max are -15.7/+21.1 but 99th pct of |stat| is ~3.6)
    vmax = np.percentile(np.abs(stat), 99)
    norm = Normalize(vmin=-vmax, vmax=vmax)
    cmap = matplotlib.colormaps["RdBu_r"]
    img = np.tile(stat, (2, 1))
    ax2.imshow(img, aspect="auto", cmap=cmap, norm=norm, extent=[0, n, 0, 1])
    ax2.set_yticks([])
    ax2.tick_params(labelsize=11, labelbottom=True)
    for spine in ax2.spines.values():
        spine.set_visible(False)

    ax2.text(0.0, -1.1, "ByJcl (higher)", transform=ax2.transAxes,
              ha="left", va="top", fontsize=12, fontweight="bold", color="firebrick")
    ax2.text(1.0, -1.1, "AJcl (higher)", transform=ax2.transAxes,
              ha="right", va="top", fontsize=12, fontweight="bold", color="steelblue")
    ax2.text(0.5, -1.1, "Rank in ordered gene list\n(DESeq2 Wald statistic, Day 5 ADR ByJcl vs. AJcl)",
              transform=ax2.transAxes, ha="center", va="top", fontsize=11)

    fig.subplots_adjust(bottom=0.22)
    for ext in ("pdf", "svg", "png"):
        outpath = os.path.join(fig_dir, f"{info['outfile']}.{ext}")
        fig.savefig(outpath, dpi=300 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)
    return [os.path.join(fig_dir, f"{info['outfile']}.{ext}") for ext in ("pdf", "svg", "png")]


def main():
    meta = load_metadata()
    stat = load_ranking()
    written = []
    for panel_key in PANELS:
        written += plot_panel(panel_key, meta, stat, FIG_DIR)
    print("Wrote:")
    for w in written:
        print(" ", w)


if __name__ == "__main__":
    main()
