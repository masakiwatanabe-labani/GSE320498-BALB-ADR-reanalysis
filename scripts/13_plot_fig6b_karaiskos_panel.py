#!/usr/bin/env python3
"""
Fig.6B replacement figure: KARAISKOS2018_PODOCYTE_TOP50 (main, proposed
replacement for the podocyte-ageing set) shown together with
KARAISKOS2018_PODOCYTE_EXCLUSIVE (adjacent panel, n=12 stringent/cell-type-
exclusive markers) -- demonstrating agreement between two independently
sized marker lists rather than dependence on one specific gene list.

Two comparisons, per the request:
  - ADR_B_vs_A   (A1-excluded, main manuscript setting; Fig.6-equivalent)
  - B_ADR_vs_Ctrl (within-ByJcl ADR response; strongest single result overall)

Data: outputs/tables/enrichment_curves/{curve,hits,ranking_vector}_*.tsv,
written by 12_karaiskos_fig6b_panel_data.R from the canonical (stat-ranked)
fgsea run -- NES/FDR annotated here are read from the same run, not
recomputed.
"""

import os
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.colors import Normalize
import numpy as np

BASE = "/usr/local/jupyter/ADR_BALB/outputs"
CURVE_DIR = os.path.join(BASE, "tables/enrichment_curves")
FIG_DIR = os.path.join(BASE, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

SET_LABELS = {
    "KARAISKOS_TOP50": "Karaiskos et al. 2018 podocyte markers, TOP50 (n=49 detected)",
    "KARAISKOS_EXCLUSIVE": "Karaiskos et al. 2018 podocyte markers, EXCLUSIVE (n=12)",
}
COMPARISON_LABELS = {
    "ADR_B_vs_A": "Day 5 ADR, ByJcl vs. AJcl (A1-excluded, main analysis)",
    "B_ADR_vs_Ctrl": "Within-ByJcl ADR response (ADR vs. Ctrl)",
}
STAT_AXIS_LABEL = {
    "ADR_B_vs_A": "DESeq2 Wald statistic, Day 5 ADR ByJcl vs. AJcl",
    "B_ADR_vs_Ctrl": "DESeq2 Wald statistic, ByJcl ADR vs. Ctrl",
}
END_LABELS = {
    "ADR_B_vs_A": ("ByJcl (higher)", "AJcl (higher)"),
    "B_ADR_vs_Ctrl": ("ADR (higher)", "Ctrl (higher)"),
}


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def load_curve(set_key, comp):
    rows = read_tsv(os.path.join(CURVE_DIR, f"curve_{set_key}_{comp}.tsv"))
    x = np.array([float(r["x"]) for r in rows])
    y = np.array([float(r["y"]) for r in rows])
    return x, y


def load_hits(set_key, comp):
    rows = read_tsv(os.path.join(CURVE_DIR, f"hits_{set_key}_{comp}.tsv"))
    return np.array([int(r["hit_position"]) for r in rows])


def load_ranking(comp):
    rows = read_tsv(os.path.join(CURVE_DIR, f"ranking_vector_{comp}_stat.tsv"))
    return np.array([float(r["stat"]) for r in rows])


def load_meta():
    rows = read_tsv(os.path.join(BASE, "tables/KARAISKOS_TOP50_vs_EXCLUSIVE_comparison.tsv"))
    return {(r["comparison"], r["gene_set"]): r for r in rows}


def fmt_p(x):
    x = float(x)
    return f"{x:.2e}" if x < 0.001 else f"{x:.4f}"


def draw_enrichment(ax0, ax1, ax2, set_key, comp, meta, title=None, title_fontsize=13):
    x, y = load_curve(set_key, comp)
    hits = load_hits(set_key, comp)
    stat = load_ranking(comp)
    n = len(stat)

    ax0.plot(x, y, color="green", linewidth=2)
    ax0.axhline(0, color="black", linewidth=1)
    ax0.set_xlim(0, n)
    ax0.tick_params(labelbottom=False, labelsize=10)
    for spine in ("top", "right"):
        ax0.spines[spine].set_visible(False)
    if title:
        ax0.set_title(title, fontsize=title_fontsize, fontweight="bold")

    key = (comp, f"KARAISKOS2018_PODOCYTE_{set_key.split('_', 1)[1]}")
    m = meta[key]
    nes = float(m["NES"])
    annot = (f"NES = {nes:+.2f}\n"
             f"P = {fmt_p(m['nominal_p'])}\n"
             f"joint FDR = {fmt_p(m['joint_FDR'])}\n"
             f"n = {m['n_genes']}")
    y_top = max(y.max(), 0)
    y_bot = min(y.min(), 0)
    if nes >= 0:
        ax0.text(0.97, 0.93, annot, transform=ax0.transAxes, ha="right", va="top",
                 fontsize=10, family="monospace",
                 bbox=dict(boxstyle="round", facecolor="white", alpha=0.85, edgecolor="gray"))
    else:
        ax0.text(0.03, 0.07, annot, transform=ax0.transAxes, ha="left", va="bottom",
                 fontsize=10, family="monospace",
                 bbox=dict(boxstyle="round", facecolor="white", alpha=0.85, edgecolor="gray"))
    pad = (y_top - y_bot) * 0.15 if (y_top - y_bot) > 0 else 0.05
    ax0.set_ylim(y_bot - pad, y_top + pad)

    ax1.vlines(hits, 0, 1, color="black", linewidth=0.8)
    ax1.set_xlim(0, n)
    ax1.set_ylim(0, 1)
    ax1.set_yticks([])
    ax1.tick_params(labelbottom=False)
    for spine in ax1.spines.values():
        spine.set_visible(False)

    vmax = np.percentile(np.abs(stat), 99)
    norm = Normalize(vmin=-vmax, vmax=vmax)
    cmap = matplotlib.colormaps["RdBu_r"]
    img = np.tile(stat, (2, 1))
    ax2.imshow(img, aspect="auto", cmap=cmap, norm=norm, extent=[0, n, 0, 1])
    ax2.set_yticks([])
    for spine in ax2.spines.values():
        spine.set_visible(False)
    left_lab, right_lab = END_LABELS[comp]
    left_col = "firebrick" if comp == "ADR_B_vs_A" else "darkorange"
    right_col = "steelblue" if comp == "ADR_B_vs_A" else "gray"
    ax2.text(0.0, -0.9, left_lab, transform=ax2.transAxes, ha="left", va="top",
              fontsize=10, fontweight="bold", color=left_col)
    ax2.text(1.0, -0.9, right_lab, transform=ax2.transAxes, ha="right", va="top",
              fontsize=10, fontweight="bold", color=right_col)


def single_panel_figure(set_key, comp, meta):
    fig = plt.figure(figsize=(8, 6))
    gs = gridspec.GridSpec(3, 1, height_ratios=[5, 0.6, 1.1], hspace=0.05)
    ax0 = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1], sharex=ax0)
    ax2 = fig.add_subplot(gs[2], sharex=ax0)
    title = f"{SET_LABELS[set_key]}\n{COMPARISON_LABELS[comp]}"
    draw_enrichment(ax0, ax1, ax2, set_key, comp, meta, title=title)
    ax0.set_ylabel("Enrichment score (ES)", fontsize=12)
    ax2.text(0.5, -0.9, f"Rank in ordered gene list\n({STAT_AXIS_LABEL[comp]})",
              transform=ax2.transAxes, ha="center", va="top", fontsize=10)
    fig.subplots_adjust(bottom=0.24)
    outname = f"Fig6B_{set_key}_{comp}"
    paths = []
    for ext in ("pdf", "svg", "png"):
        p = os.path.join(FIG_DIR, f"{outname}.{ext}")
        fig.savefig(p, dpi=300 if ext == "png" else None, bbox_inches="tight")
        paths.append(p)
    plt.close(fig)
    return paths


def composite_figure(comp, meta):
    """Main (TOP50, large) + adjacent (EXCLUSIVE, smaller) in one figure."""
    fig = plt.figure(figsize=(13, 6.5))
    outer = gridspec.GridSpec(1, 2, width_ratios=[1.3, 1], wspace=0.28)

    gs_main = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=outer[0],
                                                height_ratios=[5, 0.6, 1.1], hspace=0.05)
    ax0m = fig.add_subplot(gs_main[0])
    ax1m = fig.add_subplot(gs_main[1], sharex=ax0m)
    ax2m = fig.add_subplot(gs_main[2], sharex=ax0m)
    draw_enrichment(ax0m, ax1m, ax2m, "KARAISKOS_TOP50", comp, meta,
                     title="Main: Karaiskos TOP50 (n=49)\nfull podocyte-cluster marker panel", title_fontsize=13)
    ax0m.set_ylabel("Enrichment score (ES)", fontsize=12)

    gs_adj = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=outer[1],
                                               height_ratios=[5, 0.6, 1.1], hspace=0.05)
    ax0a = fig.add_subplot(gs_adj[0])
    ax1a = fig.add_subplot(gs_adj[1], sharex=ax0a)
    ax2a = fig.add_subplot(gs_adj[2], sharex=ax0a)
    draw_enrichment(ax0a, ax1a, ax2a, "KARAISKOS_EXCLUSIVE", comp, meta,
                     title="Confirmatory: Karaiskos EXCLUSIVE (n=12)\n91.7% nested within TOP50", title_fontsize=13)

    fig.suptitle(f"Fig. 6B (proposed): Karaiskos podocyte markers (main panel + stringent core subset), {COMPARISON_LABELS[comp]}",
                 fontsize=13.5, fontweight="bold", y=1.03)
    fig.text(0.5, 0.10, f"Rank in ordered gene list ({STAT_AXIS_LABEL[comp]})",
             ha="center", va="top", fontsize=11)
    fig.subplots_adjust(bottom=0.26, top=0.86)
    outname = f"Fig6B_composite_{comp}"
    paths = []
    for ext in ("pdf", "svg", "png"):
        p = os.path.join(FIG_DIR, f"{outname}.{ext}")
        fig.savefig(p, dpi=300 if ext == "png" else None, bbox_inches="tight")
        paths.append(p)
    plt.close(fig)
    return paths


def main():
    meta = load_meta()
    written = []
    for comp in ["ADR_B_vs_A", "B_ADR_vs_Ctrl"]:
        for set_key in ["KARAISKOS_TOP50", "KARAISKOS_EXCLUSIVE"]:
            written += single_panel_figure(set_key, comp, meta)
        written += composite_figure(comp, meta)
    print("Wrote:")
    for w in written:
        print(" ", w)


if __name__ == "__main__":
    main()
