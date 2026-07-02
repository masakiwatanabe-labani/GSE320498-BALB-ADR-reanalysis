#!/usr/bin/env python3
"""
Fig.6B replacement package: defensibility addenda requested on top of the
existing Karaiskos TOP50/EXCLUSIVE outputs (12_karaiskos_fig6b_panel_data.R,
13_plot_fig6b_karaiskos_panel.py):

1. Gene overlap between KARAISKOS2018_PODOCYTE_TOP50 and _EXCLUSIVE -- so the
   two-panel figure is described honestly (nested subset vs. independent
   sets). Appends columns to KARAISKOS_TOP50_vs_EXCLUSIVE_comparison.tsv.
2. Genome-wide rank (of ~1293-1296 jointly-tested Reactome+custom gene sets)
   for every focal gene set, in every comparison/A1-status configuration
   already run by 09_gsea_stat_ranking_canonical.R. Appends columns to
   GSEA_canonical_focal_judgment_table.tsv, backing Key finding 7 with exact
   numbers instead of a qualitative "not literally #1" statement.

Re-run after 09 and 12 (needs their full-joint and comparison TSVs on disk).
"""

import csv
import os

BASE = "/usr/local/jupyter/ADR_BALB/outputs"
TABLES = os.path.join(BASE, "tables")

KARAISKOS_EXCLUSIVE = ["Nphs2", "Cdkn1c", "Tcf21", "Enpep", "Nphs1", "Synpo",
                       "Npnt", "Wt1", "Pard3b", "Ptpro", "Iqgap2", "Mafb"]
KARAISKOS_TOP50 = ["Nphs2", "Cdkn1c", "Clic3", "Nupr1", "Dpp4", "Enpep", "Tcf21",
                    "Nphs1", "Gadd45a", "Rab3b", "Rhpn1", "Tmsb4x", "Col4a3",
                    "Rasl11a", "Mafb", "Npnt", "Arhgap24", "Adm", "Pak1", "Synpo",
                    "Foxd2os", "Golim4", "Igfbp7", "Vegfa", "Cd59a", "Sdc4",
                    "Sema3g", "Tdrd5", "Nap1l1", "Shisa3", "Eif3m", "Thsd7a",
                    "Pth1r", "Sept11", "Ctsl", "Podxl", "Cryab", "Mertk", "Htra1",
                    "Nes", "Wt1", "Npr3", "Ildr2", "Robo2", "Pard3b", "Tmem150c",
                    "Gas1", "Hoxc8", "Iqgap2", "Sema3e"]

FULL_JOINT_FILES = {
    ("baseline_B_vs_A", "A1_excluded_main"): "GSEA_baseline_B_vs_A_full_joint_canonical.tsv",
    ("ADR_B_vs_A", "A1_excluded_main"): "GSEA_ADR_B_vs_A_full_joint_canonical.tsv",
    ("A_ADR_vs_Ctrl", "A1_excluded_main"): "GSEA_A_ADR_vs_Ctrl_full_joint_canonical.tsv",
    ("B_ADR_vs_Ctrl", "A1_excluded_main"): "GSEA_B_ADR_vs_Ctrl_full_joint_canonical.tsv",
    ("ADR_B_vs_A", "A1_included_sens"): "GSEA_ADR_B_vs_A_full_joint_canonical_A1included.tsv",
    ("A_ADR_vs_Ctrl", "A1_included_sens"): "GSEA_A_ADR_vs_Ctrl_full_joint_canonical_A1included.tsv",
}
FOCAL_SETS = [
    "TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING",
    "REACTOME_INTEGRIN_SIGNALING",
    "REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS",
    "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
    "KARAISKOS2018_PODOCYTE_EXCLUSIVE",
    "KARAISKOS2018_PODOCYTE_TOP50",
]


def compute_overlap():
    ex, top = set(KARAISKOS_EXCLUSIVE), set(KARAISKOS_TOP50)
    overlap = ex & top
    n_overlap = len(overlap)
    pct = n_overlap / len(ex)
    print(f"EXCLUSIVE n={len(ex)}, TOP50 n={len(top)}, overlap n={n_overlap} "
          f"({pct:.1%} of EXCLUSIVE contained in TOP50); EXCLUSIVE-only: {sorted(ex - top)}")
    return n_overlap, pct


def append_overlap_to_comparison_table(n_overlap, pct):
    path = os.path.join(TABLES, "KARAISKOS_TOP50_vs_EXCLUSIVE_comparison.tsv")
    with open(path) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    fieldnames = [k for k in rows[0].keys()
                  if k not in ("n_overlap_with_other_set", "pct_EXCLUSIVE_contained_in_TOP50")]
    fieldnames += ["n_overlap_with_other_set", "pct_EXCLUSIVE_contained_in_TOP50"]
    for r in rows:
        r["n_overlap_with_other_set"] = n_overlap
        r["pct_EXCLUSIVE_contained_in_TOP50"] = f"{pct:.1%}"
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"Updated {path} with overlap columns")


def compute_genome_wide_ranks():
    rank_lookup = {}
    for (comp, a1), fn in FULL_JOINT_FILES.items():
        with open(os.path.join(TABLES, fn)) as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        n = len(rows)
        rows_sorted = sorted(rows, key=lambda r: float(r["padj"]) if r["padj"] not in ("NA", "") else 1.0)
        for i, r in enumerate(rows_sorted, start=1):
            if r["pathway"] in FOCAL_SETS:
                rank_lookup[(comp, a1, r["pathway"])] = (i, n)

    lookup_path = os.path.join(TABLES, "genome_wide_rank_lookup.tsv")
    with open(lookup_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["comparison", "a1_status", "pathway", "genome_wide_rank", "n_sets_tested"])
        for (comp, a1, pw), (rank, n) in rank_lookup.items():
            w.writerow([comp, a1, pw, rank, n])
    print(f"Wrote {lookup_path} ({len(rank_lookup)} rows)")
    return rank_lookup


def append_ranks_to_judgment_table(rank_lookup):
    path = os.path.join(TABLES, "GSEA_canonical_focal_judgment_table.tsv")
    with open(path) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    fieldnames = [k for k in rows[0].keys() if k not in ("genome_wide_rank", "n_sets_tested")]
    fieldnames += ["genome_wide_rank", "n_sets_tested"]
    for r in rows:
        key = (r["comparison"], r["a1_status"], r["pathway"])
        rank, n = rank_lookup.get(key, ("NA", "NA"))
        r["genome_wide_rank"] = rank
        r["n_sets_tested"] = n
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"Updated {path} with genome_wide_rank/n_sets_tested columns")


def main():
    n_overlap, pct = compute_overlap()
    append_overlap_to_comparison_table(n_overlap, pct)
    rank_lookup = compute_genome_wide_ranks()
    append_ranks_to_judgment_table(rank_lookup)


if __name__ == "__main__":
    main()
