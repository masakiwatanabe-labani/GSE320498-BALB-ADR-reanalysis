# REVISION_PACKAGE_v2

This folder holds everything produced **from this point forward** in the
reviewer-response re-analysis. The original `REVISION_PACKAGE/` (one level up)
is left exactly as-is — nothing has been moved out of it, and it remains the
complete record of everything done before this folder was created (R1-5 final
deliverables, the A1/A-Ctrl3 sensitivity checks, GSEA re-analysis, Fig. 6B
replacement, etc. — see `../README_analysis_log.md` and `../MANIFEST.md` for
that full history).

## Layout

Same convention as the parent package:

- `figures/` — new figures
- `tables/` — new tables
- `scripts/` — new scripts (numbering continues from where `../scripts/`
  left off, i.e. starting at `21`)
- `logs/` — new run logs

Each new deliverable will be listed in this file as it's added, so this stays
a single, short index instead of requiring a search through the parent
folder's much longer history.

## Contents added so far

- `Manuscript_Results_and_Discussion_additions.md` — draft Results and
  Discussion text for the manuscript itself (not the response letter),
  covering: baseline substrain DE description (replaces "modest" wording),
  Wt1/Nphs1 verification, cross-substrain ADR-response concordance, GSEA
  ranking-metric correction, Integrin/ECM robustness, the podocyte-ageing
  sign-flip and its non-robustness, the Karaiskos marker panel as a proposed
  Fig. 6B substitute, and the A1/A-Ctrl3 sample-purity QC findings. Every
  number cross-references the source table/figure in `../`.
- `figures/Fig5B_Wt1_Nphs1_expression.{png,pdf}` (from
  `scripts/21_Wt1_Nphs1_expression_panel.R`, DESeq2-normalized counts, no new
  statistical model) — proposed Fig. 5B companion panel: per-sample
  boxplot+jitter of Wt1 and Nphs1 expression at baseline (AJcl-Ctrl vs.
  ByJcl-Ctrl only, matching Fig. 5's scope), so the "no baseline difference"
  claim is checkable against the raw expression distribution, not just a
  volcano-plot label. A-Ctrl3 (flagged for tubular contamination) is shown
  as an open triangle for transparency. Bracket annotates the main-analysis
  baseline_B_vs_A FDR for each gene.
