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

- `Manuscript_Results_and_Discussion_additions.md` — **the single consolidated
  document to work from when revising the manuscript text.** Updated
  2026-07-08 to be comprehensive: a "Quick reference" table at the top maps
  every reviewer comment to the exact section/paragraph to edit; includes the
  proposed Methods paragraph (M1), proposed Fig. 5/6A legend text (sign
  convention), 9 Results paragraphs (R1-R9: baseline DE description, Wt1/Nphs1,
  ADR-response concordance, GSEA FDR reporting, Fig. 6D reproduction, Fig. 6B
  sign-flip, Karaiskos panel substitute, R2-2 ECM gene validation, A1/A-Ctrl3
  QC), and 6 Discussion paragraphs (D1-D6: ranking-metric correction, "data
  not shown" contradiction, Karaiskos justification, bulk-RNA-seq caveat,
  sample-purity limitation, general limitations). Every number
  cross-references the source table/figure in `../`.
- `figures/Fig5B_Wt1_Nphs1_expression.{png,pdf}` (from
  `scripts/21_Wt1_Nphs1_expression_panel.R`, DESeq2-normalized counts, no new
  statistical model) — proposed Fig. 5B companion panel: per-sample
  boxplot+jitter of Wt1 and Nphs1 expression at baseline (AJcl-Ctrl vs.
  ByJcl-Ctrl only, matching Fig. 5's scope), so the "no baseline difference"
  claim is checkable against the raw expression distribution, not just a
  volcano-plot label. A-Ctrl3 (flagged for tubular contamination) is shown
  as an open triangle for transparency. Bracket annotates the main-analysis
  baseline_B_vs_A FDR for each gene.
- `tables/TableS_GSEA_all_comparisons_canonical.xlsx` (from
  `scripts/22_R1-6_GSEA_supplementary_table.R`, repackages the existing
  canonical `GSEA_<comparison>_full_joint_canonical.tsv` output, no new GSEA
  run) — R1-6 supplementary table: **every** gene set tested (1,293 of the
  1,336-set joint Reactome + podocyte universe survived fgsea's size filter),
  for all 4 primary comparisons, one sheet each, with nominal p AND
  BH-adjusted FDR/q-value columns plus an explicit `significant_FDR0.05`
  YES/no column (the manuscript's own stated significance criterion — not
  nominal p, which the original Fig. 6B/6D text used). Each sheet's caption
  states the exact software/package/versions (fgsea 1.24.0, msigdbr 26.1.0,
  MSigDB 2026.1.Mm, DESeq2 1.38.3), ranking metric, gene-set universe, and
  FDR method. An `overview` sheet summarizes gene sets tested/significant
  per comparison (219/1293, 355/1293, 99/1293, 171/1293).
- `figures/Figure6C_dotplot_fgsea.{png,pdf}` (+`_legend.txt`),
  `figures/Figure6C_alt_absNES.{png,pdf}` (+`_legend.txt`),
  `figures/Figure6C_alt_ECMonly.{png,pdf}` (+`_legend.txt`),
  `tables/TableS_fgsea_ADR_B_vs_A_full.csv` (from
  `scripts/23_Fig6C_fgsea_dotplot.R`, reuses the existing canonical
  `ADR_B_vs_A` A1-excluded fgsea output, no new GSEA run) — **Fig. 6C
  rebuilt from canonical fgsea**, replacing the original edgeR +
  ReactomePA `enrichPathway()` (ORA) dot plot (x=GeneRatio, size=Count,
  color=p.adjust), which mixed statistical frameworks with Fig. 6B/6D's
  GSEA-based panels (see `../README_analysis_log.md`, "A second, independent
  methodological correction" section, for that original ORA-vs-GSEA
  history). The new panel uses the same run as Fig. 6B/D
  (x=NES, size=gene-set size, color=joint BH FDR, y sorted by NES, dashed
  line at NES=0), in 3 variants: main = top 20 of 355 FDR<0.05 gene sets by
  FDR; alt = top 20 by \|NES\|; alt = all 22 ECM/collagen/integrin/adhesion
  name-matched gene sets regardless of significance (20/22 positive NES,
  15/22 FDR<0.05). Each `_legend.txt` states the selection rule, method, and
  axis/encoding definitions verbatim for drop-in use as a figure legend.
  `logs/23_Fig6C_consistency_check.txt` confirms all 4 primary comparisons
  (underlying Fig. 5, 6A, 6B, 6C, 6D) test the identical 19,662 genes and
  1,293 gene sets, with one consistent FDR definition throughout — the
  original manuscript's mismatched "Fig.5=13,064 vs. Fig.6A=13,278
  variables" issue does not recur — and reports that
  `REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION` (rank 55/1293, FDR=9.7e-11)
  and `REACTOME_INTEGRIN_SIGNALING` (rank 170/1293, FDR=9.1e-4) are both
  significant but fall outside the top-20-by-FDR view, outranked by large
  housekeeping gene sets (translation, ribosome, antigen presentation) —
  the reason the ECM/adhesion-filtered alt panel exists.
