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
- **`manuscript_values.md`** — the single reconciliation table for embedding
  confirmed values into the manuscript text: every number needed for Fig. 5,
  6A, 6B, 6C, 6D and the ADR-response-concordance addition, mapped to its
  manuscript location, with current vs. canonical-run values and an explicit
  一致/要修正/新規追加 verdict for each row (28 rows). Two items stand out:
  (1) Fig. 6B's podocyte-ageing set has the OPPOSITE sign from the
  manuscript (NES=+1.93 vs. reported −1.42) in the manuscript's own stated
  main-analysis configuration; (2) Fig. 6D's Integrin Signaling reproduces
  in direction but not magnitude (NES=+2.03 vs. reported +1.54). Built from
  `scripts/24_manuscript_values_extraction.R` (Steps 1/2/3/5, re-derives
  every DE/GSEA value fresh from source tables, no new model) and
  `scripts/25_Fig6C_ORA_rerun_DESeq2_DEGs.R` (Step 4, the one new analysis).
- `figures/Figure6C_ORA_dotplot.{png,pdf}`, `tables/TableS_ORA_D5_full.csv`
  (from `scripts/25_Fig6C_ORA_rerun_DESeq2_DEGs.R`) — **Fig. 6C kept as ORA**
  (per the manuscript's stated method), but rerun with the DEG list swapped
  from edgeR- to DESeq2-derived (`DE_ADR_B_vs_A.tsv`, padj<0.05, A1-excluded,
  n=546 DEGs, 99.5% mapped to Entrez via `org.Mm.eg.db`; background = all
  19,662 tested genes). `ReactomePA::enrichPathway` v1.42.0, BH-adjusted:
  734 Reactome terms tested, 19 significant (BH p<0.05). Top term
  "Cell-extracellular matrix interactions" (7 genes, BH p=4.8e-4);
  "Extracellular matrix organization" itself significant (18 genes, BH
  p=0.020); 16/20 top terms are ECM/collagen/integrin/adhesion-related by
  name. This is a genuinely different, complementary panel to
  `Figure6C_dotplot_fgsea` above (ORA on a DEG list vs. GSEA on the full
  ranked list) — both are provided since the manuscript's stated method for
  Fig. 6C is specifically ORA.
- `tables/TableS_top20_absLog2FC_baseline.csv`,
  `tables/TableS_ADR_response_A_vs_B_merged.csv` — supporting tables behind
  `manuscript_values.md`.
- `logs/24_Step1_Fig5_Fig6A_values.txt`, `logs/24_Step2_Fig6B_values.txt`,
  `logs/24_Step3_Fig6D_values.txt`, `logs/24_Step5_ADR_concordance_values.txt`,
  `logs/25_Fig6C_ORA_versions_and_summary.txt`,
  `logs/25_sessionInfo_ORA_and_canonical.txt` — full-precision raw values
  and package/session versions behind every row of `manuscript_values.md`.
- `figures/FigS_ADR_response_concordance.png` — copied in from `../figures/`
  for self-containedness; Step 5 numbers re-derived fresh and match exactly.
- **`tables/Supplementary_Tables_S1-S11.xlsx`** — the submission-ready
  consolidated supplementary tables file (from
  `scripts/26_Supplementary_Tables_consolidated.R`, pure repackaging of
  already-verified tables, no new analysis except where noted). 22 sheets:
  a `TOC` sheet, then Table S1 (baseline top DE genes, 64 genes), **S2
  (substrain x treatment interaction genes, n=24, + Reactome/GO BP ORA on
  that list, 3 sub-sheets — table number fixed by the manuscript text)**,
  S3-S6 (full DE results for all 4 primary comparisons, 19,662 genes each,
  underlying Fig. 5/6A and both within-substrain ADR-response comparisons),
  S7 (canonical GSEA, all 1,293 tested gene sets x 4 comparisons, 5
  sub-sheets), S8 (Fig. 6C ORA full results, 734 Reactome terms), S9 (A-ADR1
  sensitivity, 3 sub-sheets), S10 (A-Ctrl3 sensitivity, 3 sub-sheets), S11
  (gene-level ADR-response concordance between substrains). Every sheet
  carries a caption row stating its contrast/sign convention, source, and
  significance criterion. This is the single file to attach to the
  manuscript submission for "supplementary tables." (Supersedes the earlier
  `Supplementary_Tables_S1-S10.xlsx`, removed, which had a table-number
  conflict with the manuscript text on "Table S2".)
- `scripts/27_interaction_gene_pathway_coherence.R`,
  `tables/TableS2_interaction_genes.csv`,
  `tables/TableS_interaction_ORA_reactome.csv`,
  `tables/TableS_interaction_ORA_GO_BP.csv`,
  `logs/27_interaction_pathway_coherence.txt`, `logs/27_verdict_summary.txt`
  — pathway-coherence check on the 24 substrain x treatment interaction
  genes (reuses the existing canonical interaction model verbatim, no new
  DE/interaction model fit; only new analysis = ORA on the fixed 24-gene
  list). **Verdict: (A, borderline).** Reactome ORA: 2/67 terms significant
  (BH<0.05), but both ("Smooth Muscle Contraction", "Striated Muscle
  Contraction") are driven by the identical 2-gene pair (Tpm2/Tpm4,
  paralogous tropomyosin genes) — a low bar for a 13-member term, not
  evidence of a shared program. GO BP: 0/339 significant. 22 of 24 genes
  (92%) belong to no significant term in either ontology. Only 1/24 genes
  (H2-Q6) matches a known-polymorphic-locus name pattern, and only 2/24
  show large pre-existing baseline divergence — not enough to suggest the
  interaction signal is dominated by allele/mapping artifacts. **Conclusion:
  supports keeping "rather than a globally divergent injury program" in the
  manuscript text, with the minor, explicitly-stated exception of the
  Tpm2/Tpm4 pair — not an unqualified "zero shared annotation" claim.**
- `figures/FigS_glomerular_purity_QC.{png,pdf}`,
  `tables/TableS_purity_ratios.csv` (from
  `scripts/28_glomerular_purity_QC_all12.R`, recomputes CPM directly from
  `../tables/00_merged_counts.tsv` using the exact same marker definitions
  as `18_A1_contamination_QC_figure.R`/`20_ACtrl3_contamination_and_
  sensitivity_figure.R`, no new sequencing analysis) — **single combined
  QC figure, all 12 samples in one panel**, superseding the previously
  separate ADR-only/Ctrl-only purity panels so A-ADR1's outlier status is
  visible against the full sample set at once. Podocyte:tubular marker-sum
  CPM ratio (log scale), colored by group (A-Ctrl/B-Ctrl/A-ADR/B-ADR), grey
  band = range of the other 11 samples, A-ADR1 highlighted in vermillion
  with a bold "A1 (excluded), ratio=2.51" label; A-Ctrl3 separately labeled
  "(retained)" so the figure cannot be read as implying it should also be
  dropped. All 12 recomputed ratios matched the prior scripts' values
  exactly (script halts if any mismatch is found — none were). A second
  panel is a 6-tubular-marker x 12-sample heatmap with A-ADR1's column
  outlined. One important correction versus the task's draft caption: the
  fold-difference to the *next-lowest sample across all 12* is **1.8x**
  (A-Ctrl3, 4.50), not 2.8x — the 2.8x figure that has been used elsewhere
  in this project is A-ADR1 vs. the next-lowest *within the ADR group only*
  (A-ADR2, 7.20); both are correct statements, but for different reference
  sets, and the combined 12-sample figure/caption uses the true
  all-samples comparison.
- **`ACtrl3_retention_justification.md`** — the specific numbers that justify
  *retaining* A-Ctrl3 despite tubular contamination comparable in magnitude
  to A-ADR1: genome-wide log2FC concordance (Pearson r=0.88), 0/6 focal
  gene sets flip sign with A-Ctrl3 excluded, and every specifically-discussed
  gene keeps the same significance-threshold call either way. States the
  applied rule explicitly: exclude a flagged sample only if its inclusion
  changes a reported conclusion (true for A-ADR1, false for A-Ctrl3), not
  "exclude any sample with elevated contamination." Includes a drop-in
  response-letter/Methods paragraph.
- `figures/Figure7_ECM_genes_CPM.{pdf,png}`, `tables/TableS_ECM_genes_CPM.csv`
  (from `scripts/29_Figure7_ECM_genes_CPM.R`, CPM computed directly from
  `../tables/00_merged_counts.tsv`, same definition as the other CPM-based
  QC figures; no prior "Figure7_ECM_genes_CPM" file existed anywhere in the
  project before this script — checked first) — CPM of Serpine1, Loxl1,
  Col4a1, Col4a2 (AJcl vs ByJcl, baseline and Day 5), 2x2 faceted panels,
  mean±SEM bars with individual samples overlaid, A-ADR1 shown as a
  separate orange triangle excluded from the Day-5 AJcl mean (n=2, matching
  the main analysis). PDF confirmed genuinely vector (`cairo_pdf` device;
  no embedded raster image; `pdfinfo` shows Producer=cairo, no `/Image`
  XObject). All 4 Day-5 log2FC values matched the manuscript text exactly
  (Serpine1 +1.84, Loxl1 +1.38, Col4a1 +0.97, Col4a2 +1.22). **Important
  correction to the task's assumption:** not all 4 genes are non-significant
  at baseline — Serpine1/Col4a1/Col4a2 are (FDR>0.45), but **Loxl1 is
  already significantly elevated at baseline** (log2FC=+0.68, FDR=0.0052),
  rising further by Day 5 (log2FC=+1.37); the figure annotates Loxl1's
  actual baseline FDR rather than labeling it "n.s.," and the caption states
  this explicitly rather than claiming a uniform "no baseline difference"
  across all 4 genes.
