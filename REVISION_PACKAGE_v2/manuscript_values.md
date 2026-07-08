# Manuscript values: canonical-run confirmed numbers, mapped to text location

Every "canonical run" value in this document is re-derived directly from
source tables by `scripts/24_manuscript_values_extraction.R` (Steps 1/2/3/5)
and `scripts/25_Fig6C_ORA_rerun_DESeq2_DEGs.R` (Step 4), run fresh for this
task — not copied from earlier summaries. Raw script output (with full
precision) is in `logs/24_Step1_Fig5_Fig6A_values.txt`,
`logs/24_Step2_Fig6B_values.txt`, `logs/24_Step3_Fig6D_values.txt`,
`logs/24_Step5_ADR_concordance_values.txt`, and
`logs/25_Fig6C_ORA_versions_and_summary.txt`. Where a value could not be
found, this is stated explicitly as "該当なし" rather than omitted or guessed.

**Canonical run definition (fixed):** DESeq2 v1.38.3, Wald test, main
analysis = A-ADR1 excluded; all contrasts = ByJcl vs AJcl (positive log2FC =
higher in ByJcl); GSEA = fgsea v1.24.0 multilevel, ranking = DESeq2 Wald
statistic; FDR = joint BH across all Reactome gene sets tested per
comparison. Fig. 6C ORA = the one new analysis in this task: ReactomePA
v1.42.0 `enrichPathway()`, DESeq2-derived DEG list, mouse organism.

## Reconciliation table

| # | 本文の記述箇所 | 現在の記載値 | canonical run の値 | 一致 or 要修正 |
|---|---|---|---|---|
| 1 | Fig. 5 legend, "total = N variables" | 13,064 variables | **19,662** tested genes (`baseline_B_vs_A`) | **要修正** — does not match either manuscript number |
| 2 | Fig. 6A legend, "total = N variables" | 13,278 variables | **19,662** tested genes (`ADR_B_vs_A`, A1-excluded) | **要修正** — does not match either manuscript number |
| 3 | (internal check) Fig. 5 vs Fig. 6A gene count | 13,064 ≠ 13,278 (mismatched between the two panels) | 19,662 = 19,662 (**identical**, same DESeq2 dds object for all 4 primary contrasts) | **要修正** — original mismatch is an edgeR-pipeline artifact; canonical run has no such mismatch |
| 4 | Prior draft response text ("19,662 tested genes") | — | **19,662** | **一致** — canonical run confirms this exact number |
| 5 | Results, baseline substrain divergence, "modest" | not quantified (qualitative "modest") | 103/19,662 genes padj<0.05 (0.5%); 64/19,662 padj<0.05 & \|log2FC\|>1 (0.3%) | **要修正** — replace qualitative "modest" with these counts (see `Manuscript_Results_and_Discussion_additions.md` R1) |
| 6 | Results, Wt1 baseline | not quantified ("no robust difference") | log2FC=+0.138, padj=0.869 | **要追加** — add exact values; direction of claim (non-significant) is correct |
| 7 | Results, Nphs1 baseline | not quantified | log2FC=+0.226, padj=0.794 | **要追加** — add exact values; direction of claim (non-significant) is correct |
| 8 | Results, Wt1 Day 5 | not quantified ("do not change robustly") | log2FC=−0.708, padj=**0.0498** (main, A1-excluded) | **要修正** — this is significant at the manuscript's own FDR<0.05 threshold, not merely "not robust"; report as borderline-significant, not absent |
| 9 | Results, Nphs1 Day 5 | not quantified ("do not change robustly") | log2FC=−0.781, padj=0.0638 (main, A1-excluded) | **要修正** — approaching but not reaching significance; describe as borderline, not simply "not robust" |
| 10 | (reference) Wt1/Nphs1 Day 5, A1 included | — | Wt1: log2FC=−0.611, padj=0.0683; Nphs1: log2FC=−0.628, padj=0.127 | **要追加** — optional footnote showing A1-inclusion sensitivity for these two genes |
| 11 | Results, Serpine1 Day 5 (R2-2) | direction stated qualitatively | log2FC=+1.840, padj=1.70×10⁻⁹ | **要追加** — add exact values; direction confirmed |
| 12 | Results, Loxl1 Day 5 (R2-2) | direction stated qualitatively | log2FC=+1.375, padj=4.86×10⁻¹² | **要追加** — add exact values; direction confirmed |
| 13 | Results, Col4a1 Day 5 (R2-2) | direction stated qualitatively | log2FC=+0.971, padj=4.48×10⁻³ | **要追加** — add exact values; direction confirmed |
| 14 | Results, Col4a2 Day 5 (R2-2) | direction stated qualitatively | log2FC=+1.216, padj=1.53×10⁻³ | **要追加** — add exact values; direction confirmed |
| 15 | Fig. 6B, podocyte-ageing set, `ADR_B_vs_A`, A1-excluded (manuscript's stated main analysis) | NES=−1.42 | **NES=+1.93**, nominal p=8.14×10⁻⁷, joint FDR=9.49×10⁻⁶ | **要修正 — OPPOSITE SIGN.** Canonical main-analysis configuration reproduces the manuscript's stated setting but with reversed direction |
| 16 | Fig. 6B, podocyte-ageing set, `ADR_B_vs_A`, A1 included | ("data not shown", "robustness...confirmed") | NES=**−2.04**, joint FDR=1.10×10⁻⁶ (matches manuscript SIGN only when A1 is put back in) | **要修正** — "robustness...confirmed, data not shown" is contradicted for this gene set; report actual numbers |
| 17 | Fig. 6B substitute: Karaiskos TOP50, `ADR_B_vs_A`, A1-excluded | not currently in manuscript | NES=−2.044, nominal p=6.21×10⁻⁵, joint FDR=5.25×10⁻⁴, size=49; leading edge includes Nphs1, Nphs2, Wt1, Podxl, Synpo, Mafb, Tcf21 | **新規提案** — proposed substitute/addition, see R7/D3 in `Manuscript_Results_and_Discussion_additions.md` |
| 18 | Fig. 6B substitute: Karaiskos EXCLUSIVE, `ADR_B_vs_A`, A1-excluded | not currently in manuscript | NES=−2.220, nominal p=9.48×10⁻⁶, joint FDR=9.29×10⁻⁵, size=12 | **新規提案** |
| 19 | Fig. 6D, Integrin Signaling, `ADR_B_vs_A`, A1-excluded | NES=+1.54, nominal p=0.0235 | **NES=+2.026**, nominal p=1.197×10⁻⁴, joint FDR=9.11×10⁻⁴ | **要修正** — same direction (positive/higher in ByJcl), but magnitude and nominal p both differ substantially from the manuscript-reported value; FDR (not previously reported) must also be added |
| 20 | (supplementary) Integrin Cell Surface Interactions, `ADR_B_vs_A`, A1-excluded | not in manuscript | NES=+2.167, nominal p=1.84×10⁻⁶, joint FDR=1.97×10⁻⁵, size=72 | **新規追加** — optional supporting gene set |
| 21 | (supplementary) ECM Organization, `ADR_B_vs_A`, A1-excluded | discussed qualitatively re: Fig. 6C | NES=+2.194, nominal p=4.12×10⁻¹², joint FDR=9.69×10⁻¹¹, size=248 | **要追加** — strongest and most robust focal result in the entire re-analysis |
| 22 | Fig. 6C, ORA method description | "Reactome pathway over-representation analysis of differentially expressed genes" | Method retained as ORA (ReactomePA::enrichPathway v1.42.0); **input DEG list re-derived from DESeq2** (`DE_ADR_B_vs_A.tsv`, padj<0.05, A1-excluded, n=546 DEGs; 543/546 mapped to Entrez, 99.5%) | **要修正** — method description can stay, but must state DESeq2 (not edgeR) as the DEG source, since the two engines' DEG lists are not identical gene-for-gene (see `../edgeR_vs_DESeq2_verification.md` for the ~50% Jaccard overlap at threshold in weaker-signal comparisons) |
| 23 | Fig. 6C, background/universe | not stated in manuscript | All 19,662 genes tested in `ADR_B_vs_A` (19,449/19,662 mapped to Entrez, 98.9%) | **要追加** — background was previously unstated; now explicit |
| 24 | Fig. 6C, number of significant terms | not stated in manuscript | 19 of 734 Reactome terms tested, BH-adjusted p<0.05 | **要追加** |
| 25 | Fig. 6C, top ORA term | not stated numerically | "Cell-extracellular matrix interactions"; 7/264 genes; BH-adjusted p=4.82×10⁻⁴ | **要追加** — citable sentence: *"Top term: Cell-extracellular matrix interactions; 7 genes; BH-adjusted p = 4.8 × 10⁻⁴."* |
| 26 | Fig. 6C, ECM/collagen terms present? | manuscript describes an ECM/collagen finding | **Yes** — 16 of 20 top terms are ECM/collagen/integrin/adhesion-related by name; "Extracellular matrix organization" itself is significant (18/264 genes, BH-adjusted p=0.0202); "Integrin cell surface interactions" is borderline (BH-adjusted p=0.0573, just above 0.05) | **一致** (direction of claim confirmed) — with the caveat that Integrin cell surface interactions narrowly misses BH<0.05 under ORA (it IS significant under GSEA, item 20) |
| 27 | Results, ADR-response substrain concordance (last item of R1-5) | not currently in manuscript | Pearson r=0.5354, Spearman ρ=0.4253 (n=19,662); direction concordance 62.8% genome-wide / 87.0% (sig. in either, n=2,715) / 99.6% (sig. in both, n=257); interaction term significant for 24/19,662 genes | **新規追加** — see R3 in `Manuscript_Results_and_Discussion_additions.md` |
| 28 | (QC caveat) Top 20 \|log2FC\| baseline genes | not currently in manuscript | Several top-|log2FC| genes have `padj = NA` (e.g. Trim75, Gm46851, Crabp2 — DESeq2 independent-filtering/low-count exclusions); see `tables/TableS_top20_absLog2FC_baseline.csv` | **要注意** — do not cite raw top-|log2FC| genes without checking padj is non-NA and significant; several of the largest fold changes are not statistically interpretable at low baseMean |

## Notes on items that could NOT be verified (report honestly, not silently dropped)

- The manuscript's **original numeric ORA output** for Fig. 6C (term list,
  GeneRatio, p.adjust values from the edgeR-derived DEG list) was not
  available to this re-analysis for a direct term-by-term diff — 該当なし.
  Item 26 above compares the *qualitative* claim (ECM/collagen enrichment)
  against the canonical ORA rerun, not a term-by-term numeric comparison,
  because the original ORA result table itself is not part of the deposited
  data.
- Fig. 6D's manuscript-stated **NES=+1.54** (item 19) was supplied as this
  task's input describing the current manuscript text; it was not
  independently re-derived from a primary source document by this script
  (the manuscript PDF/text itself is outside this re-analysis's inputs).

## Deliverables in this folder

- `figures/Figure6C_ORA_dotplot.{png,pdf}` — Fig. 6C ORA rerun on
  DESeq2-derived DEGs (x=GeneRatio, size=Count, color=BH-adjusted p; same
  semantics as the original panel)
- `tables/TableS_ORA_D5_full.csv` — full ORA result, all 734 Reactome terms
  tested (Description, GeneRatio, BgRatio, pvalue, p.adjust, Count, geneID)
- `figures/FigS_ADR_response_concordance.png` — copied in from
  `../figures/FigS_ADR_response_concordance.png` (produced earlier in this
  reviewer-response project); Step 5 statistics re-confirmed fresh in
  `logs/24_Step5_ADR_concordance_values.txt` and match this figure exactly
  (Pearson r=0.5354, Spearman ρ=0.4253)
- `logs/25_sessionInfo_ORA_and_canonical.txt` — full R `sessionInfo()`,
  package versions, seed (20260220)
- `logs/24_Step1_Fig5_Fig6A_values.txt`, `logs/24_Step2_Fig6B_values.txt`,
  `logs/24_Step3_Fig6D_values.txt`, `logs/24_Step5_ADR_concordance_values.txt`,
  `logs/25_Fig6C_ORA_versions_and_summary.txt` — full-precision raw output
  behind every value in the table above
- `tables/TableS_top20_absLog2FC_baseline.csv`,
  `tables/TableS_ADR_response_A_vs_B_merged.csv` — supporting tables
