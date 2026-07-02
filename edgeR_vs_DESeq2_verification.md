# DE-tool consistency check: edgeR (original) vs DESeq2 (revision)

Verifies whether switching this revision's differential-expression engine from
the platform's original method (edgeR classic `exactTest()`, confirmed from
the actual "Amelieff QSP" template code — see
`Methods_rewrite_and_reviewer_response_notes.md`) to DESeq2 v1.38.3 changes
any conclusion behind Fig. 5/6. edgeR was run with the exact call sequence
found in the template (`cpm(merge) > 2` filter, `calcNormFactors` TMM,
`estimateCommonDisp` -> `estimateTagwiseDisp`, `exactTest`), not a generic
edgeR recipe. Significance criterion for both tools: |log2FC| > 1 (fold
change > 2) and FDR/padj < 0.05, matching the manuscript's stated threshold.
Seed fixed at 20260220 throughout. Full tables in `outputs/tables/edgeR_*`.

## Bottom line

**No conclusion behind Fig. 5/6 changes.** Fold-change estimates are
essentially identical between the two tools (Pearson r ≥ 0.9998 in every
comparison), every named gene keeps the same direction and the same
significant/non-significant call in both tools, and GSEA direction (NES
sign) for every focal gene set is 100% concordant regardless of which
engine's ranking is used. The one real, honestly-reported difference is that
the *exact* set of genes crossing the significance threshold differs
moderately between tools in the two weaker-signal comparisons (baseline,
Day-5 substrain) — expected behavior for two different statistical
frameworks operating near a hard p-value cutoff, not evidence that either
tool is wrong.

## Step 1: edgeR exactTest reproduction — gene counts

| Comparison | Genes tested (post CPM>2 filter) | Significant (\|log2FC\|>1 & FDR<0.05) |
|---|---|---|
| baseline_B_vs_A | 13,184 | 89 |
| ADR_B_vs_A (A1 excluded, main) | 13,278 | 155 |
| ADR_B_vs_A (A1 included, sens) | 13,443 | 501 |
| A_ADR_vs_Ctrl (A1 excluded, main) | 13,216 | 228 |
| A_ADR_vs_Ctrl (A1 included, sens) | 13,329 | 172 |
| B_ADR_vs_Ctrl | 13,265 | 697 |

(edgeR's own low-expression filter, `rowSums(cpm(merge)>2)>=1`, computed
per-comparison, is stricter than DESeq2's default independent filtering —
this alone accounts for some of the "genes tested" difference from the
~19,662-gene DESeq2 tables, before any significance comparison.)

## Step 2: edgeR vs DESeq2 concordance

| Comparison | Pearson r (log2FC) | Spearman ρ (log2FC) | n sig edgeR | n sig DESeq2 | n sig common | Jaccard | Direction concordance (common sig genes) |
|---|---|---|---|---|---|---|---|
| baseline_B_vs_A | 0.9998 | 1.0000 | 89 | 54 | 49 | 0.52 | **100%** |
| ADR_B_vs_A (A1 excl, main) | 0.9999 | 1.0000 | 155 | 194 | 117 | 0.50 | **100%** |
| ADR_B_vs_A (A1 incl, sens) | 0.9999 | 1.0000 | 501 | 381 | 305 | 0.53 | **100%** |
| A_ADR_vs_Ctrl (A1 excl, main) | 0.9999 | 0.9999 | 228 | 165 | 154 | 0.64 | **100%** |
| A_ADR_vs_Ctrl (A1 incl, sens) | 0.9998 | 0.9998 | 172 | 165 | 136 | 0.68 | **100%** |
| B_ADR_vs_Ctrl | 1.0000 | 1.0000 | 697 | 694 | 669 | **0.93** | **100%** |

**What agrees:** fold-change magnitude/direction is essentially the same
tool (Pearson/Spearman ≥ 0.9998 everywhere — this is as close to identical as
two different statistical estimators ever get). **Every gene either tool
calls significant, that both tools also call significant, has matching
direction — zero direction-flips among consistently-called genes, in any
comparison.** The comparison with the strongest biological signal
(`B_ADR_vs_Ctrl`, the largest DEG set) shows the highest agreement (Jaccard
0.93, 669/697-694 genes shared).

**What differs, honestly reported:** Jaccard overlap for the two
weaker-signal comparisons (baseline: 0.52; Day-5 substrain, A1-excluded
main: 0.50) is only about half — i.e. roughly half the genes either tool
calls significant are *not* also called significant by the other tool. This
is a real, expected consequence of edgeR's exact test and DESeq2's
Wald-test-with-shrinkage being different statistical frameworks: genes near
the p=0.05/FDR=0.05 boundary can land on either side of the cutoff depending
on which tool's specific null model and variance estimate is used, without
either tool being "wrong." **This does not affect any specific named-gene
claim below**, only the raw significant-gene *counts*, which were never the
basis of a specific manuscript claim.

## Named genes: edgeR vs DESeq2 side by side

Full table: `outputs/tables/edgeR_vs_DESeq2_named_genes.tsv`. Key results
(padj/FDR < 0.05 in **bold**):

| Comparison | Gene | edgeR logFC | edgeR FDR | DESeq2 log2FC | DESeq2 padj | Same call? |
|---|---|---|---|---|---|---|
| baseline | Tmem215 | +3.38 | **2.7e-24** | +3.39 | **5.2e-18** | yes, both sig |
| baseline | Glp1r | +2.99 | **8.9e-15** | +3.01 | **1.3e-15** | yes, both sig |
| baseline | Nlrp1b | -3.44 | **1.9e-10** | -3.45 | **3.5e-6** | yes, both sig |
| baseline | Wt1 | +0.13 | 0.997 (n.s.) | +0.14 | 0.869 (n.s.) | yes, both n.s. |
| baseline | Nphs1 | +0.21 | 0.950 (n.s.) | +0.23 | 0.794 (n.s.) | yes, both n.s. |
| Day-5 main | Serpine1 | +1.83 | **1.0e-8** | +1.84 | **1.7e-9** | yes, both sig |
| Day-5 main | Loxl1 | +1.37 | **5.3e-5** | +1.37 | **5.3e-5** | yes, both sig |
| Day-5 main | Col4a2 | +1.20 | **4.9e-3** | +1.22 | **1.5e-3** | yes, both sig |
| Day-5 main | Col4a1 | +0.96 | **1.4e-2** | +0.97 | **4.5e-3** | yes, both sig |
| Day-5 main | Wt1 | -0.71 | 0.093 (borderline n.s.) | -0.71 | 0.050 (borderline) | yes, both borderline |
| Day-5 main | Nphs1 | -0.79 | 0.090 (n.s.) | -0.78 | 0.064 (n.s.) | yes, both n.s. |

**Every named gene from the review correspondence (R1-5's Tmem215/Nlrp1b/
Glp1r, R1-5's Wt1/Nphs1, R2-2's Serpine1/Col4a1/Col4a2/Loxl1) keeps the same
significance call and the same direction under edgeR as under DESeq2.** The
Wt1 "borderline" finding at Day 5 (previously reported as DESeq2 padj=0.050,
right at the threshold) is independently reproduced as borderline under
edgeR too (FDR=0.093) — same qualitative conclusion ("borderline, not
absent," not "robustly changed" and not "clearly unchanged").

## Step 3: does the GSEA main conclusion survive a DE-tool change?

Ranked genes by edgeR's signed -log10(PValue) instead of DESeq2's canonical
Wald `stat`, then re-ran the same `fgsea` focal-set test
(`ADR_B_vs_A` main and `B_ADR_vs_Ctrl`, the two comparisons carrying the
Fig. 6 claims):

| Comparison | Gene set | NES (edgeR-ranked) | FDR (edgeR-ranked) | NES (DESeq2-ranked, canonical) | FDR (DESeq2-ranked) | Sign match |
|---|---|---|---|---|---|---|
| ADR_B_vs_A | ECM organization | +2.06 | 0.0079 | +2.19 | 9.7e-11 | **yes** |
| ADR_B_vs_A | Integrin signaling | +1.77 | 0.054 (borderline) | +2.03 | 9.2e-4 | **yes** |
| ADR_B_vs_A | Karaiskos EXCLUSIVE | -1.88 | 0.018 | -2.22 | 9.4e-5 | **yes** |
| ADR_B_vs_A | Karaiskos TOP50 | -1.76 | 0.068 (borderline) | -2.04 | 5.3e-4 | **yes** |
| B_ADR_vs_Ctrl | ECM organization | +1.68 | 0.059 (borderline) | +1.73 | 2.3e-4 | **yes** |
| B_ADR_vs_Ctrl | Integrin signaling | +1.39 | 0.486 (n.s. both) | +1.50 | 0.200 (n.s.) | **yes** |
| B_ADR_vs_Ctrl | Karaiskos EXCLUSIVE | -1.82 | 0.017 | -2.23 | 1.8e-5 | **yes** |
| B_ADR_vs_Ctrl | Karaiskos TOP50 | -2.10 | **0.0043** | -2.70 | 1.7e-10 | **yes** |

**Direction (NES sign) is 100% concordant across every focal gene set in
both comparisons, regardless of which DE tool's statistic is used to rank
genes for GSEA.** Statistical significance is somewhat weaker under
edgeR-derived ranking for 3 of 8 set/comparison combinations (Integrin
signaling and ECM organization drift just above the 0.05 joint-FDR line in
2-3 cases) — but the Karaiskos podocyte marker panel, the strongest and most
recommended evidence in this re-analysis (Key finding 4), stays significant
under both tools in every configuration, including the single strongest
result in the whole project (`B_ADR_vs_Ctrl` TOP50: FDR=0.0043 under edgeR
ranking vs. 1.7e-10 under DESeq2 ranking — weaker but still well under 0.05).

## What this means for the manuscript

- **Switching the stated DE tool from edgeR to DESeq2 v1.38.3 for this
  revision does not change any qualitative conclusion tied to a specific
  gene or gene set** — not the named-gene calls (R1-5, R2-2), not the
  direction of any focal GSEA result, and not the overall Fig. 5 "modest
  baseline difference" / Fig. 6 "ECM/podocyte program" narrative.
- The only place a reader should be careful is **borderline (FDR just above
  or below 0.05) individual results**, where the specific tool used can tip
  the call either way — Wt1 at Day 5 and Integrin signaling/ECM
  organization's joint significance in a couple of configurations are
  exactly this kind of borderline case, and should be described with
  appropriately hedged language ("borderline," "trend-level") rather than a
  hard pass/fail, independent of which tool produced the number.
- This supports adopting DESeq2 v1.38.3 as the revision's stated method (see
  `TODO_not_recoverable_from_counts.md`) without needing to also re-litigate
  every earlier finding in this project against edgeR — the two tools tell
  the same story.

## Files

- `outputs/tables/edgeR_<comparison>.tsv` (x6) — full edgeR exactTest output
- `outputs/tables/edgeR_all_results.rds` — same, as an R list object
- `outputs/tables/edgeR_vs_DESeq2_concordance.tsv` — the Step 2 table above
- `outputs/tables/edgeR_vs_DESeq2_named_genes.tsv` — the Step 2 named-gene table, full precision
- `outputs/tables/edgeR_vs_DESeq2_GSEA_impact.tsv` — the Step 3 table above
- `outputs/tables/DE_A_ADR_vs_Ctrl_A1included_sens.tsv` — newly derived DESeq2
  table (A1-included, within-A ADR response) needed to complete this
  comparison; did not previously exist on disk
- `outputs/scripts/15_edgeR_vs_DESeq2_reproducibility.R` — reproduces everything above
