# Why A-Ctrl3 is retained (not excluded), with numbers

A-Ctrl3 shows tubular contamination of a similar *magnitude* to the sample
that WAS excluded (A-ADR1), but the two samples are not treated the same way
because their effect on the analysis is different, not because of an
inconsistent standard. This file collects the exact numbers that justify
retaining A-Ctrl3. All values reused verbatim from
`../logs/19_ACtrl3_sensitivity_concordance.txt`,
`../logs/20_ACtrl3_sensitivity_figure_summary.txt`,
`../tables/Step19_ACtrl3_sensitivity_sig_counts.tsv`, and
`../tables/Step19_ACtrl3_sensitivity_key_genes.tsv` (from
`19_ACtrl3_sensitivity_dds_and_DE.R` / `20_ACtrl3_contamination_and_
sensitivity_figure.R`); see `manuscript_values.md`-style sourcing, no new
computation for this note.

## 1. The contamination is real and comparable in magnitude to A-ADR1 (this is NOT being denied)

- A-Ctrl3 podocyte:tubular ratio = 4.50, the lowest of the 6 baseline (Ctrl)
  samples (next-lowest = 8.52, 1.9x higher).
- A-Ctrl3's summed tubular-marker CPM (3483.9) is the 2nd-highest of all 12
  sequenced samples, behind only A-ADR1 (3739.5) — see
  `../figures/FigS_glomerular_purity_QC.png` (combined 12-sample figure) and
  `../figures/FigS_ACtrl3_contamination_QC.pdf`.

## 2. But excluding it does not change any conclusion — this is the actual reason it's retained

**Genome-wide log2FC concordance, baseline_B_vs_A, A-Ctrl3 included (main, n=3 vs 3) vs. excluded (sensitivity, n=2 vs 3):**
- Pearson r = **0.88**
- Direction concordance: **88.1%** genome-wide (n=19,573 genes compared)
- Direction concordance: **100%** among genes significant in either configuration (n=109)

**All 6 focal gene sets retain the same sign (0/6 sign flips):**

| Gene set | NES, A-Ctrl3 included (main) | NES, A-Ctrl3 excluded (sensitivity) |
|---|---|---|
| Podocyte ageing (Fig. 6B set) | −2.50 (FDR=6.0×10⁻¹⁴) | −2.03 (FDR=1.4×10⁻⁶) |
| Integrin signaling (Fig. 6D set) | +1.84 (FDR=0.022) | +1.72 (FDR=0.040) |
| Integrin cell surface interactions | +2.08 (FDR=5.1×10⁻⁴) | +1.91 (FDR=5.9×10⁻³) |
| ECM organization (Fig. 6C theme) | +2.28 (FDR=7.0×10⁻¹⁰) | +2.04 (FDR=2.0×10⁻⁷) |
| Karaiskos podocyte, exclusive (n=12) | +1.71 (FDR=0.058) | +1.57 (FDR=0.17) |
| Karaiskos podocyte, top50 (n=49) | +2.41 (FDR=6.3×10⁻⁶) | +1.79 (FDR=0.036) |

**Key genes discussed in the text (baseline_B_vs_A), A-Ctrl3 included vs. excluded:**

| Gene | log2FC, included (main) | padj, included | log2FC, excluded (sens) | padj, excluded |
|---|---|---|---|---|
| Tmem215 | +3.387 | 5.20×10⁻¹⁸ | +3.098 | 5.80×10⁻¹⁵ |
| Glp1r | +3.007 | 1.31×10⁻¹⁵ | +3.273 | 6.52×10⁻¹⁵ |
| Nlrp1b | −3.447 | 3.50×10⁻⁶ | −3.500 | 2.73×10⁻⁴ |
| Wt1 | +0.138 | 0.869 | +0.079 | 0.977 |
| Nphs1 | +0.226 | 0.794 | +0.176 | 0.935 |
| Loxl1 | +0.683 | 5.16×10⁻³ | +0.708 | 0.048 |
| Serpine1, Col4a1, Col4a2 | n.s. both ways | — | n.s. both ways | — |

Direction and significance-threshold conclusion is unchanged for every gene
in this table (Tmem215/Glp1r/Nlrp1b remain highly significant either way;
Wt1/Nphs1 remain non-significant either way).

## 3. Significant-gene counts drop but the pattern doesn't reverse

| Comparison | n padj<0.05, A-Ctrl3 included | n padj<0.05, A-Ctrl3 excluded |
|---|---|---|
| baseline_B_vs_A | 103 | 56 |
| A_ADR_vs_Ctrl | 298 | 241 |

This drop (n=3 vs 3 → n=2 vs 3) is the expected effect of losing one-third
of a group's statistical power, not evidence of a qualitatively different
result — nothing reverses sign or changes its significance-threshold
conclusion (see table above).

## 4. Why A-Ctrl3 and A-ADR1 are treated differently despite comparable contamination

- A-ADR1 was excluded because doing so **changes a conclusion**: the
  podocyte-ageing gene set's sign in `ADR_B_vs_A` depends on whether A-ADR1
  is included (see `manuscript_values.md` item 15-16; `Manuscript_Results_
  and_Discussion_additions.md` R6/D2).
- A-Ctrl3 is retained because excluding it **changes no conclusion** for
  `baseline_B_vs_A`: 0/6 focal gene sets flip sign, direction concordance is
  88-100%, and every specifically-discussed gene keeps the same
  significance-threshold call either way (sections 2-3 above).
- The exclusion rule applied consistently across this re-analysis is
  therefore **"exclude a flagged sample only if its inclusion changes a
  reported conclusion,"** not **"exclude any sample with elevated tubular
  contamination."** A-Ctrl3's contamination is real and disclosed (Panel A
  of `FigS_glomerular_purity_QC.png` and the dedicated
  `FigS_ACtrl3_contamination_QC.pdf`/`FigS_ACtrl3_sensitivity_NES_
  comparison.pdf`), but by this consistent rule it does not warrant
  exclusion from the main analysis.

## Suggested one-paragraph text (response letter / Methods note)

> A-Ctrl3 shows a podocyte:tubular marker-sum ratio (4.50) and summed
> tubular-marker CPM (3483.9) comparable in magnitude to the excluded
> A-ADR1 sample. We therefore repeated the baseline substrain comparison
> with A-Ctrl3 excluded as a sensitivity check. Genome-wide log2 fold-change
> estimates were highly concordant with the main analysis (Pearson r=0.88;
> 88.1% direction-concordant genome-wide, 100% among genes significant in
> either configuration), and all six focal gene sets retained the same
> direction of enrichment (0/6 sign flips). Because A-Ctrl3's exclusion
> does not alter any reported conclusion — unlike A-ADR1, whose exclusion
> does — A-Ctrl3 was retained in the main analysis (n=3 vs. 3 for the
> baseline comparison), consistent with excluding a flagged sample only
> when doing so changes a conclusion, not simply because contamination is
> present.
