# Proposed Methods / Results / Discussion additions and revisions

Consolidated draft text for the authors to adapt and insert directly into the
manuscript. Every number below is pulled from the re-analysis output in `../`
(the main `REVISION_PACKAGE` folder); see `../README_analysis_log.md` for the
full derivation and the exact script behind each value. Figure/table
references use the actual generated file names — final figure/table
numbering is left to the authors. Bracketed notes `[...]` are instructions,
not manuscript text.

## Quick reference: what to change, where

| Reviewer comment | Manuscript section to edit | Text below | Key supporting file(s) |
|---|---|---|---|
| R1-6 (name GSEA software/version, ranking metric) | Methods | M1 | `../logs/09_gsea_canonical_versions.txt` |
| R1-5 (state contrast/sign for Fig. 5 & 6A) | Fig. 5 / Fig. 6A legends | Fig. 5/6A legend text below | `Fig5_volcano_baseline.pdf`, `Fig6A_volcano_D5.pdf` |
| R1-5 ("modest" baseline divergence) | Results | R1 | `TableS_top_baseline_genes.csv` |
| R1-5 (Wt1/Nphs1 verifiable) | Results | R2 | `Fig5_volcano_baseline.pdf`, `Fig6A_volcano_D5.pdf`, `figures/Fig5B_Wt1_Nphs1_expression.png` |
| R1-5 (baseline→ADR response comparison) | Results | R3 | `FigS_ADR_response_concordance.png` |
| R1-6 (FDR not nominal p; full gene-set table) | Results + Methods | R4 | `tables/TableS_GSEA_all_comparisons_canonical.xlsx` |
| R1-6 / Fig. 6D reproduction | Results | R5 | `GSEA_canonical_focal_judgment_table.tsv` |
| R1-7 / R2-4 (A1 sensitivity, "data not shown") | Results + Discussion | R6, D2 | `FigS_A1_sensitivity_NES_comparison.pdf` |
| Fig. 6B substitute proposal | Results + Discussion + Fig. 6B legend | R7, D3, D4 | `Fig6B_composite_*.pdf` |
| R2-2 (Serpine1/Col4a1/Col4a2/Loxl1) | Results | R8 | `DE_ADR_B_vs_A.tsv` |
| New QC finding (A1 + A-Ctrl3 purity) | Results (Methods/QC) + Discussion | R9, D5 | `FigS_A1_contamination_QC.pdf`, `FigS_ACtrl3_contamination_QC.pdf` |
| Fig. 6C rebuild (ORA → fgsea, unify with 6B/6D) | Fig. 6C + Results + legend | R10 | `figures/Figure6C_dotplot_fgsea.pdf`, `Figure6C_alt_ECMonly.pdf` |
| — | Discussion (limitations) | D6 | — |

---

## METHODS

### M1. Differential expression and gene set enrichment analysis (proposed replacement paragraph)

[Replaces the current Methods description of the DE/GSEA pipeline. Directly
addresses R1-6's request to name the software/package/version used and to
correct the stated ranking metric.]

> Differentially expressed genes were identified with DESeq2 (v1.38.3) using
> a Wald test. For gene set enrichment analysis, all detected genes were
> ranked by the DESeq2 Wald test statistic (a signed, significance-weighted
> ranking that combines the direction and magnitude of the fold change with
> its estimation uncertainty). Preranked GSEA was performed using the
> `fgsea` R package (v1.24.0, multilevel algorithm, `eps = 0`,
> minimum/maximum gene-set size 5/500) against the Reactome collection
> (MSigDB M2:CP:REACTOME, db_version 2026.1.Mm, mouse-native gene symbols
> via `msigdbr` v26.1.0; 1,333 gene sets) together with the Tabula Muris
> Senis kidney podocyte-ageing signature [and, if the Fig. 6B substitution in
> R7/D3 is adopted: a mouse glomerular single-cell podocyte marker panel from
> Karaiskos et al., 2018]. False discovery rate (Benjamini-Hochberg q-value)
> was computed jointly across every gene set tested in a given comparison
> (1,336 gene sets), not on any single gene set or focal subset in isolation.
> A fixed random seed (20260220) was set immediately before each `fgsea()`
> call for reproducibility.

[Full rationale for the ranking-metric correction — including the internal
q<p inconsistency in the originally reported Fig. 6D statistics that first
motivated this reconstruction — is in
`../Methods_rewrite_and_reviewer_response_notes.md`, and is summarized for
Discussion in D1 below.]

---

## FIGURE LEGENDS (sign convention — addresses R1-5 item 1 directly)

[Insert into the Fig. 5 and Fig. 6A legends. This is the single edit that
resolves R1-5's first point: "clarify which contrast is plotted and which
direction is positive."]

> **Figure 5 (baseline).** Volcano plot of differential gene expression
> between untreated BALB/cByJcl and BALB/cAJcl (DESeq2 Wald test). Positive
> log2 fold change indicates higher expression in BALB/cByJcl relative to
> BALB/cAJcl. *Wt1* and *Nphs1* are labeled in every panel regardless of
> significance; genes discussed in the text (*Tmem215*, *Nlrp1b*, *Glp1r*)
> are additionally highlighted with their exact log2 fold change and FDR.

> **Figure 6A (Day 5 post-adriamycin).** Volcano plot of differential gene
> expression between BALB/cByJcl and BALB/cAJcl at Day 5 after adriamycin
> administration (DESeq2 Wald test). Positive log2 fold change indicates
> higher expression in BALB/cByJcl relative to BALB/cAJcl at this time point.
> *Wt1* and *Nphs1* are labeled as in Figure 5.

---

## RESULTS

### R1. Baseline transcriptional divergence between substrains (replaces the "modest divergence" paragraph)

> At baseline (untreated), differential expression analysis identified 103 of
> 19,662 tested genes (0.5%) as significantly different between BALB/cByJcl
> and BALB/cAJcl (DESeq2 Wald test, FDR < 0.05), of which 64 genes (0.3% of
> genes tested) additionally showed a fold-change greater than 2-fold
> (|log2FC| > 1; Table S1 [`TableS_top_baseline_genes.csv`]). While this
> represents a small fraction of the transcriptome, several of the
> differentially expressed genes reached very high statistical confidence and
> substantial effect sizes, and include genes of plausible biological
> relevance to renal disease susceptibility and drug response. Notably,
> *Tmem215* (log2FC = +3.39, FDR = 5.2 × 10⁻¹⁸) and *Glp1r* (log2FC = +3.01,
> FDR = 1.3 × 10⁻¹⁵) were both markedly higher in BALB/cByJcl, while *Nlrp1b*
> (log2FC = −3.45, FDR = 3.5 × 10⁻⁶) was markedly higher in BALB/cAJcl. *Glp1r*
> encodes the glucagon-like peptide-1 receptor, which has reported
> nephroprotective signaling roles in the glomerulus; *Nlrp1b* encodes an
> inflammasome sensor implicated in podocyte injury pathways; both are
> therefore plausible contributors to a substrain-specific baseline
> susceptibility difference rather than incidental transcriptional noise. We
> revise our characterization of the baseline substrain difference
> accordingly: it is real, affects a small but non-trivial subset of genes,
> and includes individual genes with strong, disease-relevant effect sizes,
> even though it does not extend to the core podocyte-identity program (see
> R2 below).

### R2. Podocyte-identity markers Wt1 and Nphs1 across substrains and time points

> At baseline, neither *Wt1* (log2FC = +0.14, FDR = 0.87) nor *Nphs1*
> (log2FC = +0.23, FDR = 0.79) differed significantly between substrains,
> consistent with equivalent baseline podocyte identity (Figure 5B
> [`figures/Fig5B_Wt1_Nphs1_expression.png`] shows the underlying
> per-sample normalized expression for both genes). At Day 5
> post-adriamycin, *Wt1* approached but did not cross the significance
> threshold (log2FC = −0.71, FDR = 0.050), and *Nphs1* remained non-significant
> but closer to threshold than at baseline (log2FC = −0.78, FDR = 0.064). We
> therefore describe the Day-5 result as a small, borderline reduction in
> ByJcl relative to AJcl rather than an absence of change; both genes are
> labeled directly on the volcano plots (Figure 5, Figure 6A) so this
> assessment can be verified at a glance.

### R3. The adriamycin-induced transcriptional response is concordant between substrains in direction and magnitude

> To assess whether AJcl and ByJcl mount a similar transcriptional response to
> adriamycin, we compared the ADR-vs-control log2 fold change for every gene,
> computed independently within each substrain. The two substrains' responses
> were significantly correlated genome-wide (Pearson r = 0.535, 95% CI
> 0.525–0.545; Spearman ρ = 0.425; n = 19,662 genes; Fig. S_concordance
> [`FigS_ADR_response_concordance.png`]). Direction of change agreed for 62.8%
> of all genes, rising to 87.0% among genes significant in either substrain
> and 99.6% among the 257 genes significant in both. An explicit
> substrain-by-treatment interaction test identified only 24 of 19,662 genes
> (0.1%) with a significant interaction term. We conclude that the two
> substrains mount an ADR response that is highly concordant in both direction
> and overall magnitude (ordinary least-squares slope = 0.62 relative to
> identity). We note that AJcl was tested with less statistical power (n = 2,
> after excluding the flagged low-purity sample; see R9) than ByJcl (n = 3),
> so the larger raw differentially-expressed-gene count in ByJcl (2,674 vs.
> 298 genes at FDR < 0.05) partly reflects this power asymmetry rather than a
> larger true biological response; the genome-wide correlation and interaction
> test above, rather than the raw DEG counts, are the appropriate basis for
> comparing response magnitude between substrains.

### R4. Gene set enrichment results are reported by FDR/q-value, for every gene set tested

[Directly addresses R1-6's second point: the manuscript's stated
significance criterion is FDR < 0.05, but the originally reported
podocyte/integrin GSEA results were described by nominal p-value. This
paragraph, plus the supplementary table it references, closes that gap.]

> Gene set enrichment analysis was performed with the `fgsea` R package
> (v1.24.0, preranked, multilevel algorithm, minimum/maximum gene set size
> 5/500) against the Reactome collection (MSigDB M2:CP:REACTOME, 1,333 gene
> sets; MSigDB release 2026.1.Mm) plus three podocyte-focused custom gene
> sets, testing all 1,336 gene sets jointly per comparison, with
> Benjamini-Hochberg FDR correction applied across the full joint universe
> (1,293 of 1,336 sets passed the minimum/maximum size filter and were
> actually tested per comparison). Consistent with the FDR < 0.05 criterion
> stated in our Methods, we report the BH-adjusted FDR/q-value — not the
> nominal p-value — for every gene set discussed in the text, and provide the
> complete result (nominal p, FDR, normalized enrichment score, and leading
> edge genes) for all 1,293 tested gene sets, for all four primary
> comparisons, in Table S_GSEA [`tables/TableS_GSEA_all_comparisons_canonical.xlsx`],
> rather than for the six podocyte/integrin/ECM gene sets discussed in the
> main text alone.

### R5. Reactome Integrin Signaling and Extracellular Matrix Organization are robustly enriched at Day 5

> Ranking by the DESeq2 Wald statistic, the Reactome Integrin Signaling gene
> set (Fig. 6D) was significantly enriched in ByJcl relative to AJcl at Day 5
> post-adriamycin (NES = +2.03, nominal P = 1.2 × 10⁻⁴, **FDR = 9.1 × 10⁻⁴**),
> closely matching the originally reported statistics. This result was robust
> to every configuration tested, including inclusion of the flagged
> low-purity A-ADR1 sample (NES = +2.33, FDR = 4.9 × 10⁻⁵ with A-ADR1
> included) and the larger, related Reactome Integrin Cell Surface
> Interactions gene set (NES = +2.17 to +2.38 across configurations,
> FDR ≤ 2.0 × 10⁻⁵ in all cases). The Extracellular Matrix Organization gene
> set (Fig. 6C's theme) was similarly robust and was in fact the single most
> consistently significant focal gene set in the entire re-analysis, reaching
> FDR < 0.05 in every comparison and A1 configuration tested, including at
> baseline (NES = +2.28, FDR = 7.0 × 10⁻¹⁰) and in both within-substrain
> ADR-response comparisons (NES = +1.66 to +1.73, FDR ≤ 4.6 × 10⁻³).

### R6. The podocyte-ageing gene set (Fig. 6B) is sensitive to inclusion of the flagged low-purity sample

> In contrast, the Tabula Muris Senis kidney podocyte-ageing gene set used for
> Fig. 6B showed a result that depended on which AJcl Day-5 samples were
> included. With the flagged low-purity sample (A-ADR1) excluded — the
> manuscript's stated main-analysis configuration — this gene set was
> enriched in the *opposite* direction from originally reported (NES = +1.93,
> nominal P = 8.1 × 10⁻⁷, **FDR = 9.5 × 10⁻⁶**, i.e. higher, not lower, in
> ByJcl). Including A-ADR1 reversed the sign again (NES = −2.04,
> FDR = 1.1 × 10⁻⁶), matching the originally reported direction only when the
> flagged sample is retained. Both configurations individually reach high
> statistical confidence (FDR < 10⁻⁵); the instability is in which
> configuration is used, not in a lack of statistical power in either one
> (Fig. S_A1sensitivity [`FigS_A1_sensitivity_NES_comparison.pdf`]). This was
> the only one of six focal gene sets tested for which A-ADR1 status changed
> the direction of a significant result; the Integrin Signaling and
> Extracellular Matrix Organization sets (R5) and the Karaiskos podocyte
> marker sets (R7) all retained their sign regardless of A-ADR1 status.

### R7. A curated single-cell-derived podocyte marker panel confirms loss of podocyte transcriptional identity independent of sample inclusion

> As an independent check of podocyte-identity loss that does not rely on the
> ageing-associated gene set above, we tested two podocyte marker panels
> derived from Karaiskos et al. (2018, *J Am Soc Nephrol* 29:2060–2068,
> GEO GSE111107), a single-cell RNA-seq study of the mouse glomerulus: a
> 12-gene panel of the authors' most stringent, cell-type-exclusive podocyte
> markers, and a 49-gene panel of the top glomerular-podocyte cluster markers
> from the same study (detected in our count matrix; 50 in the original
> table, of which 49 mapped after accounting for the *Sept11*/*Septin11*
> nomenclature change). Both panels showed significant loss of the podocyte
> transcriptional program in ByJcl relative to AJcl at Day 5
> (12-gene panel: NES = −2.22, FDR = 9.3 × 10⁻⁵; 49-gene panel: NES = −2.04,
> FDR = 5.2 × 10⁻⁴), and — unlike the ageing gene set in R6 — **retained the
> same direction when A-ADR1 was included** (12-gene panel: NES = −1.87,
> FDR = 0.020; 49-gene panel: NES = −1.87, FDR = 0.0071). The within-ByJcl
> ADR-response comparison gave the single strongest gene-set result in the
> entire re-analysis for the 49-gene panel (NES = −2.70, FDR =
> 2.1 × 10⁻¹⁰, rank 11 of 1,293 gene sets tested genome-wide). Because these
> two marker lists overlap substantially (11 of the 12 exclusive markers are
> also members of the 49-gene panel), they should be understood as a full
> marker panel and its stringent core rather than as fully independent
> replication; nonetheless, the concordant, A1-status-independent direction
> across both panels and across comparisons supports podocyte transcriptional
> identity loss in ByJcl at Day 5 as a robust finding of this dataset,
> independent of the ageing-signature instability described in R6.

### R8. Validation of ECM-related genes discussed in the text (R2-2)

> The four ECM-related genes specifically discussed in the text — *Serpine1*,
> *Col4a1*, *Col4a2*, and *Loxl1* — were all significantly higher in ByJcl
> relative to AJcl at Day 5 post-adriamycin (*Loxl1*: FDR = 4.9 × 10⁻¹²;
> *Serpine1*: FDR = 1.7 × 10⁻⁹; *Col4a2*: FDR = 1.5 × 10⁻³; *Col4a1*:
> FDR = 4.5 × 10⁻³; DESeq2 Wald test), directly confirming the direction
> described in the text at the individual-gene level.

### R9. Quality-control assessment of glomerular preparation purity across all sequenced samples

> Because one AJcl Day-5 sample (A-ADR1) had been flagged for reduced
> glomerular purity, we assessed relative tubular contamination across all 12
> sequenced samples using the summed CPM of six canonical proximal/distal
> tubular epithelial markers (*Lrp2*, *Slc34a1*, *Aqp1*, *Slc12a1*, *Umod*,
> *Aqp2*) relative to six podocyte markers (*Nphs1*, *Nphs2*, *Wt1*, *Podxl*,
> *Synpo*, *Nes*; Fig. S_A1QC [`FigS_A1_contamination_QC.pdf`]). A-ADR1 had
> the highest tubular-marker CPM among all six Day-5 ADR samples on 5 of 6
> individual markers, a summed tubular-marker CPM 3.09-fold higher than the
> mean of the other five ADR samples, and the lowest podocyte:tubular
> marker-sum ratio of any ADR sample (2.51, versus 7.14–35.22 for the other
> five); its podocyte-marker CPM sum was not reduced relative to the other
> samples, indicating that the low ratio reflects excess tubular
> contamination rather than podocyte depletion. Extending this same
> assessment to the six untreated (baseline) samples identified one further
> sample, A-Ctrl3, with a comparable pattern: the second-highest
> tubular-marker CPM sum of all 12 sequenced samples (behind only A-ADR1) and
> the lowest podocyte:tubular ratio among the six baseline samples (4.50,
> versus 8.52–50.58 for the other five; Fig. S_ACtrl3QC
> [`FigS_ACtrl3_contamination_QC.pdf`]). This sample was not previously
> flagged by manual review. We therefore repeated the baseline substrain
> comparison (Fig. 5) with A-Ctrl3 excluded as a sensitivity check (n = 2 vs.
> 3, compared with the reported n = 3 vs. 3). Genome-wide log2 fold-change
> estimates were highly concordant between the two configurations (Pearson
> r = 0.88; 88.1% direction-concordant genome-wide, 100% among genes
> significant in either configuration), and **all six focal gene sets tested
> retained the same direction of enrichment with A-Ctrl3 excluded**
> (Fig. S_ACtrl3sens [`FigS_ACtrl3_sensitivity_NES_comparison.pdf`]). We
> conclude that, in contrast to A-ADR1's effect on the Day-5 podocyte-ageing
> result (R6), the baseline substrain-difference conclusions reported here are
> not sensitive to A-Ctrl3's inclusion or exclusion, and no change to the
> reported baseline analysis is warranted; we report the QC finding and
> sensitivity check for transparency.

### R10. Figure 6C rebuilt from the canonical fgsea run (replaces the original over-representation-analysis panel)

[Replaces the Fig. 6C dot plot and its legend. The original panel was
produced by a different statistical framework (edgeR + ReactomePA
`enrichPathway()`, over-representation analysis) than Fig. 6B/6D (preranked
GSEA), so the four Fig. 6 panels mixed methods; this brings C into the same
framework as B/D.]

> Figure 6C was regenerated from the same canonical preranked GSEA run
> (DESeq2 Wald-stat ranking, `fgsea` v1.24.0, Reactome M2:CP:REACTOME +
> podocyte-focused custom gene sets, 1,336-set joint universe, 1,293 sets
> tested after size filtering) used for Fig. 6B and 6D, rather than the
> original over-representation analysis. Of the 1,293 Reactome/podocyte
> gene sets tested in the `ADR_B_vs_A` (Day 5, A-ADR1-excluded) comparison,
> 355 (27.5%) reached joint FDR < 0.05. The 20 most significant gene sets by
> FDR were dominated by large, low-dispersion housekeeping categories
> (translation, ribosome biogenesis, proteasome assembly, antigen
> presentation), a known property of significance-weighted ranking on
> high-count genes (see Key finding 7, `../README_analysis_log.md`) rather
> than a biological claim. `REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION`
> (NES = +2.19, FDR = 9.7 × 10⁻¹¹) and `REACTOME_INTEGRIN_SIGNALING`
> (NES = +2.03, FDR = 9.1 × 10⁻⁴) are both individually significant but rank
> 55th and 170th of 1,293 gene sets respectively, and so do not appear
> in a plain top-20-by-FDR view. We therefore also filtered directly on
> ECM/collagen/integrin/adhesion gene-set names: of 22 matched Reactome gene
> sets, 20 had positive NES (higher in ByJcl) and 15 reached joint
> FDR < 0.05, indicating broad, consistent enrichment of extracellular
> matrix and cell-adhesion biology at Day 5 post-adriamycin in ByJcl, even
> though individual ECM/integrin gene sets are not the single
> most-significant hits genome-wide. The complete result for all 1,293
> tested gene sets is provided in Table S_fgsea
> [`tables/TableS_fgsea_ADR_B_vs_A_full.csv`].

[Proposed Figure 6C legend:] *"Dot plot of Reactome gene set enrichment,
Day 5 post-adriamycin, BALB/cByJcl vs. BALB/cAJcl (preranked GSEA,
`fgsea` v1.24.0, DESeq2 Wald-statistic ranking, A-ADR1 excluded). x-axis:
normalized enrichment score (NES); positive = higher in BALB/cByJcl. Point
size: number of genes in the gene set. Point color: Benjamini-Hochberg FDR,
computed jointly across all 1,293 tested Reactome gene sets (not on this
panel's subset alone). Gene sets shown: the 20 gene sets with the smallest
joint FDR among the 355 gene sets reaching FDR < 0.05."*

---

## DISCUSSION

### D1. On the corrected ranking metric and what it changes

> Re-examination of the originally reported gene set enrichment statistics
> revealed an internal inconsistency: Fig. 6D reported a nominal P-value
> (0.0235) alongside a smaller FDR-corrected q-value (< 0.001), a pairing that
> cannot arise from standard Benjamini-Hochberg correction of a magnitude-only
> (fold-change) ranking, since corrected q-values can never be smaller than
> the nominal p-value of the same test. This pattern is instead characteristic
> of a significance-weighted ranking statistic. We confirmed this
> reconstruction directly: genes ranked by signed(−log₁₀ p-value) × sign(log2
> fold change) produced an essentially identical gene order to the DESeq2
> Wald test statistic (Spearman ρ = 1.0000 across all four primary
> comparisons), and re-running preranked GSEA with the Wald statistic as the
> ranking metric reproduced the originally reported Integrin Signaling result
> closely (R5), which a literal signed-log2-fold-change ranking did not. We
> therefore used the Wald statistic as the ranking metric throughout this
> re-analysis and have corrected the Methods text accordingly (M1). This
> correction has two consequences that pull in opposite directions: it
> strengthens confidence in the Integrin Signaling/ECM Organization results
> (R5) but reinstates, rather than resolves, the sign instability of the
> podocyte-ageing result (R6).

### D2. The podocyte-ageing signature is not a robust basis for Fig. 6B, and the reported "data not shown" sensitivity check does not hold for this specific gene set

> The manuscript states that a sensitivity analysis including the flagged
> low-purity sample "confirmed the robustness of the pathway-level findings
> (data not shown)." Under the corrected ranking metric, this statement holds
> for the Integrin Signaling and Extracellular Matrix Organization results
> (R5), which are robust in direction across A-ADR1 status, but does not hold
> for the podocyte-ageing gene set used in Fig. 6B: including A-ADR1 reverses
> the direction of this specific result (R6). Because the reported main
> analysis (A-ADR1 excluded) itself yields the opposite sign from the
> published NES, and both directions are independently highly significant, we
> do not think this result — as currently presented in Fig. 6B — can be
> presented as evidence either for or against a podocyte-ageing transcriptional
> signature difference between substrains at Day 5. We recommend either
> removing this specific panel or replacing it with a gene set whose result is
> stable to sample inclusion; we present the latter option in R7 and propose
> it as a direct substitute.

### D3. A curated podocyte-identity marker panel is a more defensible basis for the podocyte-loss claim than an ageing-associated signature

> The Tabula Muris Senis podocyte-ageing gene set combines canonical podocyte
> genes with a substantial fraction of ageing-, senescence-, and immune
> (MHC)-associated genes, which may explain its sensitivity to the specific
> sample composition of a small ADR-treated group. The Karaiskos podocyte
> marker panels (R7), by contrast, are derived directly from single-cell
> RNA-seq clustering of the mouse glomerulus and consist of genes selected
> specifically for podocyte-cluster identity, not age-association. Both
> panels reproduced the direction of podocyte-identity loss in ByJcl at Day 5
> regardless of A-ADR1 status, and the larger panel produced the single
> strongest gene-set enrichment result in this entire re-analysis in the
> within-ByJcl ADR-response comparison. We recommend the authors consider
> citing the Karaiskos panel alongside, or substituting it for, the ageing
> signature currently used for Fig. 6B. We note one limitation of this
> substitution for completeness: 11 of the 12 genes in the stringent panel are
> also members of the larger 49-gene panel, so the two panels should not be
> described as fully independent lines of evidence, only as a full marker set
> and its most stringent core.

### D4. Bulk RNA-seq enrichment of a podocyte marker panel reflects the podocyte transcriptional program as a whole, not cell number specifically

[Recommended addition to the Fig. 6B legend and/or this paragraph, to
pre-empt a predictable reviewer question about bulk RNA-seq and cell-type
proportion.]

> Karaiskos marker genes are podocyte-cluster-identifying genes derived from
> single-cell data; in glomerulus-enriched bulk RNA-seq, a negative
> enrichment score for this panel indicates relative attenuation of the
> podocyte transcriptional program overall, but does not on its own
> distinguish reduced podocyte number from reduced per-cell marker expression.
> This bulk transcriptomic result should be interpreted alongside, and is
> complementary to, the histological findings (WT1-positive cell counts per
> glomerulus, NPHS1 immunoreactivity), which separately address cell number
> and per-cell expression intensity.

### D5. Sample purity as a source of variability in glomerulus-enriched bulk RNA-seq

> Beyond the single previously flagged sample, extending the same
> tubular-contamination assessment to all 12 sequenced samples identified a
> second sample (A-Ctrl3, an untreated AJcl sample) with comparable elevated
> tubular marker expression, which had not previously been flagged (R9). This
> sample's inclusion or exclusion did not change any conclusion of the
> baseline substrain comparison in this dataset, so no revision to the
> reported baseline analysis is required. We raise it nonetheless because it
> indicates that variability in glomerular enrichment purity — despite the
> Dynabead-based enrichment protocol — is not confined to the single sample
> identified by the original authors, and because the baseline comparison's
> top differentially expressed pathways are dominated by
> xenobiotic-metabolism and mitochondrial pathways characteristic of proximal
> tubule tissue (amino acid catabolism, respiratory electron transport,
> mitochondrial translation; all FDR < 10⁻¹⁴, all lower in ByJcl). We recommend
> the authors report tubular-marker CPM ratios as a routine QC metric for all
> samples in glomerulus-enrichment RNA-seq experiments of this kind, both in
> this study and in future substrain comparisons, since sample-level variation
> in enrichment purity, rather than substrain biology per se, is a plausible
> partial contributor to this specific pathway signature.

### D6. Limitations

> This re-analysis has several limitations that should be stated alongside
> the findings above. First, the AJcl Day-5 ADR-treated group has an
> irreducibly small sample size (n = 2 after excluding the flagged low-purity
> sample), which limits statistical power for that arm specifically and
> contributes to the raw differentially-expressed-gene count asymmetry
> described in R3; the genome-wide correlation and interaction-term analyses
> presented here are intended to partially mitigate, but cannot fully remove,
> this power asymmetry. Second, the enrichment-method correction described in
> D1 (gene set testing method and ranking metric) was verified for
> preranked GSEA specifically; we did not separately re-verify an
> over-representation-analysis-based reconstruction of the original pipeline
> against the preranked GSEA results reported here, since the manuscript's own
> reported statistics (normalized enrichment scores) already presuppose a
> preranked GSEA method for the panels discussed above. Third, the
> differential-expression engine used throughout this revision (DESeq2) was
> independently cross-checked against edgeR's classic exact test on the same
> count data and produced concordant conclusions for every named gene and
> focal gene set (log2 fold-change Pearson r ≥ 0.9998 in every comparison; see
> `../edgeR_vs_DESeq2_verification.md`), but we note this as a verification
> exercise rather than a claim that the original manuscript's exact
> statistical pipeline has been fully reconstructed.
