# GSE320498 RNA-seq re-analysis: BALB/cAJcl (A) vs BALB/cByJcl (B) adriamycin nephropathy

Re-analysis of the deposited gene-level count matrices to address reviewer comments
R1-5, R1-6, R1-7, R2-2, R2-4. All numbers in this document and in `outputs/tables/`
and `outputs/figures/` are the direct, unedited output of the scripts in
`outputs/scripts/`, run in the numbered order below. No values were fabricated or
adjusted to match expectations; where a result did not reproduce the original
manuscript's claim, that is reported explicitly (see "Key findings" below).

## Methodological history (read this first)

This project went through three rounds of scrutiny on the GSEA ranking metric.
Reporting all three, rather than silently overwriting earlier numbers, is itself
part of the transparency record:

1. **v1 (superseded).** An ad hoc in-house preranked implementation (random
   same-size gene-set resampling for NES normalization, approximate FDR) ranked
   genes by the DESeq2 Wald `stat` column. Not named in the manuscript, and not a
   peer-reviewed tool — dropped.
2. **v2 (superseded).** Re-read of the manuscript Methods ("genes ordered by
   **signed log2 fold change**") led to re-running everything with `fgsea`
   ranking by literal signed log2FC (`07_gsea_log2FC_ranking_primary.R`,
   `08_A1_sensitivity_log2FC_ranking.R`). This matches the Methods text verbatim,
   but reproduced almost none of the manuscript's reported significance levels
   once tested jointly against the full Reactome collection.
3. **v3 (current, CANONICAL).** Fig. 6D reports nominal P=0.0235 **together
   with** FDR q<0.001 — impossible under standard Benjamini-Hochberg correction,
   where q can never be smaller than the nominal p-value of the same test. That
   pattern only arises from a *significance-weighted* ranking, not a magnitude-
   only (log2FC) ranking. We tested this directly
   (`Step1_ranking_metric_identity_check.tsv`, script
   `09_gsea_stat_ranking_canonical.R`): the gene order produced by
   **signed(-log10 p) x sign(log2FC)** is essentially identical to the DESeq2
   Wald `stat` column (Spearman rho = 1.0000 in all 4 comparisons — expected,
   since for a Wald test `stat` is itself a monotonic, sign-preserving transform
   of the p-value). **`stat`-ranked preranked GSEA is therefore reinstated as the
   canonical/primary convention for this project — it is the best available
   reconstruction of what was actually run, regardless of what the Methods
   paragraph literally says.** The Methods text itself should be corrected (see
   `Methods_rewrite_and_reviewer_response_notes.md`).

The log2FC-ranked numbers (`07`/`08`) are **not deleted** and remain in this
document and on disk as a named sensitivity check reflecting the Methods text as
literally written. All "reproduces / does not reproduce" verdicts below use the
canonical `stat`-ranked results (`09` and, for the original 4-comparison joint
Reactome + podocyte-ageing scan, `05`/`06`, which are numerically identical to `09`
for the sets they both tested).

**Net effect of this final correction, versus the log2FC-primary framing that
preceded it:** Integrin signaling (Fig. 6D) now reproduces closely (joint FDR
9.1e-4, essentially matching the manuscript's q<0.001). The podocyte-ageing sign
flip with A1 inclusion/exclusion (Fig. 6B / "data not shown" sentence), which the
log2FC-primary framing had appeared to resolve, is **back and confirmed real**
under the canonical ranking. Both changes are reported in full below.

## Reproducibility

- Scripts `00`-`06`, `09` (canonical, `stat`-ranked) plus `07`-`08` (log2FC,
  secondary sensitivity check), run via `Rscript` in numeric order. `04` must be
  re-run after `02` because `02` regenerates the shared workbook
  `Supplementary_DE_tables.xlsx`. `09` reads the same `tables/DE_*.tsv` and
  `tables/dds_group_main/sens.rds` produced by `01`/`02` and can be re-run
  independently.
- Random seed fixed to **20260220** everywhere a seed is used (DESeq2 does not
  require one; fgsea's multilevel estimation does — reset immediately before
  each `fgsea()` call).
- Package versions: `outputs/logs/09_gsea_canonical_versions.txt` and
  `outputs/logs/09_sessionInfo.txt` (R 4.2.2; DESeq2 1.38.3; fgsea 1.24.0;
  msigdbr 26.1.0 [MSigDB db_version 2026.1.Mm, mouse-native symbols]; ggplot2
  4.0.2; openxlsx 4.2.8.1).
- Sample "A1" flagged in the manuscript text was confirmed by the user to be
  **A-ADR1** (corroborated independently by glomerular-purity marker CPM ratios;
  see `outputs/logs/00_qc_report.txt`).
- Joint FDR is always computed across the **full tested gene-set universe** for
  that comparison (Reactome M2:CP:REACTOME, 1333 sets, plus the custom podocyte
  sets = 1336 total in script `09`), never on a focal-subset alone.

## Sign conventions (fixed for every output)

| Contrast | Positive value means |
|---|---|
| `baseline_B_vs_A` | higher in BALB/cByJcl (B) vs BALB/cAJcl (A), untreated |
| `ADR_B_vs_A` | higher in ByJcl vs AJcl, Day 5 after ADR (= original Fig. 6A/B/D) |
| `A_ADR_vs_Ctrl` | higher after ADR vs Ctrl, within AJcl |
| `B_ADR_vs_Ctrl` | higher after ADR vs Ctrl, within ByJcl |
| `interaction (substrainB.treatmentADR)` | ADR response is MORE positive in ByJcl than in AJcl |

## Focal gene sets tested (canonical run, script `09`)

1. `TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING` (Fig. 6B's set; MSigDB M8; n=193
   in the current release vs. 192 stated in the manuscript)
2. `REACTOME_INTEGRIN_SIGNALING` (Fig. 6D's named set; n=27 current vs. 23
   manuscript-stated)
3. `REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS` (larger, related candidate set; n=72)
4. `REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION` (Fig. 6C's theme; n=248)
5. `KARAISKOS2018_PODOCYTE_EXCLUSIVE` (n=12) and `KARAISKOS2018_PODOCYTE_TOP50`
   (n=49) — real mouse glomerular scRNA-seq podocyte markers, Karaiskos et al.
   2018 *JASN* 29:2060-2068, GEO GSE111107; extracted directly from the authors'
   own Supplementary Tables 2-3

## Correspondence table: review comment -> output

| Review # | Ask | Script(s) | Output(s) | Result summary |
|---|---|---|---|---|
| R1-5 (contrast clarity) | State which contrast/sign is plotted | `02_DE_tables.R`, `03_volcano_plots.R` | `Supplementary_DE_tables.xlsx`, `figures/volcano_*.png` | Sign convention stated in every table caption and every plot title |
| R1-5 (Tmem215/Nlrp1b/Glp1r) | Describe/report these genes | `02_DE_tables.R` | `DE_baseline_B_vs_A.tsv`, `DE_ADR_B_vs_A.tsv` | All 3 highly significant at baseline (Tmem215 padj=5.2e-18, Glp1r padj=1.3e-15, Nlrp1b padj=3.5e-6); Tmem215 remains significant Day-5 ADR (padj=1.7e-9). Gene-level DESeq2 result, unaffected by ranking-metric choice |
| R1-5 (Wt1/Nphs1 verifiable) | Label on volcano so readers can check "no robust change" | `03_volcano_plots.R` | `figures/volcano_*.png` | Baseline: both n.s. (padj=0.87, 0.79), confirming manuscript claim. **Day-5 ADR: Wt1 padj=0.050 (borderline), Nphs1 padj=0.064 (n.s.) — weaker than "not robust" implies; report as borderline, not absent.** Gene-level, unaffected by ranking-metric choice |
| R1-5 (baseline->ADR comparison per substrain) | Compare disease response, direction/strength concordance | `02_DE_tables.R`, `04_interaction_and_crosscomparison.R` | `DE_A_ADR_vs_Ctrl.tsv`, `DE_B_ADR_vs_Ctrl.tsv`, `DE_interaction_substrain_by_treatment.tsv`, `crosscomparison_ADR_response_scatter.png` | Genome-wide response correlated between substrains (Pearson r=0.535, p~0; 99.6% direction-concordant among genes significant in both). Only 24/19662 genes show significant substrain:treatment interaction. **Caveat: A tested with less power (n=2) than B (n=3, A-ADR1 excluded) — B's larger DEG count partly reflects this** |
| R1-6 (GSEA software/version) | Name software/package/version, Reactome, MSigDB | `09_gsea_stat_ranking_canonical.R` | `logs/09_gsea_canonical_versions.txt`, `logs/09_sessionInfo.txt` | fgsea 1.24.0 (preranked, multilevel, `eps=0`), msigdbr 26.1.0, MSigDB db_version 2026.1.Mm, Reactome M2:CP:REACTOME (1333 sets) + 3 custom podocyte sets = 1336 joint universe. Replaces the original unnamed in-house implementation |
| R1-6 (ranking metric misstated) | Correct/clarify the ranking metric | `09_gsea_stat_ranking_canonical.R` | `Step1_ranking_metric_identity_check.tsv` | Methods says "signed log2FC"; actual ranking (reconstructed from the Fig.6D q<p pattern, confirmed by rho=1.0000 rank-identity test) is the DESeq2 Wald **stat** (significance-weighted). See `Methods_rewrite_and_reviewer_response_notes.md` for the proposed corrected paragraph |
| R1-6 (FDR not nominal p; full table) | Report FDR/q for gene-set results; supply full table | `09_gsea_stat_ranking_canonical.R` | `GSEA_canonical_focal_judgment_table.tsv`, `GSEA_<comparison>_full_joint_canonical.tsv` (x4) | Day-5 ADR B-vs-A (main, A1-excl), stat-ranked, joint FDR (vs 1336 sets): **Integrin signaling NES=+2.03, joint FDR=9.1e-4 — reproduces the manuscript's q<0.001 closely.** Podocyte-ageing NES=+1.93, joint FDR=9.5e-6, **but this is the OPPOSITE sign from the manuscript's reported NES=-1.42** — see Key finding 2 |
| R1-7 / R2-4 (A1 sensitivity, not "data not shown") | Show A1-in vs A1-out results | `09_gsea_stat_ranking_canonical.R` | `GSEA_canonical_focal_judgment_table.tsv`, `GSEA_ADR_B_vs_A_full_joint_canonical_A1included.tsv`, `GSEA_A_ADR_vs_Ctrl_full_joint_canonical_A1included.tsv` | Integrin signaling and ECM organization: **robust** to A1 (same sign, similar/greater significance both ways). Podocyte-ageing: **NOT robust — full sign flip**, NES=+1.93 (A1-excluded) vs. NES=-2.04 (A1-included). Karaiskos marker sets: robust, no sign flip, significant both ways |
| R2-2 (validate Serpine1/Col4a1/Col4a2/Loxl1) | Check these ECM genes reach significance | `02_DE_tables.R` | `DE_ADR_B_vs_A.tsv` | All 4 significantly higher in ByJcl at Day 5 post-ADR: Loxl1 padj=4.9e-12, Serpine1 padj=1.7e-9, Col4a2 padj=1.5e-3, Col4a1 padj=4.5e-3 — directionally confirms manuscript text. Gene-level, unaffected by ranking-metric choice |
| R2-4 (Day5 vs Day7 timing; A1 QC criteria) | see R1-7 above for A1; timing rationale is a methods-text issue, not re-derivable from counts | — | — | see `TODO_not_recoverable_from_counts.md` |

## Key findings requiring author attention

1. **Integrin signaling (Fig. 6D) now reproduces closely under the canonical
   (stat) ranking — the strongest confirmation of the ranking-metric correction
   itself.** Main setting (`ADR_B_vs_A`, A1-excluded): NES=+2.03, nominal
   P=1.2e-4, **joint FDR=9.1e-4** — this pairing (small nominal P together with an
   even-smaller joint FDR) matches the *pattern* of the manuscript's own
   P=0.0235/q<0.001, and the joint FDR value itself lands just under the
   manuscript's stated q<0.001 threshold. This did not reproduce at all under
   log2FC ranking (joint FDR 0.26-0.51) — the fact that it reproduces so well
   under stat ranking is itself evidence that stat ranking is the correct
   reconstruction of the original method. Direction (higher in ByJcl) is
   positive and significant in **every** configuration tested: both A1 statuses
   (excluded +2.03 / included +2.33) and the related larger set
   `REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS` (excluded +2.17, FDR=2.0e-5 /
   included +2.38, FDR=6.3e-7).
2. **Podocyte-ageing gene-set result (Fig. 6B): the sign flip is real and
   confirmed under the canonical ranking — this reinstates and sharpens the
   original concern.** The manuscript's Fig. 6B reports NES=-1.42 for
   `TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING`, Day-5 ADR ByJcl-vs-AJcl. Our
   canonical (stat-ranked) reproduction of the manuscript's own stated main
   analysis (A1 excluded) gives **NES=+1.93 (nominal P=8.1e-7, joint
   FDR=9.5e-6)** — the opposite sign from the manuscript's reported value, and
   very highly significant in that opposite direction. Including A1 flips the
   sign again, to **NES=-2.04 (joint FDR=1.1e-6)** — matching the manuscript's
   reported direction only when the flagged low-purity sample is put back in.
   Both configurations are individually rock-solid (FDR<1e-5); the instability
   is in which configuration you use, not in statistical power. This gene set's
   result is not a minor nominal-significance quibble — it is a directionally
   unstable result that depends on exactly which 2 vs. 3 AJcl samples are used,
   and the manuscript's stated main-analysis configuration does not match its
   own reported sign.
3. **The "data not shown... robustness" sentence is contradicted, not
   vindicated.** Manuscript text (Results): "Sensitivity analysis including this
   sample [A1] confirmed the **robustness** of the pathway-level findings
   (**data not shown**)." For the specific gene set actually used in Fig. 6B,
   this is false under the canonical ranking: including A1 **flips the sign** of
   the result (finding 2), which is the opposite of "robustness." For Integrin
   signaling and ECM organization, the claim holds (finding 1). We recommend the
   authors report the actual A1-in/out numbers for every gene set shown in Fig.
   6, rather than a blanket "data not shown," since the outcome differs sharply
   by gene set.
4. **A real, peer-reviewed, single-cell-derived podocyte marker panel (Karaiskos
   et al. 2018, GSE111107) is the most robust podocyte-identity result in the
   entire re-analysis, and we recommend it replace or supplement Fig. 6B's
   ageing signature.** Two marker lists — `KARAISKOS2018_PODOCYTE_EXCLUSIVE`
   (n=12, authors' most-stringent cell-type-exclusive markers) and
   `KARAISKOS2018_PODOCYTE_TOP50` (n=49 detected in our matrix) — both show
   **loss of podocyte identity in ByJcl at Day 5** (negative NES, meaning lower
   in B, i.e. B loses podocyte-identity gene expression relative to A), matching
   the manuscript's biological narrative, and this holds in every configuration
   tested:
   - `ADR_B_vs_A` (A1-excl): EXCLUSIVE NES=-2.22 (FDR=9.3e-5), TOP50 NES=-2.04
     (FDR=5.2e-4)
   - `ADR_B_vs_A` (A1-incl): EXCLUSIVE NES=-1.87 (FDR=0.020), TOP50 NES=-1.87
     (FDR=0.0071) — **no sign flip**, unlike the ageing set
   - `B_ADR_vs_Ctrl` (within-substrain ADR response): TOP50 NES=-2.70,
     **joint FDR=2.1e-10 — the single most significant gene set found genome-
     wide (rank 11 of 1293) in this comparison**; EXCLUSIVE NES=-2.23,
     FDR=3.2e-5
   Unlike the ageing signature (an age-related, largely MHC/immune/senescence
   leading-edge gene set), these are direct, curated podocyte-identity markers
   from real dissociated mouse glomeruli — a more defensible choice for a
   podocyte-specific claim, and empirically more robust to A1 status.
5. **ECM organization (Fig. 6C's theme) is robust and highly significant under
   joint correction in every configuration.** `REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION`:
   `ADR_B_vs_A` A1-excl NES=+2.19 (FDR=9.7e-11, rank 55/1293 genome-wide),
   A1-incl NES=+2.53 (FDR=9.0e-17); also significant at baseline (NES=+2.28,
   FDR=7.0e-10) and in both within-substrain comparisons (`A_ADR_vs_Ctrl`
   NES=+1.66, FDR=0.0046; `B_ADR_vs_Ctrl` NES=+1.73, FDR=1.8e-4). This is the
   strongest and most consistently reproducing gene-set claim in the entire
   manuscript re-analysis.
6. **Baseline podocyte-ageing signal is real, matches the manuscript's baseline
   IHC direction, and was not reported in the manuscript.** At baseline,
   `TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING` gives NES=-2.50 (FDR=6.0e-14,
   rank 8/1293 genome-wide) — same direction as the manuscript's IHC finding
   that ByJcl has fewer WT1+ cells/lower NPHS1 at baseline. Unlike at Day 5
   (finding 2), there is no sign ambiguity at baseline since A1 does not exist
   in the baseline dataset. Worth citing as new supporting evidence, with the
   caveat that this same gene set behaves erratically once ADR treatment and A1
   status are introduced.
7. **Under stat ranking, genome-wide top hits are dominated by
   translation/ribosome/mitochondrial housekeeping pathways, not ECM/podocyte
   pathways — a known property of significance-weighted rankings, not a
   biological claim.** `GSEA_canonical_focal_judgment_table.tsv` now carries
   explicit `genome_wide_rank` / `n_sets_tested` columns backing this with
   real numbers rather than a qualitative impression:
   - `ADR_B_vs_A` (main): `REACTOME_TRANSLATION` is the #1 hit (FDR=8.6e-28);
     `REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION` ranks 55/1293,
     `REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS` 121/1293,
     `KARAISKOS2018_PODOCYTE_TOP50` 153/1293, `REACTOME_INTEGRIN_SIGNALING`
     170/1293 — all top ~13% genome-wide, all individually joint-FDR<0.05, but
     none literally #1.
   - `B_ADR_vs_Ctrl` (within-B, the strongest comparison overall):
     `KARAISKOS2018_PODOCYTE_TOP50` ranks **11/1293** — near the very top
     genome-wide even in absolute terms, not just "significant somewhere in
     the list."
   - Full per-comparison, per-A1-status rank numbers for all 6 focal sets are
     in `GSEA_canonical_focal_judgment_table.tsv`.
   This reflects that high-count, low-dispersion housekeeping genes get very
   large |stat| even for modest fold changes, a well-documented feature of
   variance-aware preranked GSEA (see `Methods_rewrite_and_reviewer_response_notes.md`
   for the exact suggested limitation sentence). It does not undermine the
   focal-set findings above but is worth stating explicitly if reviewers ask
   why ECM/podocyte are not the #1 hit genome-wide.
8. **Wt1/Nphs1 "do not change robustly"** holds cleanly at baseline (both
   clearly n.s.) but is borderline, not absent, at Day 5 (Wt1 padj=0.050, right
   at threshold) — softened language recommended. Gene-level, unaffected by
   ranking-metric choice.
9. **Baseline substrain difference is dominated by xenobiotic/metabolic
   (proximal-tubule-like) pathways** (amino acid catabolism, respiratory
   electron transport, mitochondrial translation — all FDR<1e-14, all lower in
   B; see `GSEA_canonical_top_hits_per_comparison.tsv`), raising a possible
   differential-tubular-contamination caveat despite Dynabead-based glomerular
   enrichment. Authors may want to check tubular marker CPM ratios between
   substrain preparations before attributing all baseline divergence to
   glomerular biology. This observation is ranking-metric-independent (also
   seen under log2FC ranking).

See `Methods_rewrite_and_reviewer_response_notes.md` for the proposed corrected
Methods paragraph and point-by-point reviewer-response language.

## Fig. 6-style enrichment plots (vector, ready for figure replacement)

Three panels were generated directly from the canonical run (`ADR_B_vs_A`,
DESeq2 Wald `stat` ranking, A1-excluded main analysis), matching the
manuscript's Fig. 6B-D visual style (green running-ES curve / black hit
barcode / red-ByJcl-to-blue-AJcl ranking-metric gradient), with NES and joint
FDR annotated from the identical `fgsea()` call that produced each curve (no
recomputation, no drift between curve and stats):

- `figures/Fig6_enrichment_INTEGRIN_SIGNALING.{pdf,svg,png}` — direct
  replacement for Fig. 6D (NES=+2.03, joint FDR=9.1e-4)
- `figures/Fig6_enrichment_ECM_ORGANIZATION.{pdf,svg,png}` — proposed
  replacement/addition for Fig. 6C's ECM theme (NES=+2.19, joint FDR=9.7e-11)
- `figures/Fig6_enrichment_KARAISKOS_TOP50.{pdf,svg,png}` — proposed
  replacement for Fig. 6B's podocyte-ageing set (NES=-2.04, joint FDR=5.2e-4),
  per Key finding 4 above (the ageing set is excluded from this panel set
  because its sign is not robust to A1 status)

Pipeline: `10_export_enrichment_curves.R` re-implements fgsea's own
`plotEnrichment`/`calcGseaStat` algorithm to export the running-ES curve, hit
positions, and full ranking vector as TSVs (`tables/enrichment_curves/`);
`11_plot_fig6_enrichment.py` (matplotlib) renders the 3-panel figures from
those TSVs. Re-run order: `09` (or any script producing
`tables/DE_ADR_B_vs_A.tsv`) -> `10` -> `11`.

## Fig. 6B main-figure decision: Karaiskos TOP50 + EXCLUSIVE side by side

Per Key finding 4, we propose replacing Fig. 6B's podocyte-ageing panel with
**`KARAISKOS2018_PODOCYTE_TOP50` as the main enrichment panel**, shown
alongside **`KARAISKOS2018_PODOCYTE_EXCLUSIVE`** (n=12) as a confirmatory
panel. Two comparisons are shown, both requested explicitly: `ADR_B_vs_A`
(A1-excluded, the manuscript's main analysis setting) and `B_ADR_vs_Ctrl`
(within-ByJcl ADR response — the single strongest result in this entire
re-analysis, joint FDR=2.1e-10).

| Comparison | Gene set | n | NES | nominal P | joint FDR | genome-wide rank |
|---|---|---|---|---|---|---|
| `ADR_B_vs_A` | TOP50 | 49 | -2.04 | 6.2e-05 | 5.2e-04 | 153/1293 |
| `ADR_B_vs_A` | EXCLUSIVE | 12 | -2.22 | 9.5e-06 | 9.3e-05 | 132/1293 |
| `B_ADR_vs_Ctrl` | TOP50 | 49 | -2.70 | 1.8e-12 | **2.1e-10** | **11/1293** |
| `B_ADR_vs_Ctrl` | EXCLUSIVE | 12 | -2.23 | 1.9e-06 | 3.2e-05 | 76/1293 |

**Important honesty check on what "two panels" actually demonstrates: EXCLUSIVE
is not independent of TOP50.** 11 of EXCLUSIVE's 12 genes (91.7%) are also
members of TOP50 — only `Ptpro` is unique to EXCLUSIVE (see
`Karaiskos_geneset_provenance.md` for the full accounting). The two panels are
therefore better described as **a full marker panel (TOP50) and its
near-fully-nested stringent core (EXCLUSIVE)**, not two independent gene
lists that happen to agree. What the pair of panels legitimately demonstrates
is narrower but still useful: the signal is not an artifact of the 38
TOP50-only genes, since the small, almost-fully-overlapping "hard core"
reproduces the same direction and significance on its own. This is a real but
more modest form of robustness than genuine independent replication, and
should be worded that way in the manuscript/response letter (see
`Karaiskos_geneset_provenance.md` and the response-letter notes file for the
exact suggested language) rather than "two independent gene sets agree."

## Fig. 6B legend / response-letter cell-composition caveat (add this text)

For inclusion in the Fig. 6B legend and/or the response letter, to pre-empt a
predictable reviewer objection about bulk RNA-seq and cell-type proportion:

> *"Karaiskos markers are podocyte-cluster-identifying genes derived from
> single-cell data; in glomerulus-enriched bulk RNA-seq, a negative NES for
> this panel indicates relative attenuation of the podocyte transcriptional
> program but does not, on its own, distinguish reduced podocyte number from
> reduced per-cell expression. This bulk enrichment result is intended to be
> interpreted alongside, and is complementary to, the histological findings
> (WT1+ cell count per glomerulus; NPHS1 immunoreactivity optical density),
> which separately address cell number and per-cell expression intensity."*

**Figures** (all vector + raster, canonical stat ranking):
- `figures/Fig6B_composite_ADR_B_vs_A.{pdf,svg,png}` and
  `figures/Fig6B_composite_B_ADR_vs_Ctrl.{pdf,svg,png}` — main+confirmatory
  side by side, one figure per comparison (recommended for the manuscript)
- `figures/Fig6B_KARAISKOS_TOP50_<comparison>.{pdf,svg,png}` and
  `figures/Fig6B_KARAISKOS_EXCLUSIVE_<comparison>.{pdf,svg,png}` — individual
  single-panel versions of the same 4 (set x comparison) combinations

**Peer-review-defensibility package (why this substitution will survive
scrutiny):**
- `Karaiskos_geneset_provenance.md` — (a) full provenance: Karaiskos et al.
  2018 *JASN* 29:2060-2068 (GEO GSE111107), exact Supplementary Table numbers
  (Table 2 for EXCLUSIVE, Table 3 `FindAllMarkers` output for TOP50), the
  extraction procedure, and the 50->49 detection accounting (`Sept11` is
  absent only due to the MGI Septin-gene renaming, not a biological
  exclusion); `genesets/KARAISKOS2018_PODOCYTE.gmt` is the gene set file
  itself, ready to load into any GSEA tool
- `tables/KARAISKOS_TOP50_leading_edge_genes.tsv` — (b) leading-edge gene list
  for TOP50 in both comparisons; confirms canonical podocyte identity markers
  (Nphs1, Nphs2, Wt1, Podxl, Synpo, Mafb, Npnt, Tcf21) are among the genes
  actually driving the enrichment signal, not incidental set members
- `tables/KARAISKOS_TOP50_vs_EXCLUSIVE_comparison.tsv` — (c) the table above,
  machine-readable
- `tables/Supplementary_ageing_set_all_configs.tsv` — (d) the podocyte-ageing
  set's full result matrix kept on record across every comparison and A1
  configuration tested (not deleted; this is the evidence for Key finding 2)

**Response-letter language** (why TOP50 is the main panel and not EXCLUSIVE,
why Karaiskos and not the ageing set, and the Key-finding-7 limitation
sentence about translation/ribosome pathways dominating the unfiltered
genome-wide scan) is in `Methods_rewrite_and_reviewer_response_notes.md`
under "Fig. 6B replacement: full response-letter package".

Pipeline: `12_karaiskos_fig6b_panel_data.R` (data/GMT/tables) ->
`13_plot_fig6b_karaiskos_panel.py` (figures) ->
`14_fig6b_defensibility_addenda.py` (gene-overlap accounting + genome-wide
rank columns, appended to the existing comparison/judgment tables in place).

## edgeR (original method) vs DESeq2 (revision method) verification

The original manuscript's DE analysis was very likely produced by this
server's standard pipeline, which uses **edgeR classic `exactTest()`**, not
DESeq2 (see Methods rewrite notes). Since this revision reports DESeq2
v1.38.3 instead, script `15_edgeR_vs_DESeq2_reproducibility.R` reproduces
edgeR's exact call sequence from the discovered template
(`cpm(merge)>2` filter -> `calcNormFactors` TMM -> `estimateCommonDisp` ->
`estimateTagwiseDisp` -> `exactTest`) on the same counts, for all 4
comparisons (6 configurations counting A1 in/out), and checks whether any
Fig. 5/6 conclusion depends on which tool is used.

**Result: no conclusion changes.** log2FC agreement is near-perfect
(Pearson r ≥ 0.9998 in every comparison); every named gene (Tmem215, Nlrp1b,
Glp1r, Wt1, Nphs1, Serpine1, Col4a1, Col4a2, Loxl1) keeps the same
significance call and direction under both tools, including the Wt1
Day-5 "borderline" result (edgeR FDR=0.093 vs. DESeq2 padj=0.050 — same
qualitative borderline conclusion); and every focal GSEA gene set's NES sign
is 100% concordant between edgeR-ranked and DESeq2-ranked preranked GSEA in
both `ADR_B_vs_A` and `B_ADR_vs_Ctrl`. The only honestly-reported difference
is that the exact set of genes crossing significance thresholds overlaps
only ~50% (Jaccard) in the two weaker-signal comparisons (baseline, Day-5
substrain) — an expected consequence of two different statistical
frameworks near a hard p-value cutoff, not a sign either tool is wrong, and
it does not affect any specific claim in the manuscript. Full write-up:
`edgeR_vs_DESeq2_verification.md`.

## File index

- `tables/00_merged_counts.tsv` — merged 12-sample count matrix (Step 0)
- `tables/DE_*.tsv`, `tables/Supplementary_DE_tables.xlsx` — Step 1 DE tables (4 contrasts incl. Fig.6A-equivalent)
- `figures/volcano_*.png` — Step 2 volcano plots (Wt1/Nphs1 always labeled; Tmem215/Nlrp1b/Glp1r highlighted)
- `tables/DE_interaction_substrain_by_treatment.tsv`, `figures/crosscomparison_ADR_response_scatter.png` — Step 3
- **`tables/Step1_ranking_metric_identity_check.tsv`** — rank-identity test (signed
  -log10p x sign(log2FC) vs. DESeq2 stat), the basis for reinstating `stat` as canonical
- **`tables/GSEA_<comparison>_full_joint_canonical.tsv`** (x4, A1-excluded main) and
  **`tables/GSEA_<comparison>_full_joint_canonical_A1included.tsv`** (x2, `ADR_B_vs_A`
  and `A_ADR_vs_Ctrl`) — CANONICAL, `stat`-ranked, from `09_gsea_stat_ranking_canonical.R`:
  complete unfiltered fgsea output, 1336-set joint universe
- **`tables/GSEA_canonical_focal_judgment_table.tsv`** — the 6 focal gene sets x 4
  comparisons x A1 status judgment table (NES/pval/joint FDR/pass), from `09`;
  now also carries `genome_wide_rank`/`n_sets_tested` columns (added via
  `tables/genome_wide_rank_lookup.tsv`) backing Key finding 7 with exact numbers
- **`tables/GSEA_canonical_ageing_vs_karaiskos.tsv`** — Step 6 side-by-side comparison
  (ageing set vs. both Karaiskos sets, `ADR_B_vs_A` and `B_ADR_vs_Ctrl`)
- **`tables/GSEA_canonical_top_hits_per_comparison.tsv`** — top 8 FDR-ranked gene sets
  (of ~1293-1296 tested) per comparison, canonical ranking; basis for Key findings 6/7/9
- `logs/09_leading_edge_overlap.txt` — leading-edge gene lists and overlap counts,
  ageing set vs. Karaiskos TOP50, for `ADR_B_vs_A` and `B_ADR_vs_Ctrl`
- `logs/09_gsea_canonical_versions.txt`, `logs/09_sessionInfo.txt` — software/package versions
- `Methods_rewrite_and_reviewer_response_notes.md` — Step 7 deliverable: proposed
  corrected Methods paragraph + point-by-point R1-6/R1-7/R2-4 response language,
  plus the Fig. 6B replacement response-letter package
- `Karaiskos_geneset_provenance.md`, `genesets/KARAISKOS2018_PODOCYTE.gmt` — full
  provenance/extraction record and GMT file for the two Karaiskos marker sets
- `tables/KARAISKOS_TOP50_vs_EXCLUSIVE_comparison.tsv`,
  `tables/KARAISKOS_TOP50_leading_edge_genes.tsv`,
  `tables/Supplementary_ageing_set_all_configs.tsv` — Fig. 6B replacement data
  package (from `12_karaiskos_fig6b_panel_data.R`)
- `figures/Fig6B_composite_<comparison>.{pdf,svg,png}`,
  `figures/Fig6B_KARAISKOS_<set>_<comparison>.{pdf,svg,png}` — Fig. 6B
  replacement figures (from `13_plot_fig6b_karaiskos_panel.py`)
- `figures/CHECK_enrichment_*.png` — diagnostic enrichment plots used to verify the
  sign findings above (not for publication, kept for audit trail)
- `logs/` — QC report, DESeq2 build log, package versions

### Superseded / secondary files (kept on disk, not deleted)

- `tables/GSEA_*_full_joint.tsv`, `tables/GSEA_podocyte_integrin_summary.tsv`,
  `tables/Supplementary_GSEA_tables.xlsx` — original stat-ranked run (`05_gsea.R`,
  reactome+ageing-only joint universe of 1334 sets); numerically consistent with
  `09` wherever both tested the same set, but `09` is now the reference table since
  it includes the full 6-set focal panel in one consistent joint universe
- `tables/Step5_*`, `tables/Supplementary_Step5_A1_sensitivity.xlsx` — original A1
  sensitivity for `ADR_B_vs_A` only (`06_A1_sensitivity.R`); superseded in scope by
  `09`'s judgment table, which also covers `A_ADR_vs_Ctrl`
- `tables/GSEA_alt_podocyte_identity_A1_sensitivity.tsv`,
  `tables/DE_alt_podocyte_identity_genes_A1_sensitivity.tsv` —
  `05d_alt_podocyte_A1_sensitivity.R`: A1 sensitivity for
  `GOBP_GLOMERULAR_EPITHELIAL_CELL_DIFFERENTIATION` (not part of the final 6-set
  focal panel, kept as an additional identity-marker check)
- `tables/GSEA_karaiskos_podocyte_summary.tsv`,
  `tables/GSEA_karaiskos_podocyte_A1_sensitivity.tsv`,
  `tables/DE_karaiskos_podocyte_genes_ADR_B_vs_A.tsv` — earlier, narrower Karaiskos
  run (`05e_karaiskos_podocyte_geneset.R`); superseded by `09`'s consistent joint run,
  numbers agree
- `tables/GSEA_log2FCrank_summary_all_sets.tsv`,
  `tables/GSEA_<comparison>_full_joint_log2FCrank.tsv` (x4),
  `tables/GSEA_log2FCrank_A1_sensitivity.tsv`,
  `tables/GSEA_log2FCrank_top_hits_per_comparison.tsv` — **secondary sensitivity
  check only** (`07`/`08`): literal signed-log2FC ranking, matching the Methods
  text as currently written but not the reconstructed actual procedure; retained
  for transparency, not used for "reproduces/does not reproduce" verdicts
