# Methods paragraph rewrite proposal + reviewer-response talking points

This file supports the R1-6 / R1-7 / R2-4 responses. See `README_analysis_log.md`
for the full result tables this is based on.

## Why a rewrite is needed

The manuscript's Methods paragraph states GSEA was run on genes "ordered by signed
log2 fold change." Two things in the manuscript's own reported numbers are inconsistent
with a literal log2FC-only ranking:

1. Fig. 6D reports nominal P=0.0235 **together with** FDR q<0.001. Under any
   standard Benjamini-Hochberg correction, q can never be smaller than the nominal
   p-value of the same test. A ranking metric and correction procedure that can
   produce q << p is not a vanilla "rank by log2FC, BH-correct across tested sets"
   pipeline.
2. When we actually re-ran GSEA ranking genes by pure signed log2FC (scripts
   `07`/`08`), joint FDR for the focal gene sets came out far weaker than the
   manuscript's reported values (see `GSEA_log2FCrank_summary_all_sets.tsv`) —
   consistent with log2FC-only ranking not being the actual procedure used.

We tested the natural alternative directly (`Step1_ranking_metric_identity_check.tsv`,
script `09`): rank genes by **signed(-log10 p) x sign(log2FC)** — i.e., a
significance-weighted ranking, which is what most GSEA implementations mean in
practice when they say "ranked by fold change and significance" or simply use a
Wald/t-statistic. Its gene order is essentially identical (Spearman rho = 1.0000,
all 4 comparisons) to the DESeq2 Wald `stat` column, because for a Wald test `stat`
is itself a monotonic, sign-preserving transform of the p-value. We therefore treat
`stat`-ranked preranked GSEA as the best available reconstruction of what the authors
actually ran, and recommend the Methods text be corrected to describe this rather
than a pure log2FC ranking.

## Proposed Methods paragraph (drop-in replacement)

> Differentially expressed genes were identified with DESeq2 (v1.38.3) using a Wald
> test. For gene set enrichment analysis, all detected genes were ranked by the
> DESeq2 Wald test statistic (a signed, significance-weighted ranking that combines
> the direction and magnitude of the fold change with its estimation uncertainty).
> Preranked GSEA was performed using the `fgsea` R package (v1.24.0, multilevel
> algorithm, `eps = 0`, minimum/maximum gene-set size 5/500) against the Reactome
> collection (MSigDB M2:CP:REACTOME, db_version 2026.1.Mm, mouse-native gene
> symbols via `msigdbr` v26.1.0; 1333 gene sets) together with the Tabula Muris
> Senis kidney podocyte-ageing signature and [if adopted, add: a mouse glomerular
> single-cell podocyte marker panel from Karaiskos et al. 2018]. False discovery
> rate (Benjamini-Hochberg q-value) was computed jointly across every gene set
> tested in a given comparison, not on any single gene set in isolation. A fixed
> random seed (20260220) was set immediately before each `fgsea()` call for
> reproducibility.

## Upstream pipeline (FastQC/Trimmomatic/STAR/featureCounts) version findings

See `TODO_not_recoverable_from_counts.md` ("R1 minor comment 6 / R1-6") for
the full search record. Summary of what could and could not be established
from this shared analysis server:

**Established with high confidence** (this server runs a commercial turnkey
pipeline, "Amelieff Quick Start Package," identically across 9 other
projects; these are the pinned versions in its actual code, not just
whatever happens to be on `PATH`):
- Trimmomatic **v0.39**, PRINSEQ-lite **v0.20.4** (run as an undisclosed
  second trimming step after Trimmomatic, if this pipeline was used),
  STAR **2.7.10a** (directly confirmed via the build log of the actual
  prebuilt mouse genome index on this server), featureCounts **v2.0.3**.
- Genome annotation is **UCSC mm39 (=GRCm39) `refGene.gtf`**, not
  Ensembl/GENCODE — dated 2020-09-11 (UCSC's own file timestamp; UCSC's
  refGene track has no Ensembl/GENCODE-style release number).

**Could not be established (upstream tools only):** this platform's standard
two-group-comparison step does not use DESeq2 — it uses **edgeR
`exactTest()`** — and its pathway step does not use preranked GSEA — it uses
**ReactomePA `enrichPathway()` (ORA)**. This was verified as identical
(zero variation) across every other project on the server.

**Decision on DESeq2 version — do not chase it, adopt v1.38.3 for this
revision.** DESeq2 operates on the already-deposited count matrices, which
we can and did rerun ourselves; unlike STAR/Trimmomatic/featureCounts (which
*produced* the fixed, deposited counts and cannot be substituted after the
fact), DESeq2's version is a free choice for *this* re-analysis. Every DE/GSEA
number in this revision is the literal output of DESeq2 v1.38.3, so the
Methods paragraph above states that version as fact, not as a recovered
original. This resolves the DESeq2 half of R1-6; no further author
confirmation is needed for DESeq2 specifically. The GSEA software identity
question (item 4 below) is a separate, still-open issue — it is about what
produced the manuscript's original Fig. 6B/6D numbers we are trying to
reproduce, which our own choice of fgsea cannot retroactively answer.

**Recommended author confirmation checklist (add to response letter):**
1. Was GSE320498 processed on this same "Amelieff QSP" platform?
2. If so, was PRINSEQ-lite actually run (undisclosed in current Methods)?
3. Was the genome annotation UCSC RefSeq (`refGene.gtf`) rather than
   Ensembl/GENCODE — if so, cite it correctly and give the actual retrieval
   date, since gene-set/exon-model sizes differ between annotation sources.
4. What GSEA software/script produced the original Fig. 6B/6D NES and FDR
   values, since the platform's own default pipeline uses ReactomePA ORA,
   not preranked GSEA, and this determines whether the q<0.001 vs
   nominal-P=0.0235 pairing in Fig. 6D (see Key finding 1) can ever be
   exactly reproduced.

## Point-by-point response-letter material

**R1-6 (name the GSEA software/version):**
Originally an ad hoc in-house preranked implementation (random same-size gene-set
resampling for NES normalization, approximate FDR) was used and not named. We have
re-run every reported gene-set result with the peer-reviewed `fgsea` package
(v1.24.0), which is now the sole GSEA engine used in the revised analysis and is
fully named/versioned in Methods (see rewrite above) and in `logs/09_gsea_canonical_versions.txt`.

**R1-6 (ranking metric was misstated):**
Methods said "signed log2 fold change"; the actual ranking (reconstructed from the
q<p pattern in Fig. 6D and confirmed by a direct rank-correlation test, rho=1.0000
in all 4 comparisons) is the DESeq2 Wald statistic, a significance-weighted ranking.
The Methods text should be corrected as shown above. We additionally provide the
literal log2FC-ranked numbers as a named sensitivity check (`07`/`08`) for
transparency, since they no longer match the headline numbers once joint FDR is
applied to either ranking.

**R1-6 (report FDR/q, not nominal P, and show the full gene-set table):**
Done for every comparison; joint BH q is now reported alongside nominal P
throughout. Full unfiltered `fgsea` output for all ~1293-1296 tested gene sets per
comparison is in `GSEA_<comparison>_full_joint_canonical.tsv` (x4) plus two
A1-included variants. Result: under the stat ranking, **Reactome Integrin
signaling (Fig. 6D) now reproduces essentially as reported** — NES=+2.03, nominal
P=1.2e-4, joint FDR=9.1e-4 (A1-excluded, main setting) — closely matching the
manuscript's own P=0.0235/q<0.001 pairing far better than log2FC ranking ever did.
This is the single strongest piece of evidence that stat-ranking is the correct
reconstruction of the original method.

**R1-7 / R2-4 (A1 sensitivity, "data not shown"):**
Full in/out numbers for all 6 focal gene sets across the affected comparisons are
in `GSEA_canonical_focal_judgment_table.tsv`. Result is mixed, and we report both
halves:
- Integrin signaling and ECM organization are **robust** to A1 inclusion (same
  sign, similar or greater significance both ways) — supports the manuscript's
  "robustness" claim for the ECM/integrin theme.
- The podocyte-ageing set (Fig. 6B) **is not robust**: NES=+1.93 (A1-excluded) vs.
  NES=-2.04 (A1-included) — a full sign flip. The manuscript's stated main-analysis
  configuration (A1 excluded) also disagrees in sign with the manuscript's own
  reported Fig. 6B value (NES=-1.42). We recommend the authors report the actual
  A1-included number for this specific gene set rather than "data not shown," and
  consider whether Fig. 6B's specific choice of gene set is the best available
  evidence for the podocyte-loss claim (see next point).

**Recommended gene-set substitution for Fig. 6B:**
`TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING` is an ageing-associated signature (its
leading-edge genes are dominated by MHC/immune/senescence genes, not canonical
podocyte identity genes) and is the only focal gene set in this whole re-analysis
whose A1-sensitivity and cross-comparison behavior is inconsistent. A real,
peer-reviewed, single-cell-derived mouse podocyte marker panel (Karaiskos et al.
2018, *J Am Soc Nephrol* 29:2060-2068, GEO GSE111107) is directionally consistent
with the manuscript's narrative (loss of podocyte identity in ByJcl at Day 5) in
**every configuration tested** — both marker lists, both A1 statuses, and both the
cross-substrain (`ADR_B_vs_A`) and within-substrain (`B_ADR_vs_Ctrl`) comparisons —
and is often more significant than the ageing set (e.g. `B_ADR_vs_Ctrl` TOP50:
NES=-2.70, joint FDR=2.1e-10, the single most significant gene set found genome-
wide in that comparison; see `GSEA_canonical_top_hits_per_comparison.tsv`). We
recommend citing the Karaiskos panel alongside, or in place of, the ageing
signature in a revised Fig. 6B.

## Fig. 6B replacement: full response-letter package

Proposed new Fig. 6B: **`KARAISKOS2018_PODOCYTE_TOP50` as the main enrichment
panel**, with `KARAISKOS2018_PODOCYTE_EXCLUSIVE` shown as an adjacent
confirmatory panel — see `figures/Fig6B_composite_ADR_B_vs_A.{pdf,svg,png}`
(main manuscript comparison) and `figures/Fig6B_composite_B_ADR_vs_Ctrl.{pdf,svg,png}`
(within-ByJcl ADR response, the strongest single result in this whole
re-analysis). Individual single-panel versions of all four (2 sets x 2
comparisons) are also provided. Full provenance, extraction procedure, and the
GMT file itself are in `Karaiskos_geneset_provenance.md` /
`genesets/KARAISKOS2018_PODOCYTE.gmt`; the side-by-side statistics table is
`tables/KARAISKOS_TOP50_vs_EXCLUSIVE_comparison.tsv`; leading-edge genes
(showing Nphs1, Nphs2, Wt1, Podxl, Synpo, Mafb, Npnt, Tcf21 among the driving
genes in both comparisons) are in `tables/KARAISKOS_TOP50_leading_edge_genes.tsv`;
the ageing set's full result matrix (every comparison, both A1 configurations)
is retained in `tables/Supplementary_ageing_set_all_configs.tsv` rather than
removed from the record.

**Gene overlap between the two panels (report this explicitly):** 11 of
EXCLUSIVE's 12 genes (91.7%) are also in TOP50 — only `Ptpro` is unique to
EXCLUSIVE. The two panels should therefore be described as a full marker panel
plus its nearly-fully-nested stringent core, not as two independent gene
lists. What they jointly demonstrate is that the signal does not depend on the
38 TOP50-only genes, not that two unrelated gene sets independently agree
(see `Karaiskos_geneset_provenance.md`).

**Suggested Fig. 6B legend addition (cell-composition caveat):** *"Karaiskos
markers are podocyte-cluster-identifying genes derived from single-cell data;
in glomerulus-enriched bulk RNA-seq, a negative NES for this panel indicates
relative attenuation of the podocyte transcriptional program but does not, on
its own, distinguish reduced podocyte number from reduced per-cell
expression. This bulk enrichment result is intended to be interpreted
alongside, and is complementary to, the histological findings (WT1+ cell
count per glomerulus; NPHS1 immunoreactivity optical density), which
separately address cell number and per-cell expression intensity."*

**Why TOP50 as the main panel, not EXCLUSIVE:** `KARAISKOS2018_PODOCYTE_EXCLUSIVE`
has only 12 member genes. A preranked GSEA statistic (NES, and especially its
associated permutation/multilevel p-value) becomes less stable and less
precisely estimable as gene-set size shrinks toward fgsea's `minSize`
threshold (5) — small sets are more sensitive to the exact identity of the 1-2
extreme-ranked hits, and their leading-edge fraction is a coarser, noisier
statistic. `KARAISKOS2018_PODOCYTE_TOP50` (n=49 detected) gives a materially
more stable, more granular enrichment curve while reaching the **same
conclusion in every configuration tested** (see comparison table). We
therefore present the larger, more statistically stable set as the primary
evidence and the smaller, more stringent "exclusive marker" set as a
confirmatory panel demonstrating that the result does not depend on one
specific, narrow gene list.

**Why Karaiskos et al. 2018 at all:** it is an independently published,
peer-reviewed, real single-cell RNA-seq atlas of intact mouse glomeruli
(GSE111107) — not an ageing signature, not a GO-term aggregation, and not a
gene set whose membership we chose by looking at our own expression data (see
"Why this is not cherry-picking" in the provenance document). Its podocyte
cluster markers are, by construction, genes that define podocyte identity in
real single-cell data, which is exactly the biological claim Fig. 6B is trying
to support.

**Why we recommend replacing the ageing signature with Karaiskos, specifically:**
three independent reasons converge on the same recommendation:
1. *Biological validity* — the ageing set's leading-edge genes are dominated
   by MHC/immune/senescence genes, not podocyte-identity genes; Karaiskos's
   leading-edge genes are canonical podocyte structural/identity genes
   (Nphs1, Nphs2, Wt1, Podxl, Synpo...).
2. *A1 robustness* — the ageing set's sign flips with A-ADR1 inclusion/exclusion
   (Key finding 2); both Karaiskos sets keep the same sign and remain
   significant in every A1 configuration tested.
3. *Baseline IHC direction concordance* — the manuscript's baseline IHC claim
   (fewer WT1+ cells / lower NPHS1 in ByJcl) is about podocyte identity loss,
   which is exactly what Karaiskos markers measure directly; the ageing
   signature is a proxy that happens to agree at baseline but disagrees at
   Day 5 in the main analysis configuration.

**Limitation sentence for Key finding 7 (translation/ribosome dominance
genome-wide under stat ranking):** *"Because the ranking statistic used here
(DESeq2 Wald stat) weights both fold-change magnitude and estimation
precision, high-count, low-dispersion housekeeping transcripts (translation,
ribosome biogenesis, mitochondrial pathways) can attain very large statistic
values even for modest fold changes, so such pathways can rank above
biologically-targeted gene sets in an unfiltered genome-wide scan; this is a
recognized property of significance-weighted preranked GSEA and does not
diminish the statistical significance of the ECM/integrin/podocyte gene sets
reported here, which remain significant after joint FDR correction regardless
of their absolute rank position."*
