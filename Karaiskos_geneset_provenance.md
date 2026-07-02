# Karaiskos et al. 2018 podocyte marker gene sets: provenance and extraction record

Peer-review defensibility package for the proposed Fig. 6B replacement
(`KARAISKOS2018_PODOCYTE_TOP50`, main; `KARAISKOS2018_PODOCYTE_EXCLUSIVE`,
confirmatory). This document exists so the gene-set construction can be
audited independently of this project's code.

## Source

Karaiskos N, Rajewsky N, Kocks C, et al. "A Single-Cell Transcriptome Atlas of
the Mouse Glomerulus." *J Am Soc Nephrol* 2018;29(8):2060-2068.
GEO accession: **GSE111107**. Droplet-based (Drop-seq) single-cell RNA-seq of
~12,000 cells dissociated from real, intact mouse glomeruli (not a cultured
cell line, not a computationally deconvolved bulk signature). Podocytes were
identified as one of the annotated clusters via unsupervised clustering plus
canonical marker expression (the paper's own analysis, not ours).

Full text was obtained from the article's own supplementary data
(`ASN.2018030238SupplementaryData1.pdf`, via the journal's direct supplemental
link, since the PMC mirror served a bot-check page instead of the PDF) and
parsed with `pdftotext -layout`.

## Extraction procedure (exactly what was done, so it can be repeated)

1. **`KARAISKOS2018_PODOCYTE_EXCLUSIVE` (n=12).** Taken verbatim from the
   paper's **Supplementary Table 2**, which lists the authors' own most
   stringent "cell-type-exclusive" marker genes per cluster (genes detected
   essentially only in that one cluster, by the paper's own criterion — we did
   not re-derive or re-threshold this list). The 12 genes listed under
   Podocytes in that table are used unmodified: `Nphs2, Cdkn1c, Tcf21, Enpep,
   Nphs1, Synpo, Npnt, Wt1, Pard3b, Ptpro, Iqgap2, Mafb`.
2. **`KARAISKOS2018_PODOCYTE_TOP50` (n=49 of 50).** Taken from the paper's
   **Supplementary Table 3** ("full dataset"), the `FindAllMarkers`
   (Seurat-style differential-expression-vs-other-clusters) output for the
   Podocyte cluster, ranked by the paper's own `avg_logFC`. The top 50 genes
   by this ranking were extracted (regex-parsed from the 4th name/value pair
   per line, matching the paper's column order Endothelium / Immune /
   Mesangium / Podocytes / Tubules). **We did not look at our own expression
   data to choose or re-rank these genes** — the list is fixed by the
   original paper's own marker ranking before it was ever intersected with our
   DESeq2 output.
3. **Mapping to our count matrix.** Of the 50 extracted symbols, **49/50
   (98%) were found** in our DESeq2 gene universe; the fgsea `size` field
   for this set is consistently 49 across all comparisons tested. The one
   symbol not found is **`Sept11`** — this is a nomenclature artifact, not a
   biological exclusion: MGI/HGNC renamed the Septin gene family
   (`Sept1`-`Sept14` → `Septin1`-`Septin14`) specifically to avoid
   spreadsheet-autocorrect corruption of gene symbols, so the paper's original
   symbol `Sept11` does not match our reference annotation's current symbol.
   No other exclusion, filtering, or curation was applied.

## EXCLUSIVE is (almost) a subset of TOP50, not an independent replication

`KARAISKOS2018_PODOCYTE_EXCLUSIVE` (n=12) and `KARAISKOS2018_PODOCYTE_TOP50`
(n=50 raw / 49 detected) come from different supplementary tables of the same
paper and the same underlying cluster, so their overlap must be reported
honestly rather than implied to be independent: **11 of the 12 EXCLUSIVE
genes (91.7%) are also members of TOP50** (`Nphs2, Cdkn1c, Enpep, Tcf21,
Nphs1, Synpo, Npnt, Wt1, Pard3b, Iqgap2, Mafb`); only **`Ptpro`** is unique to
EXCLUSIVE and absent from TOP50. TOP50 in turn contains 38 genes not present
in EXCLUSIVE.

**Correct framing:** the two panels are **not two independent marker sets**
in the statistical sense — EXCLUSIVE is essentially a curated 12-gene core
that is nearly nested inside the larger 49-gene TOP50 list. What the
side-by-side panels actually demonstrate is that **the enrichment signal is
not being driven by, or dependent on, the 38 TOP50-only genes** — the small,
almost-fully-overlapping "hard core" reproduces the same direction and
significance on its own. This is a meaningfully weaker (but still useful)
form of robustness evidence than two genuinely independent gene lists would
provide, and should be described in those terms in the manuscript/response
letter rather than as "two independent gene sets agree." See
`tables/KARAISKOS_TOP50_vs_EXCLUSIVE_comparison.tsv` (columns
`n_overlap_with_other_set`, `pct_EXCLUSIVE_contained_in_TOP50`) for the
machine-readable record.

## Files

- `outputs/genesets/KARAISKOS2018_PODOCYTE.gmt` — both gene sets in standard
  GMT format (name, description/source, then gene list), ready to load into
  any GSEA tool.
- `outputs/tables/KARAISKOS_TOP50_vs_EXCLUSIVE_comparison.tsv` — NES/nominal
  P/joint FDR for both sets, both comparisons (`ADR_B_vs_A`, `B_ADR_vs_Ctrl`).
- `outputs/tables/KARAISKOS_TOP50_leading_edge_genes.tsv` — leading-edge gene
  list for TOP50 in both comparisons, flagged for canonical podocyte identity
  markers (Nphs1, Nphs2, Wt1, Podxl, Synpo, Ptpro, Tcf21, Mafb, Cdkn1c, Ctsl,
  Mertk, Npnt).

## Why this is not cherry-picking

The gene list membership and rank order were fixed by the original authors'
own single-cell differential-expression analysis, published independently of
and prior to this re-analysis. We selected which *published* marker list to
test (a decision made and documented in `README_analysis_log.md`, Key finding
4) but did not select, reorder, or filter *which genes* go into that list
based on how they behaved in our own ADR/substrain data. This is the same
standard used for the manuscript's own choice of
`TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING` from MSigDB.
