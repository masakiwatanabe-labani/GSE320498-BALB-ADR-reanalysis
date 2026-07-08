# Revision package manifest

Everything needed for the reviewer response, consolidated into this one
folder. Start with `README_analysis_log.md` — it is the full analysis log and
the index of everything else here.

## Read first

- `README_analysis_log.md` — full analysis log: methodological history,
  correspondence table (review comment -> output), 9 key findings, file index.
- `Methods_rewrite_and_reviewer_response_notes.md` — proposed corrected
  Methods paragraph + point-by-point R1-6/R1-7/R2-4 response language +
  Fig. 6B replacement response-letter package.
- `Karaiskos_geneset_provenance.md` — full sourcing/extraction record for the
  Karaiskos et al. 2018 podocyte marker gene sets, including the
  TOP50/EXCLUSIVE overlap accounting.
- `TODO_not_recoverable_from_counts.md` — items reviewers asked for that
  cannot be derived from the deposited count matrices alone, including the
  upstream-pipeline version search (this shared server's tool versions) and
  the confirmed public GEO record check.
- `edgeR_vs_DESeq2_verification.md` — verifies that switching this revision's
  DE engine from the platform's original method (edgeR classic `exactTest()`)
  to DESeq2 v1.38.3 changes no conclusion behind Fig. 5/6.
- `R1-5_response_paragraphs.md` — one response-letter paragraph per R1-5 item
  (contrast sign/labeling, baseline gene description, supplementary DE tables,
  ADR-response concordance), each citing the exact numbers in this package.
- `figures/FigS_A1_sensitivity_NES_comparison.pdf` — R1-7/R2-4 A1 (low-purity
  A-ADR1) sensitivity check: NES with A1 excluded vs. included, all 6 focal
  gene sets, both A1-applicable comparisons; only podocyte-ageing shows a
  sign flip that is significant on both sides (highlighted in vermillion).
- `figures/FigS_A1_contamination_QC.pdf` — QC basis for excluding A-ADR1:
  per-marker renal-tubular-epithelial CPM (A-ADR1 highest on 5/6 markers) and
  podocyte:tubular marker-sum ratio (A-ADR1 = 2.51, lowest of all 6 ADR
  samples, driven by excess tubular signal, not depleted podocyte signal).
- `ZENODO_MANIFEST.txt` — sha256 checksums + sizes for every file in this
  package, for Zenodo deposit integrity verification.

## Folders

- `tables/` — every DE table, GSEA result table (canonical stat-ranked
  primary + log2FC-ranked secondary sensitivity check), judgment/comparison
  tables, and the `enrichment_curves/` subfolder (raw curve/hit/ranking-vector
  data behind every Fig.6-style plot). Large binary DESeq2 objects (`.rds`)
  are intentionally excluded from this package — they are R-only intermediate
  objects, not reviewable data; re-generate them from `scripts/01_build_dds.R`
  if needed to re-run later pipeline steps.
- `figures/` — all figures, including the Fig.6C/6D-style single-panel plots
  (`Fig6_enrichment_*`) and the Fig.6B replacement panels
  (`Fig6B_*`, including the two composite main+confirmatory figures), each in
  PDF/SVG/PNG.
- `genesets/` — `KARAISKOS2018_PODOCYTE.gmt`, standard-format gene sets ready
  to load into any GSEA tool.
- `scripts/` — the full numbered pipeline (`00`-`18`), run in order, that
  reproduces every table and figure in this package, including the
  edgeR-vs-DESeq2 verification (`15`), the final R1-5 deliverables (`16`,
  repackages/relabels existing DE output only, no new model fit), the
  R1-7/R2-4 A1 sensitivity figure (`17`, repackages the existing canonical
  GSEA judgment table, no new GSEA run), and the A-ADR1 contamination QC
  figure (`18`, recomputes CPM directly from `tables/00_merged_counts.tsv`,
  no new sequencing analysis).
- `logs/` — package versions, sessionInfo, and QC/build logs.

## Most load-bearing individual files

- `tables/GSEA_canonical_focal_judgment_table.tsv` — the single table
  underlying nearly every "reproduces/does not reproduce" verdict (NES,
  nominal P, joint FDR, pass/fail, and genome-wide rank, for all 6 focal gene
  sets x 4 comparisons x A1 status).
- `tables/KARAISKOS_TOP50_vs_EXCLUSIVE_comparison.tsv` — Fig. 6B replacement
  statistics + gene-overlap accounting between the two marker panels.
- `figures/Fig6B_composite_ADR_B_vs_A.pdf`,
  `figures/Fig6B_composite_B_ADR_vs_Ctrl.pdf` — recommended replacement
  figures for Fig. 6B.
- `figures/Fig6_enrichment_INTEGRIN_SIGNALING.pdf` — recommended replacement
  figure for Fig. 6D.
- `figures/Fig5_volcano_baseline.pdf`, `figures/Fig6A_volcano_D5.pdf` —
  R1-5 replacement volcanoes: contrast sign stated in the title, Wt1/Nphs1
  always labeled, Tmem215/Nlrp1b/Glp1r highlighted with exact log2FC/FDR
  printed on the figure.
- `tables/TableS_top_baseline_genes.csv` — R1-5 baseline gene-description
  table (FDR<0.05 & |log2FC|>1, 64 genes), replacing the "modest" wording
  with an explicit count and named examples.
- `tables/TableS_DE_all_comparisons.xlsx` — R1-5 supplementary DE table, one
  sheet per primary contrast (baseline B vs A, Day-5 B vs A, A ADR vs Ctrl,
  B ADR vs Ctrl), each captioned with its sign convention.
- `figures/FigS_ADR_response_concordance.png` — R1-5 cross-substrain
  ADR-response concordance (Pearson r=0.54; 99.6% direction-concordant among
  genes significant in both substrains).
