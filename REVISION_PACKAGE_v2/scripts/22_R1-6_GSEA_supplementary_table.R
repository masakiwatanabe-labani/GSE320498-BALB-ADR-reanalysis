#!/usr/bin/env Rscript
# R1-6: consolidate the canonical (stat-ranked) full joint-universe GSEA
# results for all 4 primary comparisons into a single supplementary
# workbook, one sheet per comparison, each captioned with software/package
# versions, ranking metric, gene-set universe size, and the significance
# threshold (FDR/padj < 0.05, matching the manuscript's own stated
# criterion -- NOT nominal p, which the original Fig. 6B/6D text used).
#
# Source data: tables/GSEA_<comparison>_full_joint_canonical.tsv (x4), from
# 09_gsea_stat_ranking_canonical.R. No new GSEA run -- this only repackages
# the already-canonical output into one reviewer-facing file.

suppressMessages(library(openxlsx))

indir <- "/usr/local/jupyter/ADR_BALB"
outdir_main <- file.path(indir, "outputs")
outdir_v2 <- file.path(outdir_main, "REVISION_PACKAGE/REVISION_PACKAGE_v2")

versions <- readLines(file.path(outdir_main, "logs/09_gsea_canonical_versions.txt"))
cat(paste(versions, collapse = "\n"), "\n\n")

comparisons <- list(
  baseline_B_vs_A = "Baseline (untreated): BALB/cByJcl (B) vs BALB/cAJcl (A). Positive NES = enriched (higher) in ByJcl. Samples: A-Ctrl n=3 vs B-Ctrl n=3.",
  ADR_B_vs_A      = "Day-5 post-ADR: BALB/cByJcl (B) vs BALB/cAJcl (A). Positive NES = enriched (higher) in ByJcl at Day 5. Samples: A-ADR n=2 (A-ADR1 excluded, low glomerular purity; see main-analysis sensitivity check) vs B-ADR n=3. This is the comparison underlying the manuscript's Fig. 6B/6D.",
  A_ADR_vs_Ctrl   = "ADR response within BALB/cAJcl (A): ADR-treated vs Ctrl. Positive NES = enriched (higher) after ADR. Samples: A-Ctrl n=3 vs A-ADR n=2 (A-ADR1 excluded).",
  B_ADR_vs_Ctrl   = "ADR response within BALB/cByJcl (B): ADR-treated vs Ctrl. Positive NES = enriched (higher) after ADR. Samples: B-Ctrl n=3 vs B-ADR n=3."
)

method_caption <- paste0(
  "Method: preranked GSEA, fgsea R package v1.24.0 (multilevel algorithm, eps=0, minSize=5, maxSize=500), ",
  "gene sets from msigdbr v26.1.0 (MSigDB db_version 2026.1.Mm, mouse-native symbols): Reactome M2:CP:REACTOME ",
  "(1,333 sets) plus 3 custom podocyte-focused sets (TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING, ",
  "KARAISKOS2018_PODOCYTE_EXCLUSIVE, KARAISKOS2018_PODOCYTE_TOP50) = 1,336 gene sets tested jointly per ",
  "comparison. Ranking metric: DESeq2 v1.38.3 Wald test statistic (signed, significance-weighted; see Methods ",
  "rewrite notes for rationale). FDR: Benjamini-Hochberg, computed jointly across the FULL tested gene-set ",
  "universe for that comparison (columns: pval = nominal P: padj = BH-adjusted FDR/q-value). Significance ",
  "threshold: padj (FDR) < 0.05, per the manuscript's own stated criterion -- NOT the nominal p-value. ",
  "Random seed 20260220, reset immediately before every fgsea() call."
)

wb <- createWorkbook()
summary_rows <- list()

for (nm in names(comparisons)) {
  df <- read.delim(file.path(outdir_main, sprintf("tables/GSEA_%s_full_joint_canonical.tsv", nm)),
                    stringsAsFactors = FALSE)
  df <- df[order(df$pval), ]
  df$significant_FDR0.05 <- ifelse(!is.na(df$padj) & df$padj < 0.05, "YES", "no")
  df <- df[, c("pathway", "size", "ES", "NES", "pval", "padj", "significant_FDR0.05", "log2err", "leadingEdge")]
  colnames(df)[colnames(df) == "pval"] <- "nominal_pval"
  colnames(df)[colnames(df) == "padj"] <- "FDR_qvalue_joint"

  caption <- paste0(comparisons[[nm]], " | ", method_caption,
                     " | n gene sets tested = ", nrow(df), ".")

  addWorksheet(wb, nm)
  writeData(wb, nm, caption, startRow = 1, startCol = 1)
  mergeCells(wb, nm, cols = 1:9, rows = 1)
  addStyle(wb, nm, createStyle(textDecoration = "italic", wrapText = TRUE, fontSize = 9), rows = 1, cols = 1:9)
  setRowHeights(wb, nm, rows = 1, heights = 90)
  writeDataTable(wb, nm, df, startRow = 3, tableStyle = "TableStyleLight9")
  setColWidths(wb, nm, cols = 1:9, widths = c(44, 8, 10, 10, 12, 12, 14, 10, 60))
  freezePane(wb, nm, firstActiveRow = 4)

  n_sig <- sum(df$FDR_qvalue_joint < 0.05, na.rm = TRUE)
  summary_rows[[nm]] <- data.frame(comparison = nm, n_gene_sets_tested = nrow(df),
                                    n_significant_FDR0.05 = n_sig)
  cat(sprintf("%s: %d gene sets tested, %d significant at FDR<0.05\n", nm, nrow(df), n_sig))
}

overview <- do.call(rbind, summary_rows)
addWorksheet(wb, "overview")
writeData(wb, "overview",
          "Overview: gene sets tested and significant (FDR<0.05) per comparison. See README_analysis_log.md and Methods_rewrite_and_reviewer_response_notes.md (parent REVISION_PACKAGE folder) for full context.",
          startRow = 1, startCol = 1)
mergeCells(wb, "overview", cols = 1:3, rows = 1)
addStyle(wb, "overview", createStyle(textDecoration = "italic", wrapText = TRUE), rows = 1, cols = 1:3)
setRowHeights(wb, "overview", rows = 1, heights = 45)
writeDataTable(wb, "overview", overview, startRow = 3, tableStyle = "TableStyleLight9")
worksheetOrder(wb) <- c(which(names(wb) == "overview"), which(names(wb) != "overview"))

out_fp <- file.path(outdir_v2, "tables/TableS_GSEA_all_comparisons_canonical.xlsx")
saveWorkbook(wb, out_fp, overwrite = TRUE)
cat("\nSaved:", out_fp, "\n")
cat("\n=== 22_R1-6_GSEA_supplementary_table.R complete ===\n")
