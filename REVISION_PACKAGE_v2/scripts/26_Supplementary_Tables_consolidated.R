#!/usr/bin/env Rscript
# Consolidate every supplementary table established as necessary across this
# reviewer-response re-analysis (R1-5, R1-6, R1-7/R2-4, R2-2, plus the
# manuscript_values.md extraction task) into ONE submission-ready workbook,
# with a table-of-contents sheet and a caption row on every data sheet.
#
# No new statistical model: every sheet is read directly from existing,
# already-verified source files and written out with a caption. This script
# only repackages/renumbers.

suppressMessages(library(openxlsx))

indir <- "/usr/local/jupyter/ADR_BALB"
outdir_main <- file.path(indir, "outputs")
outdir_v2 <- file.path(outdir_main, "REVISION_PACKAGE/REVISION_PACKAGE_v2")

wb <- createWorkbook()

add_sheet <- function(sheet_name, caption, df, ncol_merge = NULL) {
  addWorksheet(wb, sheet_name)
  if (is.null(ncol_merge)) ncol_merge <- ncol(df)
  writeData(wb, sheet_name, caption, startRow = 1, startCol = 1)
  mergeCells(wb, sheet_name, cols = 1:ncol_merge, rows = 1)
  addStyle(wb, sheet_name, createStyle(textDecoration = "italic", wrapText = TRUE, fontSize = 9), rows = 1, cols = 1:ncol_merge)
  setRowHeights(wb, sheet_name, rows = 1, heights = 60)
  writeDataTable(wb, sheet_name, df, startRow = 3, tableStyle = "TableStyleLight9")
  freezePane(wb, sheet_name, firstActiveRow = 4)
  cat(sprintf("  %-24s %5d rows x %2d cols\n", sheet_name, nrow(df), ncol(df)))
}

sign_note <- "Positive log2FC = higher in BALB/cByJcl (B) [for baseline_B_vs_A and ADR_B_vs_A] or higher after ADR vs Ctrl [for A_ADR_vs_Ctrl / B_ADR_vs_Ctrl]. DESeq2 v1.38.3 Wald test."

# ============================================================================
# Table of contents (written last, after we know final row counts, but
# defined first for readability; placeholder rows filled in below)
# ============================================================================
toc_rows <- list()
add_toc <- function(id, title, source_script, n) {
  toc_rows[[length(toc_rows) + 1]] <<- data.frame(Table = id, Title = title, Source_script = source_script, N_rows = n)
}

cat("Building consolidated supplementary tables workbook...\n\n")

# ============================================================================
# Table S1: Baseline top DE genes (padj<0.05 & |log2FC|>1)
# ============================================================================
s1 <- read.csv(file.path(outdir_main, "tables/TableS_top_baseline_genes.csv"), stringsAsFactors = FALSE)
add_sheet("TableS1_baseline_top_genes", paste0(
  "Table S1. Differentially expressed genes at baseline (untreated), BALB/cByJcl vs BALB/cAJcl, ",
  "padj<0.05 AND |log2FC|>1 (64 of 19,662 tested genes). ", sign_note
), s1)
add_toc("S1", "Baseline top DE genes (padj<0.05 & |log2FC|>1)", "16_R1-5_final_deliverables.R", nrow(s1))

# ============================================================================
# Tables S2-S5: full DE results, 4 primary comparisons
# ============================================================================
de_specs <- list(
  list(id = "S2", sheet = "TableS2_DE_baseline_B_vs_A", file = "DE_baseline_B_vs_A.tsv",
       cap = "Table S2. Full differential expression results, baseline (untreated): BALB/cByJcl vs BALB/cAJcl (all 19,662 tested genes). Underlies Figure 5."),
  list(id = "S3", sheet = "TableS3_DE_ADR_B_vs_A", file = "DE_ADR_B_vs_A.tsv",
       cap = "Table S3. Full differential expression results, Day 5 post-adriamycin: BALB/cByJcl vs BALB/cAJcl, A-ADR1 excluded (main analysis; all 19,662 tested genes). Underlies Figure 6A."),
  list(id = "S4", sheet = "TableS4_DE_A_ADR_vs_Ctrl", file = "DE_A_ADR_vs_Ctrl.tsv",
       cap = "Table S4. Full differential expression results, ADR response within BALB/cAJcl (ADR-treated vs Ctrl, A-ADR1 excluded; all 19,662 tested genes)."),
  list(id = "S5", sheet = "TableS5_DE_B_ADR_vs_Ctrl", file = "DE_B_ADR_vs_Ctrl.tsv",
       cap = "Table S5. Full differential expression results, ADR response within BALB/cByJcl (ADR-treated vs Ctrl; all 19,662 tested genes).")
)
for (spec in de_specs) {
  df <- read.delim(file.path(outdir_main, "tables", spec$file), stringsAsFactors = FALSE)
  add_sheet(spec$sheet, paste0(spec$cap, " ", sign_note), df)
  add_toc(spec$id, spec$cap, "02_DE_tables.R", nrow(df))
}

# ============================================================================
# Table S6: full canonical GSEA results, all 4 comparisons (+ overview)
# ============================================================================
gsea_wb_path <- file.path(outdir_v2, "tables/TableS_GSEA_all_comparisons_canonical.xlsx")
gsea_sheet_map <- list(
  overview = "TableS6_GSEA_overview",
  baseline_B_vs_A = "TableS6_GSEA_baseline",
  ADR_B_vs_A = "TableS6_GSEA_D5_BvsA",
  A_ADR_vs_Ctrl = "TableS6_GSEA_A_ADRvsCtrl",
  B_ADR_vs_Ctrl = "TableS6_GSEA_B_ADRvsCtrl"
)
gsea_caption <- paste0(
  "Table S6. Canonical preranked GSEA (fgsea v1.24.0, DESeq2 Wald-statistic ranking) results, ALL 1,293 tested ",
  "Reactome + podocyte-focused gene sets (of 1,336-set joint universe), for all 4 primary comparisons. Joint ",
  "BH FDR computed across the full tested universe per comparison, not a focal subset. nominal_pval = nominal ",
  "p-value; FDR_qvalue_joint = BH-adjusted q-value (the manuscript's stated significance criterion is FDR<0.05, not nominal p)."
)
for (src_name in names(gsea_sheet_map)) {
  df <- read.xlsx(gsea_wb_path, sheet = src_name, startRow = 3)
  add_sheet(gsea_sheet_map[[src_name]], paste0(gsea_caption, " [", src_name, "]"), df)
}
add_toc("S6", "Canonical GSEA, all gene sets, all 4 comparisons (5 sub-sheets: overview + 4 comparisons)",
        "09_gsea_stat_ranking_canonical.R / 22_R1-6_GSEA_supplementary_table.R", NA)

# ============================================================================
# Table S7: Fig. 6C ORA full results (DESeq2-derived DEGs)
# ============================================================================
s7 <- read.csv(file.path(outdir_v2, "tables/TableS_ORA_D5_full.csv"), stringsAsFactors = FALSE)
add_sheet("TableS7_ORA_Fig6C_full", paste0(
  "Table S7. Figure 6C over-representation analysis (ReactomePA::enrichPathway v1.42.0, organism=mouse, BH-adjusted), ",
  "ALL 734 Reactome terms tested. Input DEG list: DE_ADR_B_vs_A.tsv, padj<0.05, A-ADR1 excluded (main analysis), n=546 ",
  "genes (543/546 = 99.5% mapped to Entrez ID). Background/universe: all 19,662 genes tested in this comparison (19,449/19,662 = 98.9% mapped)."
), s7)
add_toc("S7", "Fig. 6C ORA full results (DESeq2-derived DEG list, all 734 Reactome terms)",
        "25_Fig6C_ORA_rerun_DESeq2_DEGs.R", nrow(s7))

# ============================================================================
# Table S8: A-ADR1 (A1) sensitivity analysis
# ============================================================================
a1_wb_path <- file.path(outdir_main, "tables/Supplementary_Step5_A1_sensitivity.xlsx")
a1_caption <- "Table S8. A-ADR1 (A1, flagged low-purity sample) sensitivity analysis: main analysis (A1 excluded) vs. sensitivity (A1 included), for significant-gene counts, key genes (Wt1, Nphs1, Tmem215, Nlrp1b, Glp1r, Serpine1, Col4a1, Col4a2, Loxl1), and focal GSEA gene sets."
add_sheet("TableS8_A1sens_sig_counts", paste0(a1_caption, " [significant gene counts by comparison]"),
          read.xlsx(a1_wb_path, sheet = "overview_sig_counts"))
add_sheet("TableS8_A1sens_key_genes", paste0(a1_caption, " [key gene log2FC/padj, A1-excluded vs A1-included]"),
          read.xlsx(a1_wb_path, sheet = "key_genes_compare"))
add_sheet("TableS8_A1sens_GSEA", paste0(a1_caption, " [focal GSEA gene sets, A1-excluded vs A1-included]"),
          read.xlsx(a1_wb_path, sheet = "GSEA_compare"))
add_toc("S8", "A-ADR1 sensitivity analysis (sig. counts, key genes, focal GSEA; 3 sub-sheets)",
        "06_A1_sensitivity.R / 09_gsea_stat_ranking_canonical.R", NA)

# ============================================================================
# Table S9: A-Ctrl3 sensitivity analysis (new QC finding, not reviewer-raised)
# ============================================================================
actrl3_wb_path <- file.path(outdir_main, "tables/Supplementary_Step19_ACtrl3_sensitivity.xlsx")
actrl3_caption <- "Table S9. A-Ctrl3 (new QC finding: 2nd-highest tubular-marker contamination of all 12 samples, not previously flagged) sensitivity analysis: baseline main analysis (A-Ctrl3 included) vs. sensitivity (A-Ctrl3 excluded), for significant-gene counts, key genes, and baseline focal GSEA gene sets. All 6 focal gene sets retain the same sign with A-Ctrl3 excluded (see FigS_ACtrl3_sensitivity_NES_comparison.pdf)."
add_sheet("TableS9_ACtrl3sens_sig_counts", paste0(actrl3_caption, " [significant gene counts]"),
          read.xlsx(actrl3_wb_path, sheet = "overview_sig_counts"))
add_sheet("TableS9_ACtrl3sens_key_genes", paste0(actrl3_caption, " [key gene log2FC/padj]"),
          read.xlsx(actrl3_wb_path, sheet = "key_genes_compare"))
add_sheet("TableS9_ACtrl3sens_GSEA", paste0(actrl3_caption, " [baseline focal GSEA gene sets]"),
          read.xlsx(actrl3_wb_path, sheet = "GSEA_baseline_compare"))
add_toc("S9", "A-Ctrl3 sensitivity analysis (sig. counts, key genes, baseline focal GSEA; 3 sub-sheets)",
        "19_ACtrl3_sensitivity_dds_and_DE.R", NA)

# ============================================================================
# Table S10: ADR-response concordance between substrains (gene-level)
# ============================================================================
s10 <- read.csv(file.path(outdir_v2, "tables/TableS_ADR_response_A_vs_B_merged.csv"), stringsAsFactors = FALSE)
add_sheet("TableS10_ADR_concordance", paste0(
  "Table S10. Gene-level comparison of the ADR-response log2FC (ADR vs Ctrl) between substrains: log2FC_A ",
  "(within-AJcl) vs. log2FC_B (within-ByJcl), all 19,662 genes present in both comparisons. Pearson r=0.5354, ",
  "Spearman rho=0.4253. Underlies FigS_ADR_response_concordance.png."
), s10)
add_toc("S10", "Gene-level ADR-response concordance between substrains (merged log2FC table)",
        "24_manuscript_values_extraction.R", nrow(s10))

# ============================================================================
# Table of contents sheet (inserted first)
# ============================================================================
toc_df <- do.call(rbind, toc_rows)
addWorksheet(wb, "TOC", tabColour = "#0072B2")
writeData(wb, "TOC", "Supplementary Tables -- Table of Contents", startRow = 1, startCol = 1)
addStyle(wb, "TOC", createStyle(fontSize = 14, textDecoration = "bold"), rows = 1, cols = 1)
writeData(wb, "TOC", paste0(
  "GSE320498 BALB/cAJcl vs BALB/cByJcl adriamycin nephropathy re-analysis. ",
  "All tables reuse existing DESeq2 v1.38.3 / fgsea v1.24.0 / ReactomePA v1.42.0 output; ",
  "only Table S7 (Fig. 6C ORA) involved a new analysis run for this task. See ../manuscript_values.md ",
  "for the full text-reconciliation table and ../README_analysis_log.md for the complete methodological history."
), startRow = 2, startCol = 1)
mergeCells(wb, "TOC", cols = 1:4, rows = 2)
addStyle(wb, "TOC", createStyle(wrapText = TRUE, fontSize = 9, textDecoration = "italic"), rows = 2, cols = 1:4)
setRowHeights(wb, "TOC", rows = 2, heights = 45)
writeDataTable(wb, "TOC", toc_df, startRow = 4, tableStyle = "TableStyleLight9")
setColWidths(wb, "TOC", cols = 1:4, widths = c(6, 70, 45, 10))
worksheetOrder(wb) <- c(which(names(wb) == "TOC"), which(names(wb) != "TOC"))

out_fp <- file.path(outdir_v2, "tables/Supplementary_Tables_S1-S10.xlsx")
saveWorkbook(wb, out_fp, overwrite = TRUE)
cat("\nSaved:", out_fp, sprintf("(%d sheets total)\n", length(names(wb))))

cat("\n=== 26_Supplementary_Tables_consolidated.R complete ===\n")
