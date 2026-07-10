#!/usr/bin/env Rscript
# Assemble the final, submission-numbered Supplementary Tables (S1-S6) into
# submission_supplementary/tables/. Pure repackaging/renaming of existing,
# already-verified files -- no new analysis. Figures are handled separately
# (copied directly; several requested figures do not exist in this project
# and are reported as missing, not fabricated -- see chat report).

suppressMessages(library(openxlsx))

indir <- "/usr/local/jupyter/ADR_BALB"
outdir_main <- file.path(indir, "outputs")
outdir_v2 <- file.path(outdir_main, "REVISION_PACKAGE/REVISION_PACKAGE_v2")
sub_dir <- file.path(outdir_main, "REVISION_PACKAGE/submission_supplementary")
tab_out <- file.path(sub_dir, "tables")

check <- c("=== Step 4: content verification for each table ===", "")
addchk <- function(...) check <<- c(check, sprintf(...))

# ============================================================================
# Table S1: DE tables, all 4 primary comparisons
# ============================================================================
file.copy(file.path(outdir_main, "tables/TableS_DE_all_comparisons.xlsx"),
          file.path(tab_out, "Table_S1.xlsx"), overwrite = TRUE)
wb1 <- loadWorkbook(file.path(tab_out, "Table_S1.xlsx"))
addchk("Table S1 (DE tables): sheets = %s", paste(names(wb1), collapse = ", "))
addchk("  Source: outputs/tables/TableS_DE_all_comparisons.xlsx (already captioned, 4 sheets, verified in R1-5/R1-6 work)")

# ============================================================================
# Table S2: Karaiskos podocyte marker gene list (TOP50 + EXCLUSIVE)
# ============================================================================
gmt_lines <- readLines(file.path(outdir_main, "genesets/KARAISKOS2018_PODOCYTE.gmt"))
parse_gmt_line <- function(ln) {
  parts <- strsplit(ln, "\t")[[1]]
  list(name = parts[1], desc = parts[2], genes = parts[-(1:2)])
}
gmt_sets <- lapply(gmt_lines, parse_gmt_line)
names(gmt_sets) <- sapply(gmt_sets, function(x) x$name)

exclusive <- gmt_sets[["KARAISKOS2018_PODOCYTE_EXCLUSIVE"]]
top50 <- gmt_sets[["KARAISKOS2018_PODOCYTE_TOP50"]]

counts_rownames <- read.delim(file.path(outdir_main, "tables/00_merged_counts.tsv"), stringsAsFactors = FALSE, nrows = 0) |> colnames()
matrix_genes <- read.delim(file.path(outdir_main, "tables/00_merged_counts.tsv"), stringsAsFactors = FALSE)$Geneid
df_top50 <- data.frame(rank = seq_along(top50$genes), gene = top50$genes,
                        in_exclusive_panel = top50$genes %in% exclusive$genes,
                        detected_in_count_matrix = top50$genes %in% matrix_genes)
df_exclusive <- data.frame(gene = exclusive$genes, detected_in_count_matrix = exclusive$genes %in% matrix_genes)

wb2 <- createWorkbook()
cap2 <- paste0("Table S2. Karaiskos et al. 2018 (J Am Soc Nephrol 29:2060-2068, GEO GSE111107) podocyte marker gene lists. ",
               "Source: authors' Supplementary Table 3 (FindAllMarkers podocyte-cluster output, top 50 by avg log2FC; ",
               "49 of 50 detected in this project's count matrix -- Sept11 absent due to the Septin11 MGI symbol renaming) ",
               "and Supplementary Table 2 (12-gene cell-type-exclusive stringent core). See ../Karaiskos_geneset_provenance.md ",
               "in the parent REVISION_PACKAGE for full extraction provenance.")
addWorksheet(wb2, "TOP50")
writeData(wb2, "TOP50", cap2, startRow = 1, startCol = 1); mergeCells(wb2, "TOP50", cols = 1:4, rows = 1)
addStyle(wb2, "TOP50", createStyle(textDecoration = "italic", wrapText = TRUE, fontSize = 9), rows = 1, cols = 1:4)
setRowHeights(wb2, "TOP50", rows = 1, heights = 60)
writeDataTable(wb2, "TOP50", df_top50, startRow = 3, tableStyle = "TableStyleLight9")
addWorksheet(wb2, "EXCLUSIVE")
writeData(wb2, "EXCLUSIVE", cap2, startRow = 1, startCol = 1); mergeCells(wb2, "EXCLUSIVE", cols = 1:2, rows = 1)
addStyle(wb2, "EXCLUSIVE", createStyle(textDecoration = "italic", wrapText = TRUE, fontSize = 9), rows = 1, cols = 1)
setRowHeights(wb2, "EXCLUSIVE", rows = 1, heights = 60)
writeDataTable(wb2, "EXCLUSIVE", df_exclusive, startRow = 3, tableStyle = "TableStyleLight9")
saveWorkbook(wb2, file.path(tab_out, "Table_S2.xlsx"), overwrite = TRUE)
addchk("")
addchk("Table S2 (Karaiskos marker list): TOP50 n=%d genes, EXCLUSIVE n=%d genes", nrow(df_top50), nrow(df_exclusive))
addchk("  Source: outputs/genesets/KARAISKOS2018_PODOCYTE.gmt (no file literally named 'TableS_Karaiskos_podocyte_markers' existed -- built fresh from the canonical GMT, which is itself the verified source used throughout Fig. 6B)")

# ============================================================================
# Table S3: A1 include/exclude sensitivity, focal gene sets, both settings
# ============================================================================
judgment <- read.delim(file.path(outdir_main, "tables/GSEA_canonical_focal_judgment_table.tsv"), stringsAsFactors = FALSE)
wb3 <- createWorkbook()
cap3 <- paste0("Table S3. A-ADR1 (A1) inclusion/exclusion sensitivity analysis for all 6 focal gene sets, canonical stat-ranked ",
               "fgsea, all comparisons where A1 status is applicable (ADR_B_vs_A, A_ADR_vs_Ctrl) plus the two comparisons where ",
               "A1 does not apply (baseline_B_vs_A, B_ADR_vs_Ctrl; single A1_excluded_main row each, shown for completeness). ",
               "Joint BH FDR computed across the full tested gene-set universe (1293-1294 sets) per comparison/A1-status combination.")
addWorksheet(wb3, "A1_sensitivity_focal_sets")
writeData(wb3, "A1_sensitivity_focal_sets", cap3, startRow = 1, startCol = 1)
mergeCells(wb3, "A1_sensitivity_focal_sets", cols = 1:10, rows = 1)
addStyle(wb3, "A1_sensitivity_focal_sets", createStyle(textDecoration = "italic", wrapText = TRUE, fontSize = 9), rows = 1, cols = 1:10)
setRowHeights(wb3, "A1_sensitivity_focal_sets", rows = 1, heights = 45)
writeDataTable(wb3, "A1_sensitivity_focal_sets", judgment, startRow = 3, tableStyle = "TableStyleLight9")
saveWorkbook(wb3, file.path(tab_out, "Table_S3.xlsx"), overwrite = TRUE)

ageing_check <- judgment[judgment$pathway == "TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING" & judgment$comparison == "ADR_B_vs_A", ]
addchk("")
addchk("Table S3 (A1 sensitivity, focal sets): %d rows, %d gene sets x %d comparison/A1-status combinations", nrow(judgment), length(unique(judgment$pathway)), length(unique(paste(judgment$comparison, judgment$a1_status))))
addchk("  Podocyte-ageing set, ADR_B_vs_A: NES=%.2f (A1-excluded) / NES=%.2f (A1-included) -- expected +1.93/-2.04: %s",
       ageing_check$NES[ageing_check$a1_status == "A1_excluded_main"],
       ageing_check$NES[ageing_check$a1_status == "A1_included_sens"],
       ifelse(round(ageing_check$NES[ageing_check$a1_status == "A1_excluded_main"], 2) == 1.93 &&
              round(ageing_check$NES[ageing_check$a1_status == "A1_included_sens"], 2) == -2.04, "MATCH", "*** CHECK ***"))

# ============================================================================
# Table S4: full GSEA (fgsea) + full ORA results
# ============================================================================
fgsea_full <- read.csv(file.path(outdir_v2, "tables/TableS_fgsea_ADR_B_vs_A_full.csv"), stringsAsFactors = FALSE)
ora_full <- read.csv(file.path(outdir_v2, "tables/TableS_ORA_D5_full.csv"), stringsAsFactors = FALSE)

wb4 <- createWorkbook()
cap4_fgsea <- paste0("Table S4a. Complete canonical preranked GSEA (fgsea v1.24.0, DESeq2 Wald-statistic ranking) results, ",
                      "ADR_B_vs_A (Day 5, A1-excluded, main analysis), ALL ", nrow(fgsea_full), " tested Reactome + podocyte gene sets.")
addWorksheet(wb4, "fgsea_ADR_B_vs_A_full")
writeData(wb4, "fgsea_ADR_B_vs_A_full", cap4_fgsea, startRow = 1, startCol = 1)
mergeCells(wb4, "fgsea_ADR_B_vs_A_full", cols = 1:7, rows = 1)
addStyle(wb4, "fgsea_ADR_B_vs_A_full", createStyle(textDecoration = "italic", wrapText = TRUE, fontSize = 9), rows = 1, cols = 1:7)
setRowHeights(wb4, "fgsea_ADR_B_vs_A_full", rows = 1, heights = 40)
writeDataTable(wb4, "fgsea_ADR_B_vs_A_full", fgsea_full, startRow = 3, tableStyle = "TableStyleLight9")

cap4_ora <- paste0("Table S4b. Complete Reactome ORA (ReactomePA::enrichPathway v1.42.0) results, DESeq2-derived Day-5 DEG list ",
                    "(padj<0.05, n=546, A1-excluded), background=19,662 tested genes, ALL ", nrow(ora_full), " tested Reactome terms.")
addWorksheet(wb4, "ORA_D5_full")
writeData(wb4, "ORA_D5_full", cap4_ora, startRow = 1, startCol = 1)
mergeCells(wb4, "ORA_D5_full", cols = 1:9, rows = 1)
addStyle(wb4, "ORA_D5_full", createStyle(textDecoration = "italic", wrapText = TRUE, fontSize = 9), rows = 1, cols = 1:9)
setRowHeights(wb4, "ORA_D5_full", rows = 1, heights = 40)
writeDataTable(wb4, "ORA_D5_full", ora_full, startRow = 3, tableStyle = "TableStyleLight9")
saveWorkbook(wb4, file.path(tab_out, "Table_S4.xlsx"), overwrite = TRUE)

addchk("")
addchk("Table S4 (fgsea full + ORA full): fgsea sheet n=%d rows (expected ~1293): %s | ORA sheet n=%d rows (expected 734): %s",
       nrow(fgsea_full), ifelse(nrow(fgsea_full) %in% 1290:1296, "MATCH", "*** CHECK ***"),
       nrow(ora_full), ifelse(nrow(ora_full) == 734, "MATCH", "*** CHECK ***"))

# ============================================================================
# Table S5: focal gene set results (Karaiskos exclusive / integrin / ECM)
# ============================================================================
focal5 <- judgment[judgment$a1_status == "A1_excluded_main" &
                    judgment$pathway %in% c("KARAISKOS2018_PODOCYTE_EXCLUSIVE",
                                             "REACTOME_INTEGRIN_SIGNALING",
                                             "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION"), ]
focal5 <- focal5[order(focal5$pathway, focal5$comparison), ]
wb5 <- createWorkbook()
cap5 <- paste0("Table S5. Canonical (A1-excluded, main analysis) focal gene-set results for the three gene sets used in ",
               "Figure 6 panels B (Karaiskos exclusive)/C (ECM organization theme)/D (Integrin signaling), across all 4 primary comparisons. ",
               "Positive NES = higher in ByJcl (baseline_B_vs_A, ADR_B_vs_A) or higher after ADR (A_ADR_vs_Ctrl, B_ADR_vs_Ctrl).")
addWorksheet(wb5, "focal_gene_sets")
writeData(wb5, "focal_gene_sets", cap5, startRow = 1, startCol = 1)
mergeCells(wb5, "focal_gene_sets", cols = 1:10, rows = 1)
addStyle(wb5, "focal_gene_sets", createStyle(textDecoration = "italic", wrapText = TRUE, fontSize = 9), rows = 1, cols = 1:10)
setRowHeights(wb5, "focal_gene_sets", rows = 1, heights = 45)
writeDataTable(wb5, "focal_gene_sets", focal5, startRow = 3, tableStyle = "TableStyleLight9")
saveWorkbook(wb5, file.path(tab_out, "Table_S5.xlsx"), overwrite = TRUE)
addchk("")
addchk("Table S5 (focal gene sets: Karaiskos exclusive/Integrin/ECM): %d rows (3 gene sets x 4 comparisons)", nrow(focal5))

# ============================================================================
# Table S6: interaction gene list (substrain x treatment, 24 genes)
# ============================================================================
interaction_genes <- read.csv(file.path(outdir_v2, "tables/TableS2_interaction_genes.csv"), stringsAsFactors = FALSE)
wb6 <- createWorkbook()
cap6 <- paste0("Table S6. Genes with a significant substrain x treatment interaction (DESeq2 ~substrain+treatment+substrain:treatment, ",
               "A-ADR1 excluded, canonical run), interaction FDR<0.05 (", nrow(interaction_genes), " of 19,662 tested genes). ",
               "log2FC_A/padj_A = within-AJcl ADR response; log2FC_B/padj_B = within-ByJcl ADR response; ",
               "log2FC_baseline/padj_baseline = baseline substrain difference for the same gene. ",
               "NOTE: this table was referred to as 'Table S2' during an earlier interim analysis step (interaction-gene ",
               "pathway-coherence check); it is renumbered Table S6 here to match the final submission numbering.")
addWorksheet(wb6, "interaction_genes")
writeData(wb6, "interaction_genes", cap6, startRow = 1, startCol = 1)
mergeCells(wb6, "interaction_genes", cols = 1:11, rows = 1)
addStyle(wb6, "interaction_genes", createStyle(textDecoration = "italic", wrapText = TRUE, fontSize = 9), rows = 1, cols = 1:11)
setRowHeights(wb6, "interaction_genes", rows = 1, heights = 60)
writeDataTable(wb6, "interaction_genes", interaction_genes, startRow = 3, tableStyle = "TableStyleLight9")
saveWorkbook(wb6, file.path(tab_out, "Table_S6.xlsx"), overwrite = TRUE)
addchk("")
addchk("Table S6 (interaction genes): n=%d (expected 24): %s", nrow(interaction_genes),
       ifelse(nrow(interaction_genes) == 24, "MATCH", "*** CHECK ***"))
addchk("  NOTE: this table was called 'Table S2' in the interim interaction-coherence task; renumbered to Table S6 here.")

writeLines(check, file.path(sub_dir, "table_content_check.txt"))
cat(paste(check, collapse = "\n"), "\n")
cat("\n=== 32_submission_supplementary_package.R (tables) complete ===\n")
