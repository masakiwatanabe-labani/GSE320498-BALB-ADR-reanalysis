#!/usr/bin/env Rscript
# Step 1 (R1-5, R2-4): supplementary DE tables for the three primary contrasts.
# Uses the "main" dataset (A-ADR1 excluded; 11 samples) -- see 01_build_dds.R.
# Sign convention (fixed, per instructions):
#   baseline_B_vs_A : positive log2FC = higher in BALB/cByJcl (B) vs BALB/cAJcl (A) at baseline
#   A_ADR_vs_Ctrl   : positive log2FC = higher in ADR-treated vs Ctrl, WITHIN AJcl (A)
#   B_ADR_vs_Ctrl   : positive log2FC = higher in ADR-treated vs Ctrl, WITHIN ByJcl (B)

suppressMessages({
  library(DESeq2)
  library(openxlsx)
})

indir <- "/usr/local/jupyter/ADR_BALB"
outdir <- file.path(indir, "outputs")

dds <- readRDS(file.path(outdir, "tables/dds_group_main.rds"))

get_de_table <- function(dds, contrast, caption) {
  res <- results(dds, contrast = contrast, alpha = 0.05)
  df <- as.data.frame(res)
  df$gene <- rownames(df)
  df <- df[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
  colnames(df)[colnames(df) == "log2FoldChange"] <- "log2FC"
  df <- df[order(df$pvalue), ]
  attr(df, "caption") <- caption
  df
}

contrasts <- list(
  baseline_B_vs_A = list(
    contrast = c("group", "B_Ctrl", "A_Ctrl"),
    caption = "Baseline (untreated) comparison: BALB/cByJcl (B) vs BALB/cAJcl (A). Positive log2FC = higher expression in ByJcl. Samples: A-Ctrl1-3 (n=3) vs B-Ctrl1-3 (n=3). Design: ~group, DESeq2 Wald test, BH-adjusted p (padj/FDR). Main analysis dataset (A-ADR1 not applicable to this contrast). Reproduces/updates original manuscript Fig. 5."
  ),
  ADR_B_vs_A = list(
    contrast = c("group", "B_ADR", "A_ADR"),
    caption = "Day-5 post-ADR comparison: BALB/cByJcl (B) vs BALB/cAJcl (A), among ADR-treated animals. Positive log2FC = higher expression in ByJcl at Day 5 after ADR. Samples: A-ADR2,3 (n=2; A-ADR1 excluded, low glomerular purity, see Step 0/Step 5 sensitivity) vs B-ADR1-3 (n=3). Design: ~group, DESeq2 Wald test, BH-adjusted p (padj/FDR). Reproduces/updates original manuscript Fig. 6A; this is the contrast underlying the originally reported podocyte-signature and Integrin-signaling GSEA NES values."
  ),
  A_ADR_vs_Ctrl = list(
    contrast = c("group", "A_ADR", "A_Ctrl"),
    caption = "ADR response within BALB/cAJcl (A): ADR-treated vs Ctrl. Positive log2FC = higher expression after ADR. Samples: A-Ctrl1-3 (n=3) vs A-ADR2,3 (n=2; A-ADR1 excluded as the flagged low-purity sample, see Step 0/Step 5 sensitivity analysis). Design: ~group, DESeq2 Wald test, BH-adjusted p (padj/FDR)."
  ),
  B_ADR_vs_Ctrl = list(
    contrast = c("group", "B_ADR", "B_Ctrl"),
    caption = "ADR response within BALB/cByJcl (B): ADR-treated vs Ctrl. Positive log2FC = higher expression after ADR. Samples: B-Ctrl1-3 (n=3) vs B-ADR1-3 (n=3). Design: ~group, DESeq2 Wald test, BH-adjusted p (padj/FDR)."
  )
)

wb <- createWorkbook()
summary_lines <- c()

for (nm in names(contrasts)) {
  cfg <- contrasts[[nm]]
  df <- get_de_table(dds, cfg$contrast, cfg$caption)

  # tsv
  write.table(df, file.path(outdir, paste0("tables/DE_", nm, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # xlsx sheet: caption row(s) then table
  addWorksheet(wb, nm)
  writeData(wb, nm, cfg$caption, startRow = 1, startCol = 1)
  mergeCells(wb, nm, cols = 1:7, rows = 1)
  addStyle(wb, nm, createStyle(textDecoration = "italic", wrapText = TRUE), rows = 1, cols = 1:7)
  setRowHeights(wb, nm, rows = 1, heights = 45)
  writeDataTable(wb, nm, df, startRow = 3, tableStyle = "TableStyleLight9")
  setColWidths(wb, nm, cols = 1:7, widths = c(16, 12, 12, 10, 10, 12, 12))

  n_sig <- sum(df$padj < 0.05, na.rm = TRUE)
  summary_lines <- c(summary_lines,
    sprintf("%s: %d genes tested, %d with padj<0.05 (FDR)", nm, nrow(df), n_sig))
}

saveWorkbook(wb, file.path(outdir, "tables/Supplementary_DE_tables.xlsx"), overwrite = TRUE)

writeLines(summary_lines, file.path(outdir, "logs/02_DE_tables_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
cat("\nWrote: outputs/tables/DE_<contrast>.tsv (x3) and outputs/tables/Supplementary_DE_tables.xlsx\n")
