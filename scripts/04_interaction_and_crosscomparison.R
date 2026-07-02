#!/usr/bin/env Rscript
# Step 3 (R1-5 latter half): does the ADR response differ between substrains?
#  (a) substrain:treatment interaction DE table (genes whose ADR response
#      differs significantly between AJcl and ByJcl)
#  (b) scatter of A-response (A_ADR_vs_Ctrl log2FC) vs B-response
#      (B_ADR_vs_Ctrl log2FC), correlation + regression, concordance stats
#
# NOTE on statistical power asymmetry: A_ADR_vs_Ctrl uses n=2 (A-ADR1
# excluded) vs B_ADR_vs_Ctrl n=3. This asymmetry inflates the count of
# "B-only significant" genes for reasons unrelated to true biological
# differences, and is reported explicitly below and in the interaction
# table caption so it isn't mistaken for evidence of a substrain-specific
# response by itself.

suppressMessages({
  library(DESeq2)
  library(ggplot2)
  library(openxlsx)
})

indir <- "/usr/local/jupyter/ADR_BALB"
outdir <- file.path(indir, "outputs")

# ---- (a) interaction term DE table ----
dds_int <- readRDS(file.path(outdir, "tables/dds_int_main.rds"))
res_int <- results(dds_int, name = "substrainB.treatmentADR", alpha = 0.05)
df_int <- as.data.frame(res_int)
df_int$gene <- rownames(df_int)
df_int <- df_int[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
colnames(df_int)[colnames(df_int) == "log2FoldChange"] <- "log2FC_interaction"
df_int <- df_int[order(df_int$pvalue), ]

caption_int <- paste(
  "Substrain x treatment interaction term (substrainB.treatmentADR) from the joint model",
  "~ substrain + treatment + substrain:treatment (main dataset, A-ADR1 excluded, n=11).",
  "log2FC_interaction = (B_ADR-B_Ctrl) - (A_ADR-A_Ctrl) log2 fold change, i.e. how much MORE",
  "(positive) or LESS (negative) a gene changes with ADR in ByJcl compared with its change in",
  "AJcl. A significant interaction (padj<0.05) indicates the ADR response itself differs",
  "between substrains for that gene; it does NOT by itself indicate direction within either substrain (see the two single-substrain response tables for that)."
)

write.table(df_int, file.path(outdir, "tables/DE_interaction_substrain_by_treatment.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

n_sig_int <- sum(df_int$padj < 0.05, na.rm = TRUE)

# ---- (b) cross-substrain response comparison ----
a_resp <- read.delim(file.path(outdir, "tables/DE_A_ADR_vs_Ctrl.tsv"), stringsAsFactors = FALSE)
b_resp <- read.delim(file.path(outdir, "tables/DE_B_ADR_vs_Ctrl.tsv"), stringsAsFactors = FALSE)

m <- merge(a_resp, b_resp, by = "gene", suffixes = c("_A", "_B"))
m <- m[!is.na(m$log2FC_A) & !is.na(m$log2FC_B), ]

pear <- cor.test(m$log2FC_A, m$log2FC_B, method = "pearson")
spear <- cor.test(m$log2FC_A, m$log2FC_B, method = "spearman", exact = FALSE)
fit <- lm(log2FC_B ~ log2FC_A, data = m)
fit_sum <- summary(fit)

sig_A <- !is.na(m$padj_A) & m$padj_A < 0.05
sig_B <- !is.na(m$padj_B) & m$padj_B < 0.05
sig_either <- sig_A | sig_B
sig_both <- sig_A & sig_B

concordant <- sign(m$log2FC_A) == sign(m$log2FC_B)

pct_concordant_all      <- 100 * mean(concordant)
pct_concordant_either    <- 100 * mean(concordant[sig_either])
pct_concordant_both      <- if (sum(sig_both) > 0) 100 * mean(concordant[sig_both]) else NA

log_lines <- c(
  "=== Step 3: cross-substrain ADR-response comparison (A_ADR_vs_Ctrl vs B_ADR_vs_Ctrl) ===",
  sprintf("Genes compared (non-NA log2FC in both): %d", nrow(m)),
  sprintf("Pearson r = %.3f (95%% CI %.3f-%.3f), p = %.3e",
          pear$estimate, pear$conf.int[1], pear$conf.int[2], pear$p.value),
  sprintf("Spearman rho = %.3f, p = %.3e", spear$estimate, spear$p.value),
  sprintf("OLS regression: B_response = %.3f + %.3f * A_response, R2 = %.3f, slope p = %.3e",
          coef(fit)[1], coef(fit)[2], fit_sum$r.squared, fit_sum$coefficients[2,4]),
  sprintf("Genes significant (padj<0.05) in A_ADR_vs_Ctrl (n=2 vs 3): %d", sum(sig_A)),
  sprintf("Genes significant (padj<0.05) in B_ADR_vs_Ctrl (n=3 vs 3): %d", sum(sig_B)),
  sprintf("Genes significant in BOTH: %d", sum(sig_both)),
  sprintf("Genes significant in EITHER: %d", sum(sig_either)),
  sprintf("%% direction-concordant among ALL compared genes: %.1f%%", pct_concordant_all),
  sprintf("%% direction-concordant among genes significant in EITHER comparison: %.1f%%", pct_concordant_either),
  sprintf("%% direction-concordant among genes significant in BOTH comparisons: %s",
          if (is.na(pct_concordant_both)) "NA (no gene reached padj<0.05 in both)" else sprintf("%.1f%%", pct_concordant_both)),
  sprintf("Interaction term (substrain:treatment): %d / %d genes with padj<0.05", n_sig_int, nrow(df_int)),
  "",
  "Interpretation caveat: A_ADR_vs_Ctrl (n=2 vs 3, A-ADR1 excluded) has less statistical power",
  "than B_ADR_vs_Ctrl (n=3 vs 3), so the much larger number of significant genes in B is expected",
  "to partly reflect power asymmetry, not necessarily a larger true biological response.",
  "The genome-wide correlation/regression and the interaction-term test are the appropriate",
  "summaries for 'does response differ between substrains', rather than simply counting DE genes."
)
writeLines(log_lines, file.path(outdir, "logs/04_crosscomparison_summary.txt"))
cat(paste(log_lines, collapse = "\n"), "\n")

# scatter plot
p <- ggplot(m, aes(x = log2FC_A, y = log2FC_B)) +
  geom_point(aes(color = sig_either), size = 1.6, alpha = 0.55) +
  scale_color_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "grey70"),
                      labels = c(`TRUE` = "padj<0.05 in A or B response", `FALSE` = "n.s. in both"),
                      name = NULL) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "steelblue") +
  geom_hline(yintercept = 0, color = "grey85") + geom_vline(xintercept = 0, color = "grey85") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5,
           label = sprintf("Pearson r = %.2f, p = %.1e\nSpearman rho = %.2f\nslope = %.2f, R2 = %.2f",
                            pear$estimate, pear$p.value, spear$estimate, coef(fit)[2], fit_sum$r.squared),
           size = 4.5) +
  labs(title = "ADR response concordance between substrains",
       subtitle = "x = log2FC (A: ADR vs Ctrl), y = log2FC (B: ADR vs Ctrl)\ndashed line = identity (slope=1), solid line = OLS fit",
       x = "log2FC, AJcl ADR vs Ctrl", y = "log2FC, ByJcl ADR vs Ctrl") +
  theme_bw(base_size = 15) +
  theme(axis.title = element_text(face = "bold"))

ggsave(file.path(outdir, "figures/crosscomparison_ADR_response_scatter.png"),
       plot = p, width = 9, height = 7.5, dpi = 300)

write.table(m, file.path(outdir, "tables/crosscomparison_ADR_response_merged.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# append interaction + crosscomparison sheets to the supplementary workbook
wb_path <- file.path(outdir, "tables/Supplementary_DE_tables.xlsx")
wb <- loadWorkbook(wb_path)
for (nm in c("interaction_substrainxtreat", "crosscomparison_summary")) {
  if (nm %in% names(wb)) removeWorksheet(wb, nm)
}

addWorksheet(wb, "interaction_substrainxtreat")
writeData(wb, "interaction_substrainxtreat", caption_int, startRow = 1, startCol = 1)
mergeCells(wb, "interaction_substrainxtreat", cols = 1:7, rows = 1)
addStyle(wb, "interaction_substrainxtreat", createStyle(textDecoration = "italic", wrapText = TRUE), rows = 1, cols = 1:7)
setRowHeights(wb, "interaction_substrainxtreat", rows = 1, heights = 60)
writeDataTable(wb, "interaction_substrainxtreat", df_int, startRow = 3, tableStyle = "TableStyleLight9")

addWorksheet(wb, "crosscomparison_summary")
writeData(wb, "crosscomparison_summary", data.frame(summary = log_lines))

saveWorkbook(wb, wb_path, overwrite = TRUE)
cat("\nUpdated:", wb_path, "with interaction table and crosscomparison summary sheets\n")
