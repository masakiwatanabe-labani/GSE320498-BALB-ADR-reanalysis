#!/usr/bin/env Rscript
# DESeq2 statistics (baseline and Day 5, ByJcl vs AJcl) for the 4 genes shown
# in Figure 7 (Serpine1, Loxl1, Col4a1, Col4a2). Pure extraction from the
# canonical DE tables (DE_baseline_B_vs_A.tsv, DE_ADR_B_vs_A.tsv) -- no new
# statistical model, same numbers already annotated on Figure 7 itself and
# in its caption, just collected into one file alongside
# TableS_ECM_genes_CPM.csv and TableS_ECM_genes_group_summary.csv.

indir <- "/usr/local/jupyter/ADR_BALB"
outdir_main <- file.path(indir, "outputs")
outdir_v2 <- file.path(outdir_main, "REVISION_PACKAGE/REVISION_PACKAGE_v2")

genes <- c("Serpine1", "Loxl1", "Col4a1", "Col4a2")
cols <- c("gene", "baseMean", "log2FC", "lfcSE", "stat", "pvalue", "padj")

baseline <- read.delim(file.path(outdir_main, "tables/DE_baseline_B_vs_A.tsv"), stringsAsFactors = FALSE)
d5 <- read.delim(file.path(outdir_main, "tables/DE_ADR_B_vs_A.tsv"), stringsAsFactors = FALSE)

base_sub <- baseline[baseline$gene %in% genes, cols]
base_sub$timepoint <- "Baseline"
base_sub$comparison <- "baseline_B_vs_A (ByJcl vs AJcl, n=3 vs 3)"

d5_sub <- d5[d5$gene %in% genes, cols]
d5_sub$timepoint <- "Day 5"
d5_sub$comparison <- "ADR_B_vs_A (ByJcl vs AJcl, Day 5 post-ADR, A-ADR1 excluded, n=3 vs 2)"

stats_df <- rbind(base_sub, d5_sub)
stats_df$timepoint <- factor(stats_df$timepoint, levels = c("Baseline", "Day 5"))
stats_df <- stats_df[order(match(stats_df$gene, genes), stats_df$timepoint), ]
stats_df$sign_convention <- "positive log2FC = higher in ByJcl"
stats_df$significant_FDR0.05 <- ifelse(stats_df$padj < 0.05, "YES", "no")

# round for readability, keep full precision available via pvalue/padj sci notation
stats_df$baseMean <- round(stats_df$baseMean, 1)
stats_df$log2FC <- round(stats_df$log2FC, 4)
stats_df$lfcSE <- round(stats_df$lfcSE, 4)
stats_df$stat <- round(stats_df$stat, 3)
stats_df$pvalue <- signif(stats_df$pvalue, 4)
stats_df$padj <- signif(stats_df$padj, 4)

out_fp <- file.path(outdir_v2, "tables/TableS_ECM_genes_statistics.csv")
write.csv(stats_df[, c("gene", "timepoint", "comparison", "baseMean", "log2FC", "lfcSE",
                        "stat", "pvalue", "padj", "significant_FDR0.05", "sign_convention")],
          out_fp, row.names = FALSE)
cat("Saved:", out_fp, "\n\n")
print(stats_df[, c("gene", "timepoint", "log2FC", "padj", "significant_FDR0.05")], row.names = FALSE)

cat("\nNote: Loxl1 is the only one of the 4 genes already significant at baseline\n")
cat("(padj=0.0052), rising further by Day 5 (padj=4.86e-12) -- not a purely\n")
cat("Day-5-specific onset. Serpine1/Col4a1/Col4a2 are all non-significant at\n")
cat("baseline (padj>0.45) and significant by Day 5.\n")

cat("\n=== 35_Figure7_statistics.R complete ===\n")
