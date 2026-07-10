#!/usr/bin/env Rscript
# Group mean +/- SEM summary underlying Figure 7's bars (Serpine1, Loxl1,
# Col4a1, Col4a2; AJcl vs ByJcl, baseline vs Day 5). Pure summarization of
# the already-exported per-sample CPM table (TableS_ECM_genes_CPM.csv, from
# 29_Figure7_ECM_genes_CPM.R) -- no new data, no new model.
#
# A-ADR1 is excluded from the AJcl Day-5 mean/SEM (n=2: A-ADR2, A-ADR3),
# matching the main analysis and the figure itself, where A-ADR1 is shown
# separately as an orange triangle rather than folded into the bar.

indir <- "/usr/local/jupyter/ADR_BALB"
outdir_v2 <- file.path(indir, "outputs/REVISION_PACKAGE/REVISION_PACKAGE_v2")

df <- read.csv(file.path(outdir_v2, "tables/TableS_ECM_genes_CPM.csv"), stringsAsFactors = FALSE)
genes <- c("Serpine1", "Loxl1", "Col4a1", "Col4a2")
df_main <- df[df$sample != "A-ADR1", ]  # A-ADR1 excluded from group means, main analysis convention

summary_rows <- list()
for (g in genes) {
  agg <- aggregate(df_main[[g]], by = list(substrain = df_main$substrain, timepoint = df_main$timepoint),
                    FUN = function(x) c(mean = mean(x), sem = sd(x) / sqrt(length(x)), n = length(x)))
  res <- do.call(data.frame, agg)
  colnames(res) <- c("substrain", "timepoint", "mean_CPM", "SEM", "n")
  res$gene <- g
  res$mean_CPM <- round(res$mean_CPM, 2)
  res$SEM <- round(res$SEM, 2)
  summary_rows[[g]] <- res[, c("gene", "substrain", "timepoint", "mean_CPM", "SEM", "n")]
}
summary_df <- do.call(rbind, summary_rows)
summary_df$timepoint <- factor(summary_df$timepoint, levels = c("Baseline", "Day 5"))
summary_df$substrain <- factor(summary_df$substrain, levels = c("AJcl", "ByJcl"))
summary_df <- summary_df[order(match(summary_df$gene, genes), summary_df$timepoint, summary_df$substrain), ]

out_fp <- file.path(outdir_v2, "tables/TableS_ECM_genes_group_summary.csv")
write.csv(summary_df, out_fp, row.names = FALSE)
cat("Saved:", out_fp, "\n\n")
print(summary_df, row.names = FALSE)

# A-ADR1 individual values, for reference (shown separately in the figure)
a1_row <- df[df$sample == "A-ADR1", c("sample", genes)]
cat("\nA-ADR1 (excluded from AJcl Day-5 mean above; shown as separate triangle in Figure 7):\n")
print(a1_row, row.names = FALSE)

cat("\n=== 33_Figure7_group_summary.R complete ===\n")
