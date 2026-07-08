#!/usr/bin/env Rscript
# Companion expression panel for Wt1/Nphs1 (proposed Fig. 5B): the volcano
# plot (Fig. 5) already labels Wt1/Nphs1 with their exact log2FC/FDR, but a
# reader cannot see the underlying per-sample expression distribution from a
# volcano point alone. This script plots DESeq2-normalized counts per sample,
# grouped by substrain, for the two baseline (Ctrl) groups only -- matching
# Fig. 5's scope -- so the "no baseline difference" claim is directly
# checkable against the raw normalized expression.
#
# Uses the existing dds_group_sens.rds (DESeq2 size-factor normalized)
# purely for plotting values -- no new statistical model. The significance
# bracket is the already-established main-analysis padj (baseline_B_vs_A)
# from DE_baseline_B_vs_A.tsv. A-Ctrl3 (flagged for elevated tubular
# contamination; see FigS_ACtrl3_contamination_QC.pdf) is shown as an open
# triangle for transparency, not excluded from the plot.

suppressMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
})

indir <- "/usr/local/jupyter/ADR_BALB"
outdir_main <- file.path(indir, "outputs")
outdir_v2 <- file.path(outdir_main, "REVISION_PACKAGE/REVISION_PACKAGE_v2")

dds_sens <- readRDS(file.path(outdir_main, "tables/dds_group_sens.rds"))
nc <- counts(dds_sens, normalized = TRUE)
genes <- c("Wt1", "Nphs1")

coldata <- as.data.frame(colData(dds_sens))
coldata$sample <- rownames(coldata)

ctrl_samples <- rownames(coldata)[coldata$treatment == "Ctrl"]
df <- as.data.frame(t(nc[genes, ctrl_samples, drop = FALSE]))
df$sample <- rownames(df)
df <- merge(df, coldata[, c("sample", "substrain", "treatment", "group")], by = "sample")

long <- df %>%
  tidyr::pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "norm_count")

long$substrain_label <- factor(ifelse(long$substrain == "A", "AJcl", "ByJcl"), levels = c("AJcl", "ByJcl"))
long$flagged <- long$sample == "A-Ctrl3"
long$gene <- factor(long$gene, levels = genes)

# ---- significance annotation, from the established main-analysis DE table ----
base_de <- read.delim(file.path(outdir_main, "tables/DE_baseline_B_vs_A.tsv"), stringsAsFactors = FALSE)
get_padj <- function(de, g) de$padj[de$gene == g]
fmt_p <- function(p) sprintf("FDR = %.3f%s", p, ifelse(p < 0.05, " *", ""))

max_by_gene <- long %>% group_by(gene) %>% summarise(ymax = max(norm_count) * 1.18, .groups = "drop")
bracket_df <- do.call(rbind, lapply(genes, function(g) {
  data.frame(gene = g, x = 1, xend = 2, y = max_by_gene$ymax[max_by_gene$gene == g],
             label = fmt_p(get_padj(base_de, g)), xmid = 1.5)
}))
bracket_df$gene <- factor(bracket_df$gene, levels = genes)

p <- ggplot(long, aes(x = substrain_label, y = norm_count)) +
  geom_boxplot(aes(fill = substrain_label), outlier.shape = NA, width = 0.55, alpha = 0.35) +
  geom_jitter(aes(fill = substrain_label, shape = flagged), width = 0.1, size = 3.5, stroke = 0.9, color = "black") +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 24), guide = "none") +
  scale_fill_manual(values = c("AJcl" = "#0072B2", "ByJcl" = "#CC79A7")) +
  facet_wrap(~gene, scales = "free_y", nrow = 1) +
  geom_segment(data = bracket_df, aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.5) +
  geom_text(data = bracket_df, aes(x = xmid, y = y, label = label), inherit.aes = FALSE, vjust = -0.6, size = 4, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.22)), labels = scales::comma) +
  labs(
    title = "Wt1 and Nphs1 normalized expression at baseline",
    subtitle = paste0(
      "Proposed Fig. 5B companion panel. Triangle: A-Ctrl3 (flagged for elevated tubular\n",
      "contamination, sensitivity-checked but included). Bracket: main-analysis FDR (DESeq2 Wald test)."
    ),
    x = NULL, y = "DESeq2-normalized counts", fill = "Substrain"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, lineheight = 1.15),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(size = 11.5),
    axis.text.y = element_text(size = 10.5),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  )

for (ext in c("png", "pdf")) {
  fp <- file.path(outdir_v2, "figures", paste0("Fig5B_Wt1_Nphs1_expression.", ext))
  ggsave(fp, plot = p, width = 7.5, height = 6, dpi = 300)
  cat("Saved:", fp, "\n")
}

cat("\n=== 21_Wt1_Nphs1_expression_panel.R complete ===\n")
