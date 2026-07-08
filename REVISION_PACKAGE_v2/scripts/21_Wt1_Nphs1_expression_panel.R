#!/usr/bin/env Rscript
# Companion expression panel for Wt1/Nphs1 (proposed Fig. 5B / Fig. 6A
# companion): the volcano plots (Fig. 5, Fig. 6A) already label Wt1/Nphs1
# with their exact log2FC/FDR, but a reader cannot see the underlying
# per-sample expression distribution from a volcano point alone. This script
# plots DESeq2-normalized counts per sample, grouped by substrain x
# treatment, for both genes at both time points in one figure, so the
# baseline (Fig. 5B) and Day-5 (Fig. 6A companion) claims are both directly
# checkable against the raw normalized expression.
#
# Uses the existing dds_group_sens.rds (all 12 samples, DESeq2 size-factor
# normalized) purely for plotting values -- no new statistical model. The
# significance brackets shown are the already-established main-analysis
# padj values (baseline: dds_group_main / A-ADR1 excluded; Day 5: same) from
# DE_baseline_B_vs_A.tsv and DE_ADR_B_vs_A.tsv. A-ADR1 and A-Ctrl3 (both
# flagged for elevated tubular contamination; see FigS_A1_contamination_QC.pdf
# and FigS_ACtrl3_contamination_QC.pdf) are shown as open/outlined points for
# transparency, not excluded from the plot.

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

df <- as.data.frame(t(nc[genes, , drop = FALSE]))
df$sample <- rownames(df)
df <- merge(df, coldata[, c("sample", "substrain", "treatment", "group")], by = "sample")

long <- df %>%
  tidyr::pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "norm_count")

long$group_label <- factor(
  paste0(ifelse(long$substrain == "A", "AJcl", "ByJcl"), "-", long$treatment),
  levels = c("AJcl-Ctrl", "ByJcl-Ctrl", "AJcl-ADR", "ByJcl-ADR")
)
long$substrain_label <- factor(ifelse(long$substrain == "A", "AJcl", "ByJcl"), levels = c("AJcl", "ByJcl"))
long$flagged <- long$sample %in% c("A-ADR1", "A-Ctrl3")
long$flag_label <- ifelse(long$sample == "A-ADR1", "A-ADR1 (excluded, main analysis)",
                    ifelse(long$sample == "A-Ctrl3", "A-Ctrl3 (flagged, sensitivity-checked)", "included"))
long$gene <- factor(long$gene, levels = genes)

# ---- significance annotations, from the established main-analysis DE tables ----
base_de <- read.delim(file.path(outdir_main, "tables/DE_baseline_B_vs_A.tsv"), stringsAsFactors = FALSE)
d5_de   <- read.delim(file.path(outdir_main, "tables/DE_ADR_B_vs_A.tsv"), stringsAsFactors = FALSE)
get_padj <- function(de, g) de$padj[de$gene == g]

fmt_p <- function(p) sprintf("FDR = %.3f%s", p, ifelse(p < 0.05, " *", ""))
ann <- data.frame(
  gene = rep(genes, each = 2),
  x = rep(c(1.5, 3.5), times = 2),
  xend = rep(c(1.5, 3.5), times = 2)
)

max_by_gene <- long %>% group_by(gene) %>% summarise(ymax = max(norm_count) * 1.18, .groups = "drop")

bracket_df <- do.call(rbind, lapply(genes, function(g) {
  ymax <- max_by_gene$ymax[max_by_gene$gene == g]
  rbind(
    data.frame(gene = g, x = 1, xend = 2, y = ymax, label = fmt_p(get_padj(base_de, g)), xmid = 1.5),
    data.frame(gene = g, x = 3, xend = 4, y = ymax, label = fmt_p(get_padj(d5_de, g)), xmid = 3.5)
  )
}))
bracket_df$gene <- factor(bracket_df$gene, levels = genes)

p <- ggplot(long, aes(x = group_label, y = norm_count)) +
  geom_boxplot(aes(fill = substrain_label), outlier.shape = NA, width = 0.55, alpha = 0.35) +
  geom_jitter(aes(fill = substrain_label, shape = flagged), width = 0.12, size = 3, stroke = 0.9, color = "black") +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 24), guide = "none") +
  scale_fill_manual(values = c("AJcl" = "#0072B2", "ByJcl" = "#CC79A7")) +
  facet_wrap(~gene, scales = "free_y", nrow = 1) +
  geom_segment(data = bracket_df, aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.5) +
  geom_text(data = bracket_df, aes(x = xmid, y = y, label = label), inherit.aes = FALSE, vjust = -0.6, size = 3.6, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22)), labels = scales::comma) +
  labs(
    title = "Wt1 and Nphs1 normalized expression, all groups and time points",
    subtitle = paste0(
      "Proposed Fig. 5B / Fig. 6A companion panel. Triangles: A-ADR1 (excluded from main analysis) and\n",
      "A-Ctrl3 (flagged, sensitivity-checked but included). Brackets show main-analysis FDR (DESeq2 Wald\n",
      "test): left bracket = baseline (Fig. 5), right bracket = Day 5 post-ADR (Fig. 6A)."
    ),
    x = NULL, y = "DESeq2-normalized counts", fill = "Substrain"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 8.7, lineheight = 1.15),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(size = 10.5, angle = 20, hjust = 1),
    axis.text.y = element_text(size = 10.5),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  )

for (ext in c("png", "pdf")) {
  fp <- file.path(outdir_v2, "figures", paste0("Fig5B_Wt1_Nphs1_expression.", ext))
  ggsave(fp, plot = p, width = 11, height = 6.5, dpi = 300)
  cat("Saved:", fp, "\n")
}

cat("\n=== 21_Wt1_Nphs1_expression_panel.R complete ===\n")
