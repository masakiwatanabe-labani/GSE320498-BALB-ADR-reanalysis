#!/usr/bin/env Rscript
# Step 2 (R1-5, R1-6, readability): volcano plots for the three main contrasts.
# Wt1 and Nphs1 are ALWAYS labeled (so readers can verify the "no robust
# change" claim). Tmem215, Nlrp1b, Glp1r are additionally highlighted+labeled.
# Sign convention stated explicitly in each plot title (see 02_DE_tables.R).

suppressMessages({
  library(ggplot2)
  library(ggrepel)
})

indir <- "/usr/local/jupyter/ADR_BALB"
outdir <- file.path(indir, "outputs")

always_label <- c("Wt1", "Nphs1")
highlight_label <- c("Tmem215", "Nlrp1b", "Glp1r")
all_label <- c(always_label, highlight_label)

plot_volcano <- function(df, title, subtitle, file_out) {
  df$negLog10FDR <- -log10(df$padj)
  df$negLog10FDR[is.infinite(df$negLog10FDR)] <- max(df$negLog10FDR[is.finite(df$negLog10FDR)], na.rm = TRUE) + 1
  df$sig <- ifelse(!is.na(df$padj) & df$padj < 0.05, "FDR < 0.05", "n.s.")
  df$label_class <- ifelse(df$gene %in% always_label, "Wt1/Nphs1 (podocyte markers)",
                     ifelse(df$gene %in% highlight_label, "Highlighted (Tmem215/Nlrp1b/Glp1r)", "other"))

  df_bg <- df[!(df$gene %in% all_label), ]
  df_hi <- df[df$gene %in% all_label, ]

  p <- ggplot(df, aes(x = log2FC, y = negLog10FDR)) +
    geom_point(data = df_bg, aes(color = sig), size = 1.3, alpha = 0.55) +
    scale_color_manual(values = c("FDR < 0.05" = "#D55E00", "n.s." = "grey70")) +
    geom_point(data = df_hi, aes(fill = label_class), shape = 21, color = "black",
               size = 3.6, stroke = 0.8) +
    scale_fill_manual(values = c("Wt1/Nphs1 (podocyte markers)" = "#0072B2",
                                  "Highlighted (Tmem215/Nlrp1b/Glp1r)" = "#009E73")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = 0, linetype = "solid", color = "grey85") +
    ggrepel::geom_text_repel(data = df_hi, aes(label = gene), size = 5.2, fontface = "bold",
                              max.overlaps = Inf, box.padding = 0.6, min.segment.length = 0,
                              segment.color = "black") +
    labs(title = title, subtitle = subtitle,
         x = expression(log[2]~"fold change"), y = expression(-log[10]~"(FDR, BH-adjusted p)"),
         color = "Significance", fill = "Labeled genes") +
    theme_bw(base_size = 16) +
    theme(plot.title = element_text(size = 15, face = "bold"),
          plot.subtitle = element_text(size = 12),
          axis.title = element_text(size = 16, face = "bold"),
          axis.text = element_text(size = 13),
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 12))

  ggsave(file_out, plot = p, width = 8.5, height = 7, dpi = 300)
  cat("Saved:", file_out, "\n")
}

comparisons <- list(
  baseline_B_vs_A = list(
    title = "Baseline: BALB/cByJcl vs BALB/cAJcl (positive log2FC = higher in ByJcl)",
    subtitle = "DESeq2 Wald test; A-Ctrl n=3 vs B-Ctrl n=3; main dataset"
  ),
  ADR_B_vs_A = list(
    title = "Day 5 post-ADR: BALB/cByJcl vs BALB/cAJcl (positive log2FC = higher in ByJcl)",
    subtitle = "A-ADR n=2 (A-ADR1 excluded) vs B-ADR n=3; reproduces/updates original Fig. 6A"
  ),
  A_ADR_vs_Ctrl = list(
    title = "BALB/cAJcl: ADR vs Ctrl (positive log2FC = higher after ADR, within AJcl)",
    subtitle = "A-Ctrl n=3 vs A-ADR n=2 (A-ADR1 excluded, see Step 5 sensitivity)"
  ),
  B_ADR_vs_Ctrl = list(
    title = "BALB/cByJcl: ADR vs Ctrl (positive log2FC = higher after ADR, within ByJcl)",
    subtitle = "B-Ctrl n=3 vs B-ADR n=3"
  )
)

for (nm in names(comparisons)) {
  df <- read.delim(file.path(outdir, paste0("tables/DE_", nm, ".tsv")), stringsAsFactors = FALSE)
  cfg <- comparisons[[nm]]
  plot_volcano(df, cfg$title, cfg$subtitle, file.path(outdir, paste0("figures/volcano_", nm, ".png")))
}
