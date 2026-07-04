#!/usr/bin/env Rscript
# R1-5 (transcriptional-analysis description) final deliverables.
#
# R1-5 asks for: (1) explicit contrast sign on Fig.5/6A volcanoes with
# Wt1/Nphs1 always labeled and Tmem215/Nlrp1b/Glp1r highlighted with their
# numeric log2FC/FDR, (2) a gene-level description table replacing the
# "modest" wording for the baseline comparison, (3) a supplementary DE table
# for all 4 primary contrasts with sign convention captions, and (4) a
# cross-substrain ADR-response concordance figure/summary.
#
# Items (1)-(4) were already produced by 02_DE_tables.R, 03_volcano_plots.R,
# and 04_interaction_and_crosscomparison.R earlier in this pipeline. This
# script does not recompute anything -- it packages those exact numbers into
# the specific final filenames requested (Fig5_volcano_baseline.*,
# Fig6A_volcano_D5.*, TableS_top_baseline_genes.csv,
# TableS_DE_all_comparisons.xlsx, FigS_ADR_response_concordance.*), and adds
# the one thing that did not exist yet: numeric log2FC/FDR annotated directly
# on the volcano labels (not just in the caption) and the gene-description
# text summary for the baseline comparison.

suppressMessages({
  library(ggplot2)
  library(ggrepel)
  library(openxlsx)
})

indir <- "/usr/local/jupyter/ADR_BALB"
outdir <- file.path(indir, "outputs")

always_label <- c("Wt1", "Nphs1")
highlight_label <- c("Tmem215", "Nlrp1b", "Glp1r")
all_label <- c(always_label, highlight_label)

fmt_p <- function(p) {
  ifelse(is.na(p), "NA", formatC(p, format = "e", digits = 2))
}

# ---------------------------------------------------------------------------
# Step 1: Fig5 / Fig6A volcano plots, with per-gene log2FC/FDR printed on the
# in-figure labels themselves (not only in the caption).
# ---------------------------------------------------------------------------
plot_volcano_annotated <- function(df, title, subtitle, file_stem) {
  df$negLog10FDR <- -log10(df$padj)
  df$negLog10FDR[is.infinite(df$negLog10FDR)] <- max(df$negLog10FDR[is.finite(df$negLog10FDR)], na.rm = TRUE) + 1
  df$sig <- ifelse(!is.na(df$padj) & df$padj < 0.05, "FDR < 0.05", "n.s.")
  df$label_class <- ifelse(df$gene %in% always_label, "Wt1/Nphs1 (podocyte markers)",
                     ifelse(df$gene %in% highlight_label, "Highlighted (Tmem215/Nlrp1b/Glp1r)", "other"))
  df$point_label <- ifelse(df$gene %in% all_label,
                            sprintf("%s\nlog2FC=%.2f, FDR=%s", df$gene, df$log2FC, fmt_p(df$padj)),
                            NA)

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
    ggrepel::geom_text_repel(data = df_hi, aes(label = point_label), size = 3.8, fontface = "bold",
                              max.overlaps = Inf, box.padding = 1.6, point.padding = 0.6,
                              min.segment.length = 0, segment.color = "black", lineheight = 0.9,
                              force = 20, force_pull = 0.3, max.time = 5, max.iter = 50000,
                              seed = 42) +
    scale_y_continuous(expand = expansion(mult = c(0.22, 0.06))) +
    labs(title = title, subtitle = subtitle,
         x = expression(log[2]~"fold change"), y = expression(-log[10]~"(FDR, BH-adjusted p)"),
         color = "Significance", fill = "Labeled genes") +
    theme_bw(base_size = 16) +
    theme(plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 10),
          axis.title = element_text(size = 16, face = "bold"),
          axis.text = element_text(size = 13),
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 12))

  for (ext in c("png", "pdf")) {
    fp <- paste0(file_stem, ".", ext)
    ggsave(fp, plot = p, width = 11, height = 8, dpi = 300)
    cat("Saved:", fp, "\n")
  }
}

baseline <- read.delim(file.path(outdir, "tables/DE_baseline_B_vs_A.tsv"), stringsAsFactors = FALSE)
d5 <- read.delim(file.path(outdir, "tables/DE_ADR_B_vs_A.tsv"), stringsAsFactors = FALSE)

plot_volcano_annotated(
  baseline,
  "Baseline: BALB/cByJcl vs BALB/cAJcl (positive log2FC = higher in ByJcl)",
  "DESeq2 Wald test; A-Ctrl n=3 vs B-Ctrl n=3; main dataset. Reproduces/updates original manuscript Fig. 5.",
  file.path(outdir, "figures/Fig5_volcano_baseline")
)

plot_volcano_annotated(
  d5,
  "Day 5 post-ADR: BALB/cByJcl vs BALB/cAJcl (positive log2FC = higher in ByJcl)",
  "DESeq2 Wald test; A-ADR n=2 (A-ADR1 excluded) vs B-ADR n=3. Reproduces/updates original manuscript Fig. 6A.",
  file.path(outdir, "figures/Fig6A_volcano_D5")
)

# ---------------------------------------------------------------------------
# Step 2: baseline gene-description table + text summary (replaces "modest")
# ---------------------------------------------------------------------------
sig_thresh <- baseline[!is.na(baseline$padj) & baseline$padj < 0.05 & abs(baseline$log2FC) > 1, ]
sig_thresh <- sig_thresh[order(sig_thresh$padj), ]
sig_thresh$is_R1_5_highlighted_gene <- sig_thresh$gene %in% highlight_label
tableS_top <- sig_thresh[, c("gene", "log2FC", "padj", "baseMean", "is_R1_5_highlighted_gene")]
tableS_top$log2FC <- round(tableS_top$log2FC, 3)
tableS_top$padj <- signif(tableS_top$padj, 3)
tableS_top$baseMean <- round(tableS_top$baseMean, 1)

write.csv(tableS_top, file.path(outdir, "tables/TableS_top_baseline_genes.csv"), row.names = FALSE)
cat("Saved: tables/TableS_top_baseline_genes.csv (", nrow(tableS_top), "genes, padj<0.05 & |log2FC|>1 )\n")

get_row <- function(df, g) df[df$gene == g, ]
r_tmem215 <- get_row(baseline, "Tmem215")
r_nlrp1b  <- get_row(baseline, "Nlrp1b")
r_glp1r   <- get_row(baseline, "Glp1r")
r_wt1     <- get_row(baseline, "Wt1")
r_nphs1   <- get_row(baseline, "Nphs1")

n_total <- nrow(baseline)
n_sig_padj <- sum(!is.na(baseline$padj) & baseline$padj < 0.05, na.rm = TRUE)
n_sig_strict <- nrow(sig_thresh)

summary_lines <- c(
  "=== Step 2: baseline (BALB/cByJcl vs BALB/cAJcl) gene-description summary ===",
  sprintf("Genes tested (non-NA baseMean): %d", n_total),
  sprintf("Genes with padj<0.05: %d", n_sig_padj),
  sprintf("Genes with padj<0.05 AND |log2FC|>1 (strict, used for TableS_top_baseline_genes.csv): %d", n_sig_strict),
  "",
  "Genes specifically flagged by R1-5 (baseline comparison, positive log2FC = higher in ByJcl):",
  sprintf("  Tmem215: log2FC=%.2f, padj=%s (baseMean=%.0f) -- higher in ByJcl, one of the most significant genes genome-wide",
          r_tmem215$log2FC, fmt_p(r_tmem215$padj), r_tmem215$baseMean),
  sprintf("  Glp1r:   log2FC=%.2f, padj=%s (baseMean=%.0f) -- higher in ByJcl",
          r_glp1r$log2FC, fmt_p(r_glp1r$padj), r_glp1r$baseMean),
  sprintf("  Nlrp1b:  log2FC=%.2f, padj=%s (baseMean=%.0f) -- higher in AJcl (negative log2FC)",
          r_nlrp1b$log2FC, fmt_p(r_nlrp1b$padj), r_nlrp1b$baseMean),
  "",
  "Podocyte-identity genes (should NOT differ if substrains share baseline podocyte reserve):",
  sprintf("  Wt1:   log2FC=%.2f, padj=%s -- not significant, consistent with no baseline difference",
          r_wt1$log2FC, fmt_p(r_wt1$padj)),
  sprintf("  Nphs1: log2FC=%.2f, padj=%s -- not significant, consistent with no baseline difference",
          r_nphs1$log2FC, fmt_p(r_nphs1$padj)),
  "",
  sprintf("Text-summary replacement for 'modest': of %d genes tested at baseline, %d (%.1f%%) reach padj<0.05, and %d (%.1f%%) reach the stricter padj<0.05 & |log2FC|>1 threshold used in TableS_top_baseline_genes.csv -- the baseline substrain difference is real and detectable at many loci (including Tmem215, Glp1r, Nlrp1b at high significance) but is not accompanied by any change in core podocyte-identity markers (Wt1, Nphs1).",
          n_total, n_sig_padj, 100 * n_sig_padj / n_total, n_sig_strict, 100 * n_sig_strict / n_total)
)
writeLines(summary_lines, file.path(outdir, "logs/16_R1-5_baseline_gene_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

# ---------------------------------------------------------------------------
# Step 3: supplementary DE table, all 4 primary contrasts, final filename
# ---------------------------------------------------------------------------
wb_final <- loadWorkbook(file.path(outdir, "tables/Supplementary_DE_tables.xlsx"))
final_sheets <- c("baseline_B_vs_A", "ADR_B_vs_A", "A_ADR_vs_Ctrl", "B_ADR_vs_Ctrl")
extra_sheets <- setdiff(names(wb_final), final_sheets)
for (nm in extra_sheets) removeWorksheet(wb_final, nm)
saveWorkbook(wb_final, file.path(outdir, "tables/TableS_DE_all_comparisons.xlsx"), overwrite = TRUE)
cat("Saved: tables/TableS_DE_all_comparisons.xlsx (4 sheets: ", paste(final_sheets, collapse = ", "), ")\n")

# ---------------------------------------------------------------------------
# Step 4: cross-substrain ADR-response concordance figure, final filename
# ---------------------------------------------------------------------------
invisible(file.copy(file.path(outdir, "figures/crosscomparison_ADR_response_scatter.png"),
                     file.path(outdir, "figures/FigS_ADR_response_concordance.png"), overwrite = TRUE))
cat("Saved: figures/FigS_ADR_response_concordance.png (copy of crosscomparison_ADR_response_scatter.png)\n")

cat("\n=== 16_R1-5_final_deliverables.R complete ===\n")
