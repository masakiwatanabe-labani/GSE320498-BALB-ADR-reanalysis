#!/usr/bin/env Rscript
# Figure 5, submission-ready final version: baseline (untreated) glomerular
# transcriptome volcano, BALB/cByJcl vs BALB/cAJcl, canonical DESeq2 run.
# EnhancedVolcano-style 4-category coloring (NS / log2FC only / FDR only /
# FDR and log2FC), Wt1/Nphs1 always labeled and visually distinguished
# (R1-5), 8 additional named high-effect genes labeled, all fonts >=8.5pt,
# vector PDF output. No new statistical model -- reads DE_baseline_B_vs_A.tsv
# (canonical: DESeq2 v1.38.3, all groups n=3) verbatim.

suppressMessages({
  library(ggplot2)
  library(ggrepel)
})

indir <- "/usr/local/jupyter/ADR_BALB"
outdir_main <- file.path(indir, "outputs")
outdir_v2 <- file.path(outdir_main, "REVISION_PACKAGE/REVISION_PACKAGE_v2")
fig_dir <- file.path(outdir_v2, "figures")
log_dir <- file.path(outdir_v2, "logs")

pt_mm <- function(pt) pt / 2.845276  # convert desired point size to ggplot's mm size unit

df <- read.delim(file.path(outdir_main, "tables/DE_baseline_B_vs_A.tsv"), stringsAsFactors = FALSE)
n_total <- nrow(df)

always_label <- c("Wt1", "Nphs1")
highlight_label <- c("Oas1a", "Oas1g", "H2-Q6", "Galk2", "Ppp1r3a", "Nlrp1b", "Tmem215", "Glp1r")
all_label <- c(always_label, highlight_label)

df$negLog10FDR <- -log10(df$padj)
df$category <- with(df, ifelse(is.na(padj), "NS",
                          ifelse(padj < 0.05 & abs(log2FC) > 1, "FDR and log2FC",
                          ifelse(padj < 0.05, "FDR only",
                          ifelse(abs(log2FC) > 1, "log2FC only", "NS")))))
df$category <- factor(df$category, levels = c("NS", "log2FC only", "FDR only", "FDR and log2FC"))

df$is_always <- df$gene %in% always_label
df$is_highlight <- df$gene %in% highlight_label
fmt_p <- function(p) ifelse(is.na(p), "NA", formatC(p, format = "e", digits = 2))
# gene-name-only labels (log2FC/FDR numbers deliberately omitted -- they are
# reported in the manuscript text and Table S1; keeping them off the figure
# avoids redundancy). fmt_p() is retained for the consistency-check log only.
df$point_label <- ifelse(df$gene %in% all_label, df$gene, NA)

# ---------------------------------------------------------------------------
# consistency checks (reported, not silently assumed)
# ---------------------------------------------------------------------------
check <- c(sprintf("Total tested genes (rows in DE_baseline_B_vs_A.tsv): %d (expected 19,662): %s",
                    n_total, ifelse(n_total == 19662, "MATCH", "*** MISMATCH ***")))

check <- c(check, sprintf("padj<0.05 significant genes: %d", sum(!is.na(df$padj) & df$padj < 0.05)))

wt1 <- df[df$gene == "Wt1", ]
nphs1 <- df[df$gene == "Nphs1", ]
check <- c(check,
  sprintf("Wt1: log2FC=%.4f, padj=%.4f (expected log2FC=+0.14, FDR=0.87): %s",
          wt1$log2FC, wt1$padj, ifelse(abs(wt1$log2FC - 0.14) < 0.01 && abs(wt1$padj - 0.87) < 0.01, "MATCH", "*** MISMATCH ***")),
  sprintf("Nphs1: log2FC=%.4f, padj=%.4f (expected log2FC=+0.23, FDR=0.79): %s",
          nphs1$log2FC, nphs1$padj, ifelse(abs(nphs1$log2FC - 0.23) < 0.01 && abs(nphs1$padj - 0.79) < 0.01, "MATCH", "*** MISMATCH ***"))
)

na_labeled <- df[df$gene %in% all_label & is.na(df$padj), ]
check <- c(check, sprintf("Labeled genes with padj=NA (should be zero): %d%s", nrow(na_labeled),
                           ifelse(nrow(na_labeled) > 0, paste0(" *** ", paste(na_labeled$gene, collapse=","), " ***"), " -- OK")))

known_na_genes <- c("Trim75","Gm46851","Cryga","Dreh","Crabp2","Apol7c","Platr21","Zfp988","Fcrl1")
known_ns_genes <- c("Gm7592","Ccl11","Nr5a2","Atcay")
accidentally_labeled <- df[df$gene %in% c(known_na_genes, known_ns_genes) & df$gene %in% all_label, ]
check <- c(check, sprintf("Forbidden genes (NA-padj or n.s. top-|log2FC|) accidentally labeled: %d%s",
                           nrow(accidentally_labeled),
                           ifelse(nrow(accidentally_labeled) > 0, " *** ERROR ***", " -- OK")))

fgf7 <- df[df$gene == "Fgf7", ]
check <- c(check, sprintf(
  "NOTE: Fgf7 also qualifies as a top-|log2FC|, non-NA, significant gene (log2FC=%.2f, padj=%.2e) but was NOT in the requested label list -- not labeled, per exact spec given.",
  fgf7$log2FC, fgf7$padj))

writeLines(check, file.path(log_dir, "30_Figure5_final_consistency_check.txt"))
cat(paste(check, collapse = "\n"), "\n\n")

# ---------------------------------------------------------------------------
# plot
# ---------------------------------------------------------------------------
cat_colors <- c("NS" = "grey75", "log2FC only" = "#0072B2", "FDR only" = "#009E73", "FDR and log2FC" = "#D55E00")

df_bg <- df[!(df$gene %in% all_label), ]
df_hi <- df[df$gene %in% highlight_label, ]
df_always <- df[df$gene %in% always_label, ]

df_always_lbl <- df[df$gene %in% always_label, ]
df_always_lbl <- df_always_lbl[order(df_always_lbl$gene), ]
df_always_lbl$nudge_x <- ifelse(df_always_lbl$gene == "Wt1", -4.2, 3.6)
df_always_lbl$nudge_y <- c(6, 9)
df_always_lbl$text_color <- "#9B4A73"

df_highlight_lbl <- df[df$gene %in% highlight_label, ]
df_highlight_lbl$nudge_x <- 0
df_highlight_lbl$nudge_y <- 0
df_highlight_lbl$text_color <- "black"

df_labels_all <- rbind(df_highlight_lbl, df_always_lbl)

x_range <- range(df$log2FC, na.rm = TRUE)

p <- ggplot(df, aes(x = log2FC, y = negLog10FDR)) +
  geom_point(data = df_bg, aes(color = category), size = 1.1, alpha = 0.6) +
  scale_color_manual(values = cat_colors, name = NULL, drop = FALSE) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey88", linewidth = 0.3) +
  geom_point(data = df_hi, aes(fill = "Named high-effect genes"), shape = 23, color = "black", size = 2.6, stroke = 0.6) +
  geom_point(data = df_always, aes(fill = "Wt1 / Nphs1 (podocyte markers)"), shape = 24, color = "black", size = 3.4, stroke = 0.8) +
  scale_fill_manual(values = c("Named high-effect genes" = "#F0E442", "Wt1 / Nphs1 (podocyte markers)" = "#CC79A7"), name = NULL) +
  ggrepel::geom_text_repel(
    data = df_labels_all, aes(label = point_label, colour = I(text_color), segment.colour = I(text_color)),
    nudge_x = df_labels_all$nudge_x, nudge_y = df_labels_all$nudge_y,
    size = pt_mm(9.0), fontface = "bold", lineheight = 0.9,
    max.overlaps = Inf, box.padding = 0.9, point.padding = 0.6,
    min.segment.length = 0, segment.size = 0.35,
    force = 10, force_pull = 0.25, max.time = 5, max.iter = 60000, seed = 20260220
  ) +
  annotate("text", x = x_range[1], y = -Inf, label = "← higher in AJcl", hjust = 0, vjust = -0.6, size = pt_mm(9.0), fontface = "italic", color = "grey30") +
  annotate("text", x = x_range[2], y = -Inf, label = "higher in ByJcl →", hjust = 1, vjust = -0.6, size = pt_mm(9.0), fontface = "italic", color = "grey30") +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.08))) +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.06))) +
  labs(
    title = "Figure 5. Baseline glomerular transcriptome: BALB/cByJcl vs BALB/cAJcl",
    subtitle = "DESeq2 Wald test, n=3 per substrain. Positive log2FC = higher in ByJcl. Dashed lines: FDR=0.05, |log2FC|=1.",
    x = expression(log[2]~"fold change (ByJcl vs AJcl)"),
    y = expression(-log[10]~"(FDR, BH-adjusted "*italic(p)*")"),
    caption = sprintf("total = %s variables", format(n_total, big.mark = ","))
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(size = 12.5, face = "bold"),
    plot.subtitle = element_text(size = 9.5),
    plot.caption = element_text(size = 9, hjust = 0, face = "italic", margin = margin(t = 6)),
    axis.title = element_text(size = 11.5, face = "bold"),
    axis.text = element_text(size = 9.5),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.y = unit(1, "pt")
  ) +
  guides(color = guide_legend(override.aes = list(size = 2.4, alpha = 1), nrow = 1),
         fill = guide_legend(override.aes = list(size = 3), nrow = 1))

for (ext in c("pdf", "png")) {
  fp <- file.path(fig_dir, paste0("Figure5_volcano_baseline.", ext))
  ggsave(fp, plot = p, width = 9, height = 7.5, device = if (ext == "pdf") cairo_pdf else ext, dpi = 300)
  cat("Saved:", fp, "\n")
}

caption_text <- paste0(
  "Baseline glomerular transcriptomic differences between AJcl and ByJcl are limited. Volcano plot of ",
  "differential expression (glomerulus-enriched RNA-seq) from untreated AJcl and ByJcl mice (n = 3 per ",
  "substrain), analyzed with DESeq2 (ByJcl vs AJcl; positive log2 fold change = higher in ByJcl). ",
  "Horizontal line, FDR = 0.05; vertical lines, |log2FC| = 1. Podocyte markers Wt1 and Nphs1 are labeled ",
  "and show no significant differential expression; ", sum(!is.na(df$padj) & df$padj < 0.05), " of ",
  format(n_total, big.mark = ","), " tested genes were significant (FDR < 0.05)."
)
writeLines(caption_text, file.path(log_dir, "30_Figure5_caption_draft.txt"))
cat("\n=== Figure caption (draft) ===\n", caption_text, "\n")

cat("\n=== 30_Figure5_volcano_baseline_final.R complete ===\n")
