#!/usr/bin/env Rscript
# Figure 7: CPM of 4 ECM/adhesion-related genes (Serpine1, Loxl1, Col4a1,
# Col4a2) discussed in the text, at baseline and Day 5, both substrains --
# "no baseline difference -> ByJcl-specific rise at Day 5" at a glance.
#
# Plain publication-style bar graph (mean + SEM, no jittered individual
# points, no A-ADR1 marker) per user request -- the earlier version's
# per-sample dots and separately-flagged A-ADR1 triangle were judged too
# elaborate for a main-text figure. A-ADR1 is dropped from this figure
# entirely (not shown at all, not just excluded from the mean); it remains
# visible, flagged, in TableS_ECM_genes_CPM.csv and in the dedicated QC
# figure (FigS_glomerular_purity_QC.pdf) and Figure 6A volcano.
#
# CPM computed directly from tables/00_merged_counts.tsv (library-size CPM,
# same definition as 18/20/28_*_QC*.R), all 12 samples (A-ADR1 dropped before
# plotting). Significance annotations (asterisks) are read directly from
# DE_ADR_B_vs_A.tsv / DE_baseline_B_vs_A.tsv -- no new DE model.
#
# NOTE: no prior "Figure7_ECM_genes_CPM" file was found anywhere in the
# project at the time this script was first written (checked outputs/ and
# REVISION_PACKAGE/ and REVISION_PACKAGE_v2/ recursively) -- this is a new
# figure, not a re-derivation of an existing one. CPM definition matches the
# project's other CPM-based QC figures for consistency.

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

indir <- "/usr/local/jupyter/ADR_BALB"
outdir_main <- file.path(indir, "outputs")
outdir_v2 <- file.path(outdir_main, "REVISION_PACKAGE/REVISION_PACKAGE_v2")
fig_dir <- file.path(outdir_v2, "figures")
tab_dir <- file.path(outdir_v2, "tables")
log_dir <- file.path(outdir_v2, "logs")

pt_mm <- function(pt) pt / 2.845276

genes <- c("Serpine1", "Loxl1", "Col4a1", "Col4a2")

merged <- read.delim(file.path(outdir_main, "tables/00_merged_counts.tsv"),
                      check.names = FALSE, stringsAsFactors = FALSE)
rownames(merged) <- merged$Geneid
all_samples <- c("A-Ctrl1","A-Ctrl2","A-Ctrl3","B-Ctrl1","B-Ctrl2","B-Ctrl3",
                  "A-ADR1","A-ADR2","A-ADR3","B-ADR1","B-ADR2","B-ADR3")
counts_only <- merged[, all_samples]
cpm <- sweep(counts_only, 2, colSums(counts_only), FUN = "/") * 1e6

cpm_df <- as.data.frame(t(cpm[genes, all_samples, drop = FALSE]))
cpm_df$sample <- rownames(cpm_df)
cpm_df$substrain <- ifelse(grepl("^A-", cpm_df$sample), "AJcl", "ByJcl")
cpm_df$timepoint <- ifelse(grepl("Ctrl", cpm_df$sample), "Baseline", "Day 5")
cpm_df$is_A1 <- cpm_df$sample == "A-ADR1"

long <- cpm_df %>% pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "cpm")
long$gene <- factor(long$gene, levels = genes)
long$timepoint <- factor(long$timepoint, levels = c("Baseline", "Day 5"))
long$substrain <- factor(long$substrain, levels = c("AJcl", "ByJcl"))
long$group_x <- interaction(long$timepoint, long$substrain, sep = "\n")

write.csv(cpm_df[, c("sample", "substrain", "timepoint", genes)],
          file.path(tab_dir, "TableS_ECM_genes_CPM.csv"), row.names = FALSE)
cat("Saved:", file.path(tab_dir, "TableS_ECM_genes_CPM.csv"), "\n\n")

# mean +/- SEM per gene x timepoint x substrain, EXCLUDING A-ADR1 (matches
# the main-analysis sample set: AJcl Day-5 n=2)
summary_df <- long %>%
  filter(!is_A1) %>%
  group_by(gene, timepoint, substrain) %>%
  summarise(mean_cpm = mean(cpm), sem = sd(cpm) / sqrt(n()), n = n(), .groups = "drop")

# ---------------------------------------------------------------------------
# text-vs-canonical value check
# ---------------------------------------------------------------------------
baseline_de <- read.delim(file.path(outdir_main, "tables/DE_baseline_B_vs_A.tsv"), stringsAsFactors = FALSE)
d5_de <- read.delim(file.path(outdir_main, "tables/DE_ADR_B_vs_A.tsv"), stringsAsFactors = FALSE)
manuscript_d5 <- data.frame(gene = genes, log2FC_text = c(1.84, 1.38, 0.97, 1.22))

check_lines <- c("=== Text-vs-canonical check: Day 5 (ByJcl vs AJcl, A1-excluded, main) ===")
for (g in genes) {
  r <- d5_de[d5_de$gene == g, ]
  txt <- manuscript_d5$log2FC_text[manuscript_d5$gene == g]
  match_ok <- abs(r$log2FC - txt) < 0.02
  check_lines <- c(check_lines, sprintf("  %-10s canonical log2FC=%.4f (text=%.2f, %s)  padj=%.3e",
                                         g, r$log2FC, txt, ifelse(match_ok, "MATCH", "*** MISMATCH ***"), r$padj))
}
check_lines <- c(check_lines, "", "=== Baseline significance check (verify 'no baseline difference' claim) ===")
any_sig_baseline <- FALSE
for (g in genes) {
  r <- baseline_de[baseline_de$gene == g, ]
  sig <- !is.na(r$padj) & r$padj < 0.05
  if (sig) any_sig_baseline <- TRUE
  check_lines <- c(check_lines, sprintf("  %-10s baseline log2FC=%+.4f padj=%.4g -- %s",
                                         g, r$log2FC, r$padj, ifelse(sig, "*** SIGNIFICANT AT BASELINE ***", "n.s.")))
}
if (any_sig_baseline) {
  check_lines <- c(check_lines, "",
    "*** IMPORTANT: NOT all 4 genes are non-significant at baseline. Loxl1 is significantly",
    "higher in ByJcl already at baseline (log2FC=+0.68, padj=0.0052), though at a much smaller",
    "effect size than at Day 5 (log2FC=+1.37). The 'no baseline difference' framing holds cleanly",
    "for Serpine1, Col4a1, and Col4a2 (all padj>0.45 at baseline) but must be qualified for Loxl1:",
    "it shows a modest, already-significant baseline elevation that further increases by Day 5,",
    "not a purely Day-5-specific onset. The figure and caption below reflect this (Loxl1's baseline",
    "panel is annotated with its actual padj, not simply 'n.s.')."
  )
}
writeLines(check_lines, file.path(log_dir, "29_Figure7_value_check.txt"))
cat(paste(check_lines, collapse = "\n"), "\n\n")

# ---------------------------------------------------------------------------
# figure -- plain bar graph, mean + SEM, A-ADR1 dropped entirely
# ---------------------------------------------------------------------------
plot_df <- long %>% filter(!is_A1)   # A-ADR1 excluded from the figure altogether

sig_star <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "n.s."))))
}

dodge_off <- 0.75 / 4  # position_dodge(width = 0.75), 2 groups per x -> bar centers at +/- width/4

summary_df <- plot_df %>%
  group_by(gene, timepoint, substrain) %>%
  summarise(mean_cpm = mean(cpm), sem = sd(cpm) / sqrt(n()), n = n(), .groups = "drop")

bracket_df <- summary_df %>%
  group_by(gene, timepoint) %>%
  summarise(y0 = max(mean_cpm + sem), .groups = "drop") %>%
  mutate(x_num = ifelse(timepoint == "Baseline", 1, 2),
         xmin = x_num - dodge_off, xmax = x_num + dodge_off,
         y = y0 * 1.10)
bracket_df$padj <- ifelse(bracket_df$timepoint == "Baseline",
                           sapply(bracket_df$gene, function(g) baseline_de$padj[baseline_de$gene == as.character(g)]),
                           sapply(bracket_df$gene, function(g) d5_de$padj[d5_de$gene == as.character(g)]))
bracket_df$star <- sig_star(bracket_df$padj)
# "n.s." is wider/taller than "*"/"**"/"***" and sits right above the
# bracket line by default (vjust=0.5 straddles the line) -- give it its own,
# larger vjust so the label clears the line instead of overlapping it, and
# a slightly smaller font so 3-4 characters don't crowd the narrow bracket span
bracket_df$text_vjust <- ifelse(bracket_df$star == "n.s.", -0.7, -0.35)
bracket_df$text_size <- ifelse(bracket_df$star == "n.s.", pt_mm(10), pt_mm(13))

y_max_by_gene <- bracket_df %>% group_by(gene) %>% summarise(ymax_all = max(y) * 1.18, .groups = "drop")

p <- ggplot(plot_df, aes(x = timepoint, y = cpm, fill = substrain)) +
  stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.75),
               width = 0.65, color = "black", linewidth = 0.45) +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(width = 0.75),
               width = 0.22, linewidth = 0.6) +
  scale_fill_manual(values = c("AJcl" = "#8FB8DE", "ByJcl" = "#0072B2"), name = "Substrain") +
  facet_wrap(~gene, scales = "free_y", ncol = 2) +
  geom_segment(data = bracket_df, aes(x = xmin, xend = xmax, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.5) +
  geom_segment(data = bracket_df, aes(x = xmin, xend = xmin, y = y * 0.97, yend = y), inherit.aes = FALSE, linewidth = 0.5) +
  geom_segment(data = bracket_df, aes(x = xmax, xend = xmax, y = y * 0.97, yend = y), inherit.aes = FALSE, linewidth = 0.5) +
  geom_text(data = bracket_df, aes(x = x_num, y = y, label = star, vjust = text_vjust, size = text_size),
            inherit.aes = FALSE, fontface = "bold") +
  scale_size_identity() +
  geom_blank(data = y_max_by_gene, aes(x = 1, y = ymax_all), inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) +
  labs(
    title = "ADR-induced upregulation of ECM/adhesion transcripts in ByJcl glomeruli",
    subtitle = "CPM (mean + SEM), AJcl vs ByJcl, baseline and Day 5 post-ADR\n(A-ADR1 excluded; n=3 per group except AJcl Day 5, n=2)",
    x = NULL, y = "CPM"
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(size = 15.5, face = "bold"),
    plot.subtitle = element_text(size = 10.5),
    axis.title.y = element_text(size = 14.5, face = "bold"),
    axis.text.x = element_text(size = 12.5),
    axis.text.y = element_text(size = 11.5),
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12.5),
    legend.text = element_text(size = 12),
    legend.position = "bottom"
  )

for (ext in c("pdf", "png")) {
  fp <- file.path(fig_dir, paste0("Figure7_ECM_genes_CPM.", ext))
  ggsave(fp, plot = p, width = 9, height = 7.5, device = if (ext == "pdf") cairo_pdf else ext, dpi = 300)
  cat("Saved:", fp, "\n")
}

caption <- paste0(
  "ADR-induced upregulation of ECM/adhesion-related transcripts in ByJcl glomeruli. CPM (mean + SEM) ",
  "of Serpine1, Loxl1, Col4a1, and Col4a2 in AJcl and ByJcl glomeruli at baseline and Day 5 after ADR ",
  "(A-ADR1 excluded for reduced glomerular purity; see FigS_glomerular_purity_QC.pdf; n=3 per group ",
  "except AJcl Day 5, n=2). All four transcripts were significantly higher in ByJcl at Day 5 ",
  "(Serpine1 FDR = 1.7x10-9; Loxl1 FDR = 4.9x10-12; Col4a1 FDR = 4.5x10-3; Col4a2 FDR = 1.5x10-3; DESeq2, ",
  "ByJcl vs AJcl). At baseline, Serpine1, Col4a1, and Col4a2 showed no significant difference ",
  "(all FDR>0.45); Loxl1 showed a small but statistically significant baseline elevation in ByJcl ",
  "(log2FC=+0.68, FDR=0.0052) that increased further by Day 5 (log2FC=+1.37), and so is not purely a ",
  "Day-5-specific onset. Brackets: * FDR<0.05, ** FDR<0.01, *** FDR<0.001, n.s. not significant ",
  "(DESeq2 Wald test, ByJcl vs AJcl within timepoint)."
)
writeLines(caption, file.path(log_dir, "29_Figure7_caption_draft.txt"))
cat("\n=== Figure caption (draft) ===\n", caption, "\n")

cat("\n=== 29_Figure7_ECM_genes_CPM.R complete ===\n")
