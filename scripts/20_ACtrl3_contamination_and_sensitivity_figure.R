#!/usr/bin/env Rscript
# Control-side companion to 18_A1_contamination_QC_figure.R and
# 17_A1_sensitivity_supplementary_figure.R: (A) tubular-contamination QC panel
# for the 6 Ctrl (baseline) samples, highlighting A-Ctrl3, and (B) NES
# dumbbell comparison for baseline_B_vs_A with A-Ctrl3 included (main) vs
# excluded (sensitivity), using the judgment table from
# 19_ACtrl3_sensitivity_dds_and_DE.R.

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
})

indir <- "/usr/local/jupyter/ADR_BALB"
outdir <- file.path(indir, "outputs")

# ---------------------------------------------------------------------------
# Part 1: contamination QC figure, Ctrl samples, A-Ctrl3 highlighted
# ---------------------------------------------------------------------------
merged <- read.delim(file.path(outdir, "tables/00_merged_counts.tsv"),
                      check.names = FALSE, stringsAsFactors = FALSE)
rownames(merged) <- merged$Geneid

ctrl_samples <- c("A-Ctrl1", "A-Ctrl2", "A-Ctrl3", "B-Ctrl1", "B-Ctrl2", "B-Ctrl3")
all_samples <- c(ctrl_samples, "A-ADR1","A-ADR2","A-ADR3","B-ADR1","B-ADR2","B-ADR3")
counts_only <- merged[, all_samples]
cpm <- sweep(counts_only, 2, colSums(counts_only), FUN = "/") * 1e6

podo_markers   <- c("Nphs1", "Nphs2", "Wt1", "Podxl", "Synpo", "Nes")
tubule_markers <- c("Lrp2", "Slc34a1", "Aqp1", "Slc12a1", "Umod", "Aqp2")

tub_df <- as.data.frame(t(cpm[tubule_markers, ctrl_samples]))
tub_df$sample <- rownames(tub_df)
tub_long <- tub_df %>% pivot_longer(-sample, names_to = "gene", values_to = "cpm")
tub_long$sample <- factor(tub_long$sample, levels = ctrl_samples)
tub_long$is_flagged <- tub_long$sample == "A-Ctrl3"
tub_long$gene <- factor(tub_long$gene, levels = tubule_markers)

panelA <- ggplot(tub_long, aes(x = sample, y = cpm, fill = is_flagged)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  facet_wrap(~gene, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#0072B2"), guide = "none") +
  labs(title = "A. Renal-tubular-epithelial marker CPM, per Ctrl (baseline) sample",
       subtitle = "A-Ctrl3 (vermillion) is elevated on tubular markers, comparable in magnitude to A-ADR1",
       x = NULL, y = "CPM") +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        strip.text = element_text(face = "bold", size = 10),
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9))

podo_sum_ctrl <- colSums(cpm[podo_markers, ctrl_samples])
tubule_sum_ctrl <- colSums(cpm[tubule_markers, ctrl_samples])
ratio_ctrl <- podo_sum_ctrl / tubule_sum_ctrl
ratio_df <- data.frame(sample = factor(names(ratio_ctrl), levels = ctrl_samples),
                        ratio = as.numeric(ratio_ctrl))
ratio_df$is_flagged <- ratio_df$sample == "A-Ctrl3"

next_lowest <- sort(ratio_df$ratio[!ratio_df$is_flagged])[1]
flagged_ratio <- ratio_df$ratio[ratio_df$is_flagged]

panelB <- ggplot(ratio_df, aes(x = sample, y = ratio, fill = is_flagged)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", ratio)), vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#0072B2"), guide = "none") +
  scale_y_log10(expand = expansion(mult = c(0.02, 0.18))) +
  labs(title = "B. Podocyte:tubular marker-sum ratio per Ctrl sample (log scale)",
       subtitle = sprintf("A-Ctrl3 ratio = %.2f, the lowest of all 6 Ctrl samples (next-lowest = %.2f, %.1fx higher)",
                           flagged_ratio, next_lowest, next_lowest / flagged_ratio),
       x = NULL, y = "Podocyte:tubular\nmarker-sum ratio") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9),
        axis.title.y = element_text(size = 11))

p_qc <- arrangeGrob(panelA, panelB, ncol = 1, heights = c(1.3, 1))
for (ext in c("png", "pdf")) {
  fp <- file.path(outdir, paste0("figures/FigS_ACtrl3_contamination_QC.", ext))
  ggsave(fp, plot = p_qc, width = 9, height = 10.5, dpi = 300)
  cat("Saved:", fp, "\n")
}

# ---------------------------------------------------------------------------
# Part 2: NES dumbbell, baseline_B_vs_A, A-Ctrl3 included vs excluded
# ---------------------------------------------------------------------------
judg <- read.delim(file.path(outdir, "tables/GSEA_ACtrl3_sensitivity_judgment_table.tsv"),
                    stringsAsFactors = FALSE)

pathway_labels <- c(
  TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING   = "Podocyte ageing (Fig. 6B set)",
  REACTOME_INTEGRIN_SIGNALING                 = "Integrin signaling (Fig. 6D set)",
  REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS = "Integrin cell surface interactions",
  REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION  = "ECM organization (Fig. 6C theme)",
  KARAISKOS2018_PODOCYTE_EXCLUSIVE            = "Karaiskos podocyte, exclusive (n=12)",
  KARAISKOS2018_PODOCYTE_TOP50                = "Karaiskos podocyte, top50 (n=49)"
)
pathway_order <- rev(names(pathway_labels))

dat <- judg %>%
  mutate(
    status_label = ifelse(actrl3_status == "ACtrl3_included_main",
                           "A-Ctrl3 included (main analysis)", "A-Ctrl3 excluded (sensitivity)"),
    pathway_label = factor(pathway_labels[pathway], levels = pathway_labels[pathway_order])
  )

wide <- dat %>%
  select(pathway, pathway_label, status_label, NES, FDR_joint) %>%
  pivot_wider(names_from = status_label, values_from = c(NES, FDR_joint))
names(wide) <- gsub("A-Ctrl3 included \\(main analysis\\)", "included", names(wide))
names(wide) <- gsub("A-Ctrl3 excluded \\(sensitivity\\)", "excluded", names(wide))
wide$sign_flip <- sign(wide$NES_included) != sign(wide$NES_excluded)
wide$flip_label <- ifelse(wide$sign_flip, "Sign flip", "Robust (same sign both ways)")

dat_long <- dat %>%
  select(pathway_label, status_label, NES) %>%
  left_join(unique(wide[, c("pathway_label", "flip_label")]), by = "pathway_label")
dat_long$point_lab <- sprintf("%.2f", dat_long$NES)

p_nes <- ggplot(dat_long, aes(x = NES, y = pathway_label)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey85") +
  geom_line(aes(group = pathway_label, color = flip_label), linewidth = 1.1) +
  geom_point(aes(shape = status_label, fill = status_label), size = 4, color = "black", stroke = 0.6) +
  ggrepel::geom_text_repel(aes(label = point_lab, color = status_label),
            size = 3.6, fontface = "bold", show.legend = FALSE, seed = 42,
            box.padding = 0.35, point.padding = 0.25, min.segment.length = 0.4,
            max.overlaps = Inf, direction = "both", force = 3) +
  scale_shape_manual(values = c("A-Ctrl3 included (main analysis)" = 21, "A-Ctrl3 excluded (sensitivity)" = 24)) +
  scale_fill_manual(values = c("A-Ctrl3 included (main analysis)" = "#0072B2", "A-Ctrl3 excluded (sensitivity)" = "#E69F00")) +
  scale_color_manual(
    values = c("Robust (same sign both ways)" = "grey60", "Sign flip" = "#D55E00",
               "A-Ctrl3 included (main analysis)" = "#0072B2", "A-Ctrl3 excluded (sensitivity)" = "#E69F00"),
    breaks = c("Robust (same sign both ways)", "Sign flip")
  ) +
  labs(
    title = "A-Ctrl3 sensitivity check: baseline_B_vs_A gene-set enrichment (NES),\nA-Ctrl3 included vs. excluded",
    subtitle = "Canonical stat-ranked fgsea. All 6 focal gene sets keep the same sign with or without A-Ctrl3.",
    x = "Normalized Enrichment Score (NES)", y = NULL,
    shape = "A-Ctrl3 status", fill = "A-Ctrl3 status", color = "Robustness"
  ) +
  guides(fill = guide_legend(override.aes = list(color = "black")),
         shape = guide_legend(override.aes = list(color = "black"))) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(size = 12.5, face = "bold"),
        plot.subtitle = element_text(size = 9.5),
        axis.title = element_text(size = 13, face = "bold"),
        axis.text = element_text(size = 10.5),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10.5),
        legend.position = "bottom", legend.box = "vertical")

for (ext in c("png", "pdf")) {
  fp <- file.path(outdir, paste0("figures/FigS_ACtrl3_sensitivity_NES_comparison.", ext))
  ggsave(fp, plot = p_nes, width = 9.5, height = 7, dpi = 300)
  cat("Saved:", fp, "\n")
}

summary_lines <- c(
  "=== A-Ctrl3 sensitivity check summary ===",
  "",
  "QC (tubular-marker CPM, Ctrl samples):",
  sprintf("  Podocyte:tubular ratio, A-Ctrl3 = %.2f (lowest of 6 Ctrl samples, next-lowest = %.2f, %.1fx higher)",
          flagged_ratio, next_lowest, next_lowest / flagged_ratio),
  "",
  "baseline_B_vs_A GSEA NES, A-Ctrl3 included (main) vs excluded (sensitivity):",
  sprintf("  %-38s NES_included=%+.2f (FDR=%.1e)  NES_excluded=%+.2f (FDR=%.1e)  %s",
          as.character(wide$pathway_label), wide$NES_included, wide$FDR_joint_included,
          wide$NES_excluded, wide$FDR_joint_excluded, wide$flip_label),
  "",
  sprintf("Sign flips: %d/6 focal gene sets. All 6 gene sets retain the same direction with A-Ctrl3 excluded.",
          sum(wide$sign_flip)),
  "Genome-wide baseline_B_vs_A log2FC concordance (see logs/19_ACtrl3_sensitivity_concordance.txt): Pearson r=0.88, 88.1% direction-concordant genome-wide, 100% among genes significant in either configuration.",
  "Conclusion: A-Ctrl3 shows tubular contamination comparable in magnitude to A-ADR1, but unlike A1's effect on the podocyte-ageing gene set in ADR_B_vs_A, excluding A-Ctrl3 does not flip the sign of any focal gene set for baseline_B_vs_A -- the baseline substrain-difference conclusions in this re-analysis are robust to A-Ctrl3's inclusion/exclusion."
)
writeLines(summary_lines, file.path(outdir, "logs/20_ACtrl3_sensitivity_figure_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

cat("\n=== 20_ACtrl3_contamination_and_sensitivity_figure.R complete ===\n")
