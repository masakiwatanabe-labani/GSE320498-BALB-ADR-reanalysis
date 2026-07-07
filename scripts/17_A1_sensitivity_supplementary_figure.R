#!/usr/bin/env Rscript
# R1-7 / R2-4: supplementary figure for the A1 (low-purity A-ADR1 sample)
# sensitivity analysis. The manuscript states "sensitivity analysis including
# this sample confirmed the robustness of the pathway-level findings (data
# not shown)". This script makes that comparison visible: NES with A1
# excluded (main analysis) vs included (sensitivity), for every focal gene
# set, in every comparison where A1 status is applicable (ADR_B_vs_A, the
# Fig. 6 comparison, and A_ADR_vs_Ctrl). Data source is the existing
# 09_gsea_stat_ranking_canonical.R output
# (tables/GSEA_canonical_focal_judgment_table.tsv) -- no new GSEA run.

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

indir <- "/usr/local/jupyter/ADR_BALB"
outdir <- file.path(indir, "outputs")

judg <- read.delim(file.path(outdir, "tables/GSEA_canonical_focal_judgment_table.tsv"),
                    stringsAsFactors = FALSE)
judg <- judg[!is.na(judg$comparison) & judg$comparison != "", ]

pathway_labels <- c(
  TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING   = "Podocyte ageing (Fig. 6B set)",
  REACTOME_INTEGRIN_SIGNALING                 = "Integrin signaling (Fig. 6D set)",
  REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS = "Integrin cell surface interactions",
  REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION  = "ECM organization (Fig. 6C theme)",
  KARAISKOS2018_PODOCYTE_EXCLUSIVE            = "Karaiskos podocyte, exclusive (n=12)",
  KARAISKOS2018_PODOCYTE_TOP50                = "Karaiskos podocyte, top50 (n=49)"
)
pathway_order <- rev(names(pathway_labels))

comparison_labels <- c(
  ADR_B_vs_A    = "ADR_B_vs_A: Day 5, ByJcl vs AJcl (Fig. 6 comparison)",
  A_ADR_vs_Ctrl = "A_ADR_vs_Ctrl: ADR response within AJcl"
)

dat <- judg %>%
  filter(comparison %in% names(comparison_labels), pathway %in% names(pathway_labels)) %>%
  mutate(
    a1_status_label = ifelse(a1_status == "A1_excluded_main",
                              "A1-excluded (main analysis)", "A1-included (sensitivity)"),
    pathway_label = factor(pathway_labels[pathway], levels = pathway_labels[pathway_order]),
    comparison_label = factor(comparison_labels[comparison], levels = comparison_labels)
  )

wide <- dat %>%
  select(comparison, comparison_label, pathway, pathway_label, a1_status_label, NES, FDR_joint) %>%
  pivot_wider(names_from = a1_status_label, values_from = c(NES, FDR_joint))

names(wide) <- gsub("A1-excluded \\(main analysis\\)", "excluded", names(wide))
names(wide) <- gsub("A1-included \\(sensitivity\\)", "included", names(wide))

wide$sign_flip <- sign(wide$NES_excluded) != sign(wide$NES_included)
wide$sig_either <- (!is.na(wide$FDR_joint_excluded) & wide$FDR_joint_excluded < 0.05) |
                    (!is.na(wide$FDR_joint_included) & wide$FDR_joint_included < 0.05)
wide$flip_label <- ifelse(
  !wide$sign_flip, "Robust to A1 status",
  ifelse(wide$sig_either, "Sign flip, significant in >=1 configuration (of concern)",
         "Sign flip, not significant either way (noise-level)")
)

dat_long <- dat %>%
  select(comparison_label, pathway_label, a1_status_label, NES, FDR_joint) %>%
  left_join(unique(wide[, c("comparison_label", "pathway_label", "sign_flip", "flip_label")]),
             by = c("comparison_label", "pathway_label"))

fmt_fdr <- function(f) ifelse(is.na(f), "NA", formatC(f, format = "e", digits = 1))
dat_long$point_lab <- sprintf("%.2f", dat_long$NES)

flip_colors <- c(
  "Robust (same sign both ways)" = "grey60",
  "Sign flip -- noise-level (n.s. both ways)" = "grey60",
  "Sign flip -- of concern (sig. in >=1 way)" = "#D55E00"
)
flip_linetypes <- c(
  "Robust (same sign both ways)" = "solid",
  "Sign flip -- noise-level (n.s. both ways)" = "dotted",
  "Sign flip -- of concern (sig. in >=1 way)" = "solid"
)
flip_linewidths <- c(
  "Robust (same sign both ways)" = 1.1,
  "Sign flip -- noise-level (n.s. both ways)" = 0.7,
  "Sign flip -- of concern (sig. in >=1 way)" = 1.8
)
relabel_flip <- c(
  "Robust to A1 status" = "Robust (same sign both ways)",
  "Sign flip, not significant either way (noise-level)" = "Sign flip -- noise-level (n.s. both ways)",
  "Sign flip, significant in >=1 configuration (of concern)" = "Sign flip -- of concern (sig. in >=1 way)"
)
dat_long$flip_label <- factor(relabel_flip[as.character(dat_long$flip_label)], levels = names(flip_colors))

p <- ggplot(dat_long, aes(x = NES, y = pathway_label)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey85") +
  geom_line(aes(group = pathway_label, color = flip_label, linetype = flip_label, linewidth = flip_label)) +
  geom_point(aes(shape = a1_status_label, fill = a1_status_label), size = 4, color = "black", stroke = 0.6) +
  ggrepel::geom_text_repel(aes(label = point_lab, color = a1_status_label),
            size = 3.4, fontface = "bold", show.legend = FALSE, seed = 42,
            box.padding = 0.35, point.padding = 0.25, min.segment.length = 0.4,
            max.overlaps = Inf, direction = "both", force = 3) +
  scale_shape_manual(values = c("A1-excluded (main analysis)" = 21, "A1-included (sensitivity)" = 24)) +
  scale_fill_manual(values = c("A1-excluded (main analysis)" = "#0072B2", "A1-included (sensitivity)" = "#E69F00")) +
  scale_linetype_manual(values = flip_linetypes, guide = "none") +
  scale_linewidth_manual(values = flip_linewidths, guide = "none") +
  scale_color_manual(
    values = c(flip_colors, "A1-excluded (main analysis)" = "#0072B2", "A1-included (sensitivity)" = "#E69F00"),
    breaks = names(flip_colors)
  ) +
  facet_wrap(~comparison_label, ncol = 1, scales = "free_y") +
  labs(
    title = "A1 sensitivity check: gene-set enrichment (NES), flagged low-purity sample\n(A-ADR1) excluded vs. included",
    subtitle = paste0(
      "Canonical stat-ranked fgsea (joint FDR vs. full 1293-1294 set universe).\n",
      "Positive NES = higher in ByJcl (ADR_B_vs_A) or higher after ADR (A_ADR_vs_Ctrl), per project sign convention."
    ),
    x = "Normalized Enrichment Score (NES)", y = NULL,
    shape = "A1 status", fill = "A1 status", color = "Robustness to A1 status"
  ) +
  guides(fill = guide_legend(override.aes = list(color = "black")),
         shape = guide_legend(override.aes = list(color = "black")),
         color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(linewidth = 1.3))) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 11),
    strip.text = element_text(size = 10.5, face = "bold"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.position = "bottom",
    legend.box = "vertical"
  )

for (ext in c("png", "pdf")) {
  fp <- file.path(outdir, paste0("figures/FigS_A1_sensitivity_NES_comparison.", ext))
  ggsave(fp, plot = p, width = 12, height = 10, dpi = 300)
  cat("Saved:", fp, "\n")
}

# text summary for the response letter / caption
summary_lines <- c(
  "=== A1 sensitivity check (NES, A1-excluded vs A1-included) ===",
  ""
)
for (cl in levels(dat_long$comparison_label)) {
  summary_lines <- c(summary_lines, sprintf("-- %s --", cl))
  sub <- wide[wide$comparison_label == cl, ]
  sub <- sub[order(match(sub$pathway_label, levels(dat_long$pathway_label))), ]
  for (i in seq_len(nrow(sub))) {
    r <- sub[i, ]
    summary_lines <- c(summary_lines, sprintf(
      "  %-38s NES_excluded=%+.2f (FDR=%s)  NES_included=%+.2f (FDR=%s)  %s",
      as.character(r$pathway_label), r$NES_excluded, fmt_fdr(r$FDR_joint_excluded),
      r$NES_included, fmt_fdr(r$FDR_joint_included), r$flip_label
    ))
  }
  summary_lines <- c(summary_lines, "")
}
n_flip <- sum(wide$sign_flip)
n_total <- nrow(wide)
summary_lines <- c(summary_lines, sprintf(
  "Overall: %d/%d gene-set x comparison pairs tested for A1 sensitivity; %d sign flip(s) with A1 status (%s).",
  n_total, n_total, n_flip,
  paste(sprintf("%s / %s", wide$pathway_label[wide$sign_flip], wide$comparison_label[wide$sign_flip]), collapse = "; ")
))
writeLines(summary_lines, file.path(outdir, "logs/17_A1_sensitivity_figure_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

cat("\n=== 17_A1_sensitivity_supplementary_figure.R complete ===\n")
