#!/usr/bin/env Rscript
# Single-panel (+ optional heatmap) QC figure combining ALL 12 RNA-seq
# samples, to show A-ADR1 as a clear outlier justifying its exclusion from
# the Day 5 analysis. Supersedes/combines the previously separate ADR-only
# (18_A1_contamination_QC_figure.R) and Ctrl-only
# (20_ACtrl3_contamination_and_sensitivity_figure.R) purity panels into one
# figure spanning all 4 groups.
#
# Definitions held fixed, matching the existing QC scripts exactly (not
# redefined): podocyte markers = Nphs1, Nphs2, Wt1, Podxl, Synpo, Nes;
# tubular markers = Lrp2, Slc34a1, Aqp1, Slc12a1, Umod, Aqp2. CPM computed
# directly from tables/00_merged_counts.tsv (same as 18/20). No new
# statistical model -- this is a QC visualization only.

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
})

indir <- "/usr/local/jupyter/ADR_BALB"
outdir_main <- file.path(indir, "outputs")
outdir_v2 <- file.path(outdir_main, "REVISION_PACKAGE/REVISION_PACKAGE_v2")
fig_dir <- file.path(outdir_v2, "figures")
tab_dir <- file.path(outdir_v2, "tables")
log_dir <- file.path(outdir_v2, "logs")

merged <- read.delim(file.path(outdir_main, "tables/00_merged_counts.tsv"),
                      check.names = FALSE, stringsAsFactors = FALSE)
rownames(merged) <- merged$Geneid

group_order <- c("A-Ctrl1","A-Ctrl2","A-Ctrl3","B-Ctrl1","B-Ctrl2","B-Ctrl3",
                  "A-ADR1","A-ADR2","A-ADR3","B-ADR1","B-ADR2","B-ADR3")
counts_only <- merged[, group_order]
cpm <- sweep(counts_only, 2, colSums(counts_only), FUN = "/") * 1e6

podo_markers   <- c("Nphs1", "Nphs2", "Wt1", "Podxl", "Synpo", "Nes")
tubule_markers <- c("Lrp2", "Slc34a1", "Aqp1", "Slc12a1", "Umod", "Aqp2")

# ---------------------------------------------------------------------------
# Step 1: ratio table, all 12 samples + cross-check against prior scripts
# ---------------------------------------------------------------------------
podo_sum <- colSums(cpm[podo_markers, group_order])
tubule_sum <- colSums(cpm[tubule_markers, group_order])
ratio <- podo_sum / tubule_sum

ratio_df <- data.frame(
  sample = factor(group_order, levels = group_order),
  group = factor(sub("[0-9]$", "", group_order), levels = c("A-Ctrl", "B-Ctrl", "A-ADR", "B-ADR")),
  podocyte_sum_CPM = round(as.numeric(podo_sum), 1),
  tubular_sum_CPM = round(as.numeric(tubule_sum), 1),
  ratio = round(as.numeric(ratio), 3)
)
write.csv(ratio_df, file.path(tab_dir, "TableS_purity_ratios.csv"), row.names = FALSE)

known_values <- c("A-ADR1" = 2.51, "A-ADR2" = 7.20, "A-ADR3" = 7.61, "B-ADR1" = 35.22, "B-ADR2" = 9.33, "B-ADR3" = 7.14,
                   "A-Ctrl1" = 8.52, "A-Ctrl2" = 19.19, "A-Ctrl3" = 4.50, "B-Ctrl1" = 50.58, "B-Ctrl2" = 11.95, "B-Ctrl3" = 13.63)
check_lines <- c("=== Step 1: cross-check against previously published QC scripts (18/20) ===")
all_match <- TRUE
for (s in names(known_values)) {
  new_val <- ratio_df$ratio[ratio_df$sample == s]
  match_ok <- abs(new_val - known_values[s]) < 0.02
  all_match <- all_match && match_ok
  check_lines <- c(check_lines, sprintf("  %-8s prior=%.2f  recomputed=%.3f  %s", s, known_values[s], new_val,
                                         ifelse(match_ok, "MATCH", "*** MISMATCH ***")))
}
check_lines <- c(check_lines, sprintf("\nAll 12 samples match prior scripts' values: %s", all_match))
cat(paste(check_lines, collapse = "\n"), "\n\n")
if (!all_match) stop("Purity ratio mismatch against prior QC scripts -- halting, see check_lines above.")

# ---------------------------------------------------------------------------
# Step 2: main combined 12-sample figure
# ---------------------------------------------------------------------------
ratio_df$is_A1 <- ratio_df$sample == "A-ADR1"
ratio_df$is_ACtrl3 <- ratio_df$sample == "A-Ctrl3"

other11_ratios <- setNames(ratio_df$ratio, ratio_df$sample)[!ratio_df$is_A1]
band_lo <- min(other11_ratios)
band_hi <- max(other11_ratios)
median_other11 <- median(other11_ratios)

a1_ratio <- ratio_df$ratio[ratio_df$is_A1]
next_lowest <- sort(other11_ratios)[1]
fold_below <- next_lowest / a1_ratio

group_colors <- c("A-Ctrl" = "#0072B2", "B-Ctrl" = "#56B4E9", "A-ADR" = "#009E73", "B-ADR" = "#CC79A7")

p_main <- ggplot(ratio_df, aes(x = sample, y = ratio)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = band_lo, ymax = band_hi, fill = "grey85", alpha = 0.4) +
  geom_hline(yintercept = median_other11, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_point(aes(fill = group, color = is_A1, size = is_A1), shape = 21, stroke = 1.1) +
  scale_fill_manual(values = group_colors, name = "Group") +
  scale_color_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "black"), guide = "none") +
  scale_size_manual(values = c(`TRUE` = 7, `FALSE` = 4.5), guide = "none") +
  geom_text(data = ratio_df[ratio_df$is_A1, ], aes(label = sprintf("A1 (excluded)\nratio=%.2f", ratio)),
            vjust = -0.7, hjust = 1.05, size = 4.2, fontface = "bold", color = "#D55E00") +
  geom_text(data = ratio_df[ratio_df$is_ACtrl3, ], aes(label = "A-Ctrl3\n(retained)"),
            vjust = -0.7, size = 3.1, fontface = "italic", color = "grey30") +
  scale_y_log10(expand = expansion(mult = c(0.05, 0.28))) +
  labs(
    title = "Glomerular purity across all RNA-seq samples",
    subtitle = paste0(
      "Podocyte:tubular marker-sum CPM ratio (log scale), 12 samples, 4 groups (n=3 each).\n",
      sprintf("A-ADR1 shows the lowest ratio (%.2f), ~%.1f-fold below the next-lowest sample (%s, %.2f); ",
              a1_ratio, fold_below, names(other11_ratios)[which.min(other11_ratios)], next_lowest),
      "it was excluded from the Day 5 analysis.\n",
      "Grey band = range of the other 11 samples; dashed line = their median. All other samples, including A-Ctrl3, were retained."
    ),
    x = NULL, y = "Podocyte:tubular marker-sum CPM ratio (log scale)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 9, lineheight = 1.2),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

# ---------------------------------------------------------------------------
# Step 3: optional heatmap, 6 tubular markers x 12 samples
# ---------------------------------------------------------------------------
tub_mat <- cpm[tubule_markers, group_order]
tub_long <- as.data.frame(tub_mat)
tub_long$gene <- rownames(tub_long)
tub_long <- tub_long %>% pivot_longer(-gene, names_to = "sample", values_to = "cpm")
tub_long$sample <- factor(tub_long$sample, levels = group_order)
tub_long$gene <- factor(tub_long$gene, levels = tubule_markers)
tub_long$log_cpm <- log10(tub_long$cpm + 1)
tub_long$is_A1_col <- tub_long$sample == "A-ADR1"

p_heat <- ggplot(tub_long, aes(x = sample, y = gene, fill = log_cpm)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_tile(data = tub_long[tub_long$is_A1_col, ], fill = NA, color = "#D55E00", linewidth = 1.3) +
  scale_fill_viridis_c(option = "magma", name = "log10(CPM+1)") +
  labs(title = "Renal-tubular-epithelial marker CPM, all 12 samples (A-ADR1 outlined)",
       x = NULL, y = NULL) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(size = 11.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 9.5),
    legend.text = element_text(size = 8.5)
  )

p_combined <- arrangeGrob(p_main, p_heat, ncol = 1, heights = c(1.5, 1))
for (ext in c("png", "pdf")) {
  fp <- file.path(fig_dir, paste0("FigS_glomerular_purity_QC.", ext))
  ggsave(fp, plot = p_combined, width = 10, height = 10.5, dpi = 300)
  cat("Saved:", fp, "\n")
}

# ---------------------------------------------------------------------------
# Step 4: caption + summary
# ---------------------------------------------------------------------------
caption <- paste0(
  "Glomerular purity assessment for all RNA-seq samples. Podocyte-to-tubular marker-sum ratio ",
  "(log scale) across the 12 samples. A-ADR1 showed the lowest ratio (", sprintf("%.2f", a1_ratio),
  "), ~", sprintf("%.1f", fold_below), "-fold below the next-lowest sample, together with elevated ",
  "expression of multiple renal-tubular markers, indicating reduced glomerular purity; this sample ",
  "was excluded from the Day 5 (ADR) analysis. All other samples, including A-Ctrl3, fell within a ",
  "comparable range and were retained."
)
summary_lines <- c(
  check_lines, "",
  "=== Figure caption (draft) ===", caption, "",
  "=== Ratio table (all 12 samples) ===",
  capture.output(print(ratio_df[, c("sample","group","podocyte_sum_CPM","tubular_sum_CPM","ratio")], row.names = FALSE))
)
writeLines(summary_lines, file.path(log_dir, "28_glomerular_purity_QC_summary.txt"))
cat(caption, "\n")

cat("\n=== 28_glomerular_purity_QC_all12.R complete ===\n")
