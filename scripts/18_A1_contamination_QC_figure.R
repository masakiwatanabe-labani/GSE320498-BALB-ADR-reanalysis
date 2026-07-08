#!/usr/bin/env Rscript
# Supplementary QC figure justifying the A-ADR1 exclusion (main-analysis
# sensitivity choice): shows that A-ADR1 has elevated non-glomerular (renal
# tubular epithelial) marker expression relative to the other 5 ADR samples,
# i.e. more tubular contamination in the glomeruli-enriched preparation --
# not a reduced podocyte signal, but an excess tubular one.
#
# Recomputes CPM directly from tables/00_merged_counts.tsv using the same
# marker sets and logic as 00_load_merge_qc.R (no new sequencing analysis).

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

indir <- "/usr/local/jupyter/ADR_BALB"
outdir <- file.path(indir, "outputs")

merged <- read.delim(file.path(outdir, "tables/00_merged_counts.tsv"),
                      check.names = FALSE, stringsAsFactors = FALSE)
rownames(merged) <- merged$Geneid

adr_samples <- c("A-ADR1", "A-ADR2", "A-ADR3", "B-ADR1", "B-ADR2", "B-ADR3")
all_samples <- c("A-Ctrl1","A-Ctrl2","A-Ctrl3","B-Ctrl1","B-Ctrl2","B-Ctrl3", adr_samples)
counts_only <- merged[, all_samples]
cpm <- sweep(counts_only, 2, colSums(counts_only), FUN = "/") * 1e6

podo_markers   <- c("Nphs1", "Nphs2", "Wt1", "Podxl", "Synpo", "Nes")
tubule_markers <- c("Lrp2", "Slc34a1", "Aqp1", "Slc12a1", "Umod", "Aqp2")

flag_color <- c("A-ADR1" = "#D55E00")
other_color <- "#0072B2"
sample_colors <- setNames(rep(other_color, length(adr_samples)), adr_samples)
sample_colors["A-ADR1"] <- flag_color

# ---------------------------------------------------------------------------
# Panel A: per-marker tubular-gene CPM, all 6 ADR samples, small multiples
# ---------------------------------------------------------------------------
tub_df <- as.data.frame(t(cpm[tubule_markers, adr_samples]))
tub_df$sample <- rownames(tub_df)
tub_long <- tub_df %>% pivot_longer(-sample, names_to = "gene", values_to = "cpm")
tub_long$sample <- factor(tub_long$sample, levels = adr_samples)
tub_long$is_A1 <- tub_long$sample == "A-ADR1"
tub_long$gene <- factor(tub_long$gene, levels = tubule_markers)

panelA <- ggplot(tub_long, aes(x = sample, y = cpm, fill = is_A1)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  facet_wrap(~gene, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#0072B2"), guide = "none") +
  labs(title = "A. Renal-tubular-epithelial marker CPM, per ADR sample",
       subtitle = "A-ADR1 (vermillion) is elevated on 5 of 6 tubular markers -- consistent with reduced glomerular purity",
       x = NULL, y = "CPM") +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        strip.text = element_text(face = "bold", size = 10),
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9))

# ---------------------------------------------------------------------------
# Panel B: podocyte:tubular marker-sum ratio per ADR sample (log scale)
# ---------------------------------------------------------------------------
podo_sum <- colSums(cpm[podo_markers, adr_samples])
tubule_sum <- colSums(cpm[tubule_markers, adr_samples])
ratio <- podo_sum / tubule_sum

ratio_df <- data.frame(sample = factor(names(ratio), levels = adr_samples),
                        podo_sum = as.numeric(podo_sum),
                        tubule_sum = as.numeric(tubule_sum),
                        ratio = as.numeric(ratio))
ratio_df$is_A1 <- ratio_df$sample == "A-ADR1"

panelB <- ggplot(ratio_df, aes(x = sample, y = ratio, fill = is_A1)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", ratio)), vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#0072B2"), guide = "none") +
  scale_y_log10(expand = expansion(mult = c(0.02, 0.18))) +
  labs(title = "B. Podocyte:tubular marker-sum ratio per ADR sample (log scale)",
       subtitle = sprintf("A-ADR1 ratio = %.2f, the lowest of all 6 ADR samples (next-lowest A-ADR2 = %.2f, %.1fx higher)",
                           ratio_df$ratio[ratio_df$sample == "A-ADR1"],
                           sort(ratio_df$ratio[ratio_df$sample != "A-ADR1"])[1],
                           sort(ratio_df$ratio[ratio_df$sample != "A-ADR1"])[1] / ratio_df$ratio[ratio_df$sample == "A-ADR1"]),
       x = NULL, y = "Podocyte:tubular\nmarker-sum ratio") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9),
        axis.title.y = element_text(size = 11))

p <- gridExtra::arrangeGrob(panelA, panelB, ncol = 1, heights = c(1.3, 1))

for (ext in c("png", "pdf")) {
  fp <- file.path(outdir, paste0("figures/FigS_A1_contamination_QC.", ext))
  ggsave(fp, plot = p, width = 9, height = 10.5, dpi = 300)
  cat("Saved:", fp, "\n")
}

# ---------------------------------------------------------------------------
# text summary
# ---------------------------------------------------------------------------
tub_sum_others_mean <- mean(tubule_sum[names(tubule_sum) != "A-ADR1"])
fold_tub <- tubule_sum["A-ADR1"] / tub_sum_others_mean
fold_ratio <- sort(ratio[names(ratio) != "A-ADR1"])[1] / ratio["A-ADR1"]

summary_lines <- c(
  "=== A-ADR1 contamination QC summary (recomputed from tables/00_merged_counts.tsv) ===",
  "",
  "Summed tubular-marker CPM per ADR sample:",
  sprintf("  %s: %.1f", names(tubule_sum), tubule_sum),
  "",
  "Summed podocyte-marker CPM per ADR sample:",
  sprintf("  %s: %.1f", names(podo_sum), podo_sum),
  "",
  "Podocyte:tubular marker-sum ratio per ADR sample:",
  sprintf("  %s: %.2f", names(ratio), ratio),
  "",
  sprintf("A-ADR1 summed tubular-marker CPM (%.1f) is %.2fx the mean of the other 5 ADR samples (%.1f).",
          tubule_sum["A-ADR1"], fold_tub, tub_sum_others_mean),
  sprintf("A-ADR1 podocyte:tubular ratio (%.2f) is the lowest of all 6 ADR samples, %.2fx lower than the next-lowest (A-ADR2, %.2f).",
          ratio["A-ADR1"], fold_ratio, sort(ratio[names(ratio) != "A-ADR1"])[1]),
  "A-ADR1 summed podocyte-marker CPM is NOT the lowest among ADR samples (mid-range) -- the low ratio is driven by excess tubular signal, not depleted podocyte signal.",
  "Per-marker: A-ADR1 has the highest CPM among the 6 ADR samples for 5/6 tubular markers (all except Lrp2, where A-ADR2 is marginally higher)."
)
writeLines(summary_lines, file.path(outdir, "logs/18_A1_contamination_QC_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

cat("\n=== 18_A1_contamination_QC_figure.R complete ===\n")
