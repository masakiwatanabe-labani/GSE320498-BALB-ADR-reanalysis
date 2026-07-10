#!/usr/bin/env Rscript
# Figure 6, submission-ready, all 4 panels from one canonical run:
#   A: Day-5 volcano (DESeq2, ADR_B_vs_A, A1-excluded) -- Wt1/Nphs1 +
#      Serpine1/Loxl1/Col4a1/Col4a2 labeled, gene name only
#   B: Karaiskos TOP50 podocyte-marker preranked GSEA (replaces the
#      Tabula Muris Senis ageing set)
#   C: Reactome ORA on the DESeq2-derived Day-5 DEG list (replaces the
#      edgeR-derived ORA), significant ECM/collagen terms only
#   D: Reactome Integrin Signaling preranked GSEA
#
# No new statistical model anywhere in this script: A reads
# DE_ADR_B_vs_A.tsv verbatim; B/D reuse the already-exported fgsea running-ES
# curve/hit/ranking-vector data from 10_export_enrichment_curves.R (verified
# against tables/enrichment_curves/panel_metadata.tsv, which itself is the
# canonical fgsea output); C reuses the ORA table from
# 25_Fig6C_ORA_rerun_DESeq2_DEGs.R. This script only visualizes.

suppressMessages({
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

indir <- "/usr/local/jupyter/ADR_BALB"
outdir_main <- file.path(indir, "outputs")
outdir_v2 <- file.path(outdir_main, "REVISION_PACKAGE/REVISION_PACKAGE_v2")
fig_dir <- file.path(outdir_v2, "figures")
tab_dir <- file.path(outdir_v2, "tables")
log_dir <- file.path(outdir_v2, "logs")
curve_dir <- file.path(outdir_main, "tables/enrichment_curves")

pt_mm <- function(pt) pt / 2.845276

check <- c("=== Figure 6 (A-D) consistency check ===", "")

# ============================================================================
# shared: canonical gene-count check
# ============================================================================
d5 <- read.delim(file.path(outdir_main, "tables/DE_ADR_B_vs_A.tsv"), stringsAsFactors = FALSE)
n_total <- nrow(d5)
check <- c(check, sprintf("Panel A source (DE_ADR_B_vs_A.tsv) tested genes: %d (expected 19,662): %s",
                           n_total, ifelse(n_total == 19662, "MATCH", "*** MISMATCH ***")))

rank_vec <- read.delim(file.path(curve_dir, "ranking_vector_ADR_B_vs_A_stat.tsv"), stringsAsFactors = FALSE)
check <- c(check, sprintf("Panels B/D ranking vector (ranking_vector_ADR_B_vs_A_stat.tsv) length: %d (expected 19,662): %s",
                           nrow(rank_vec), ifelse(nrow(rank_vec) == 19662, "MATCH", "*** MISMATCH ***")))

panel_meta <- read.delim(file.path(curve_dir, "panel_metadata.tsv"), stringsAsFactors = FALSE)
check <- c(check, sprintf("Panel metadata universe (n_genes_total column): %s (all rows expected 19,662): %s",
                           paste(unique(panel_meta$n_genes_total), collapse = ","),
                           ifelse(all(panel_meta$n_genes_total == 19662), "MATCH", "*** MISMATCH ***")))

ora_full <- read.csv(file.path(tab_dir, "TableS_ORA_D5_full.csv"), stringsAsFactors = FALSE)
check <- c(check, sprintf("Panel C ORA universe (from 25_Fig6C_ORA_rerun_DESeq2_DEGs.R log): 19,662 (recorded, not re-derivable from this CSV alone) -- see logs/25_Fig6C_ORA_versions_and_summary.txt"))

# ============================================================================
# PANEL A: Day-5 volcano
# ============================================================================
always_label <- c("Wt1", "Nphs1")
ecm_label <- c("Serpine1", "Loxl1", "Col4a1", "Col4a2")
all_label_A <- c(always_label, ecm_label)

d5$negLog10FDR <- -log10(d5$padj)
d5$category <- with(d5, ifelse(is.na(padj), "NS",
                          ifelse(padj < 0.05 & abs(log2FC) > 1, "FDR and log2FC",
                          ifelse(padj < 0.05, "FDR only",
                          ifelse(abs(log2FC) > 1, "log2FC only", "NS")))))
d5$category <- factor(d5$category, levels = c("NS", "log2FC only", "FDR only", "FDR and log2FC"))
d5$point_label <- ifelse(d5$gene %in% all_label_A, d5$gene, NA)

wt1_d5 <- d5[d5$gene == "Wt1", ]; nphs1_d5 <- d5[d5$gene == "Nphs1", ]
check <- c(check, "",
  sprintf("Panel A Wt1: log2FC=%.4f, padj=%.4f (expected padj~0.0498): %s", wt1_d5$log2FC, wt1_d5$padj,
          ifelse(abs(wt1_d5$padj - 0.0498) < 0.002, "MATCH", "*** MISMATCH ***")),
  sprintf("Panel A Nphs1: log2FC=%.4f, padj=%.4f (expected padj~0.064): %s", nphs1_d5$log2FC, nphs1_d5$padj,
          ifelse(abs(nphs1_d5$padj - 0.064) < 0.003, "MATCH", "*** MISMATCH ***")))
na_labeled_A <- d5[d5$gene %in% all_label_A & is.na(d5$padj), ]
check <- c(check, sprintf("Panel A labeled genes with padj=NA (should be 0): %d -- %s", nrow(na_labeled_A),
                           ifelse(nrow(na_labeled_A) == 0, "OK", "*** ERROR ***")))

cat_colors <- c("NS" = "grey75", "log2FC only" = "#0072B2", "FDR only" = "#009E73", "FDR and log2FC" = "#D55E00")
d5_bg <- d5[!(d5$gene %in% all_label_A), ]
d5_hi <- d5[d5$gene %in% ecm_label, ]
d5_always <- d5[d5$gene %in% always_label, ]

d5_always_lbl <- d5[d5$gene %in% always_label, ]; d5_always_lbl <- d5_always_lbl[order(d5_always_lbl$gene), ]
d5_always_lbl$nudge_x <- ifelse(d5_always_lbl$gene == "Wt1", -4.8, 4.2)
d5_always_lbl$nudge_y <- c(7, 11)
d5_always_lbl$text_color <- "#9B4A73"
d5_hi_lbl <- d5[d5$gene %in% ecm_label, ]
d5_hi_lbl$nudge_x <- 0; d5_hi_lbl$nudge_y <- 0
d5_hi_lbl$text_color <- "black"
d5_labels_all <- rbind(d5_hi_lbl, d5_always_lbl)

x_range_A <- range(d5$log2FC, na.rm = TRUE)

panelA <- ggplot(d5, aes(x = log2FC, y = negLog10FDR)) +
  geom_point(data = d5_bg, aes(color = category), size = 0.9, alpha = 0.6) +
  scale_color_manual(values = cat_colors, name = NULL, drop = FALSE) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.35) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40", linewidth = 0.35) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey88", linewidth = 0.3) +
  geom_point(data = d5_hi, aes(fill = "ECM genes"), shape = 23, color = "black", size = 2.3, stroke = 0.55) +
  geom_point(data = d5_always, aes(fill = "Wt1 / Nphs1"), shape = 24, color = "black", size = 3.0, stroke = 0.7) +
  scale_fill_manual(values = c("ECM genes" = "#F0E442", "Wt1 / Nphs1" = "#CC79A7"), name = NULL) +
  ggrepel::geom_text_repel(
    data = d5_labels_all, aes(label = point_label, colour = I(text_color), segment.colour = I(text_color)),
    nudge_x = d5_labels_all$nudge_x, nudge_y = d5_labels_all$nudge_y,
    size = pt_mm(11.5), fontface = "bold", max.overlaps = Inf, box.padding = 1.0, point.padding = 0.55,
    min.segment.length = 0, segment.size = 0.35, force = 14, force_pull = 0.2, max.time = 5, max.iter = 60000, seed = 20260220
  ) +
  annotate("text", x = x_range_A[1], y = -Inf, label = "← higher in AJcl", hjust = 0, vjust = -0.6, size = pt_mm(10.5), fontface = "italic", color = "grey30") +
  annotate("text", x = x_range_A[2], y = -Inf, label = "higher in ByJcl →", hjust = 1, vjust = -0.6, size = pt_mm(10.5), fontface = "italic", color = "grey30") +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.1))) +
  scale_x_continuous(expand = expansion(mult = c(0.07, 0.07))) +
  labs(title = "Day 5 post-ADR: BALB/cByJcl vs BALB/cAJcl (A-ADR1 excluded)",
       x = expression(log[2]~"FC (ByJcl vs AJcl)"), y = expression(-log[10]~"(FDR)"),
       caption = sprintf("total = %s variables", format(n_total, big.mark = ","))) +
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 2.6)),
         fill = guide_legend(nrow = 2, override.aes = list(size = 3.2))) +
  theme_bw(base_size = 13.5) +
  theme(plot.title = element_text(size = 14.5, face = "bold"),
        plot.caption = element_text(size = 10.5, hjust = 0, face = "italic"),
        axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 11.5),
        legend.text = element_text(size = 11), legend.title = element_text(size = 11),
        legend.key.size = unit(0.45, "cm"),
        legend.position = "bottom", legend.margin = margin(0, 0, 0, 0), legend.box = "vertical")

# ============================================================================
# helper: classic 3-track GSEA panel (ES curve / hit barcode / ranking strip)
# ============================================================================
make_gsea_panel <- function(curve_file, hits_file, nes, fdr, gene_set_label, title) {
  curve <- read.delim(file.path(curve_dir, curve_file), stringsAsFactors = FALSE)
  hits <- read.delim(file.path(curve_dir, hits_file), stringsAsFactors = FALSE)
  n_max <- max(curve$x, na.rm = TRUE)

  fdr_str <- formatC(fdr, format = "e", digits = 2)
  ann_lines <- sprintf("NES = %+.2f\nFDR = %s", nes, fdr_str)

  # anchor the NES/FDR box to the plot's top-margin (expansion headroom),
  # not to curve data min/max -- robust regardless of curve shape/sign
  # (a data-anchored position previously clipped for one-directional curves
  # whose min/max sits right at the panel edge, e.g. Panel D).
  p_es <- ggplot(curve, aes(x = x, y = y)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_line(color = "#1B7837", linewidth = 0.9) +
    annotate("label", x = n_max * 0.6, y = Inf, vjust = 1.3,
             label = ann_lines, size = pt_mm(12), fontface = "bold", hjust = 0) +
    labs(title = title, y = "Enrichment\nscore (ES)") +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.22))) +
    theme_bw(base_size = 13.5) +
    theme(plot.title = element_text(size = 14.5, face = "bold"),
          axis.title.y = element_text(size = 12.5, face = "bold"), axis.text.y = element_text(size = 11),
          axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.margin = margin(4, 4, 0, 4))

  p_hits <- ggplot(hits, aes(x = hit_position)) +
    geom_segment(aes(xend = hit_position, y = 0, yend = 1), color = "black", linewidth = 0.4) +
    scale_x_continuous(limits = c(0, n_max)) +
    scale_y_continuous(breaks = NULL) +
    labs(y = NULL) +
    theme_void(base_size = 13.5) +
    theme(axis.text = element_blank(), plot.margin = margin(0, 4, 0, 4),
          panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.3))

  rv <- rank_vec
  n_bins <- 250
  rv$bin <- pmin(n_bins, ceiling(rv$rank_position / (nrow(rv) / n_bins)))
  bin_df <- aggregate(stat ~ bin, data = rv, FUN = mean)
  bin_df$xmin <- (bin_df$bin - 1) * (n_max / n_bins)
  bin_df$xmax <- bin_df$bin * (n_max / n_bins)
  max_abs_stat <- max(abs(bin_df$stat))

  p_rank <- ggplot(bin_df) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1, fill = stat)) +
    scale_fill_gradient2(low = "#08519C", mid = "white", high = "#A50026", midpoint = 0,
                          limits = c(-max_abs_stat, max_abs_stat), guide = "none") +
    scale_x_continuous(limits = c(0, n_max), expand = c(0, 0)) +
    scale_y_continuous(breaks = NULL, expand = c(0, 0)) +
    labs(x = "Rank in ordered gene list\n(DESeq2 Wald statistic, Day 5 ADR ByJcl vs. AJcl)", y = NULL) +
    theme_void(base_size = 13.5) +
    theme(axis.title.x = element_text(size = 10.5, margin = margin(t = 3)),
          plot.margin = margin(0, 4, 4, 4))

  combined <- p_es / p_hits / p_rank + plot_layout(heights = c(6, 0.7, 1))
  combined
}

panelB <- make_gsea_panel(
  "curve_KARAISKOS_TOP50_ADR_B_vs_A.tsv", "hits_KARAISKOS_TOP50_ADR_B_vs_A.tsv",
  nes = -2.0438, fdr = 0.0005246, gene_set_label = "KARAISKOS2018_PODOCYTE_TOP50",
  title = "Karaiskos et al. 2018 podocyte markers (TOP50)"
)
panelD <- make_gsea_panel(
  "curve_INTEGRIN_SIGNALING.tsv", "hits_INTEGRIN_SIGNALING.tsv",
  nes = 2.0260, fdr = 0.0009105, gene_set_label = "REACTOME_INTEGRIN_SIGNALING",
  title = "Reactome Integrin Signaling"
)

check <- c(check, "",
  "Panel B (Karaiskos TOP50) NES=-2.0438, FDR=5.246e-04 (~5.2e-4 at 2sf, 5.25e-4 at 3sf).",
  "  NOTE: task text stated 'FDR = 5.3e-4' -- precise canonical value is 5.25e-04 (rounds to 5.2e-4",
  "  at 2 significant figures, matching the pre-existing Fig6_enrichment_KARAISKOS_TOP50.png box",
  "  exactly, which shows 'joint FDR (q) = 5.25e-04'). Figure legend below uses 5.25x10^-4.",
  sprintf("Panel D (Integrin Signaling) NES=+2.0260, FDR=9.105e-04 (~9.1e-4): MATCHES task text exactly.")
)

# ============================================================================
# PANEL C: ORA dot plot, significant ECM/collagen terms only
# ============================================================================
ecm_pattern <- "ECM|COLLAGEN|EXTRACELLULAR MATRIX|EXTRACELLULAR_MATRIX|INTEGRIN|CELL ADHESION|LAMININ|BASEMENT MEMBRANE"
ecm_ora <- ora_full[grepl(ecm_pattern, ora_full$Description, ignore.case = TRUE), ]
ecm_ora_sig <- ecm_ora[!is.na(ecm_ora$p.adjust) & ecm_ora$p.adjust < 0.05, ]
ecm_ora_sig <- ecm_ora_sig[order(ecm_ora_sig$p.adjust), ]

integrin_csi <- ecm_ora[grepl("Integrin cell surface", ecm_ora$Description), ]
check <- c(check, "",
  sprintf("Panel C: %d ECM/collagen-pattern ORA terms significant (BH<0.05) of %d matched.",
          nrow(ecm_ora_sig), nrow(ecm_ora)),
  sprintf("Panel C: 'Integrin cell surface interactions' BH p.adjust=%.4f -- %s (excluded from panel per spec).",
          integrin_csi$p.adjust, ifelse(integrin_csi$p.adjust >= 0.05, "correctly NOT significant", "*** now significant, check spec ***"))
)

parse_ratio <- function(x) sapply(strsplit(x, "/"), function(v) as.numeric(v[1]) / as.numeric(v[2]))
ecm_ora_sig$GeneRatioNum <- parse_ratio(ecm_ora_sig$GeneRatio)
ecm_ora_sig$Description <- factor(ecm_ora_sig$Description, levels = rev(ecm_ora_sig$Description[order(ecm_ora_sig$GeneRatioNum)]))

panelC <- ggplot(ecm_ora_sig, aes(x = GeneRatioNum, y = Description)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_viridis_c(option = "viridis", direction = -1, name = "BH-adj.\np", trans = "log10",
                         guide = guide_colorbar(raster = FALSE, nbin = 20, frame.colour = "black", frame.linewidth = 0.3, ticks.colour = "black")) +
  scale_size_continuous(name = "Count", range = c(3, 9.5)) +
  labs(title = "Reactome ORA: significant ECM/collagen terms",
       subtitle = "DESeq2-derived DEGs, padj<0.05, n=546",
       x = "GeneRatio", y = NULL) +
  theme_bw(base_size = 13.5) +
  theme(plot.title = element_text(size = 14.5, face = "bold"),
        plot.subtitle = element_text(size = 10.5),
        axis.title.x = element_text(size = 14, face = "bold"), axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11), legend.text = element_text(size = 11),
        legend.title = element_text(size = 11))

# ============================================================================
# combine, tag panels, save
# ============================================================================
combined <- (panelA | patchwork::wrap_elements(panelB)) / (panelC | patchwork::wrap_elements(panelD)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 17, face = "bold"))

ggsave(file.path(fig_dir, "Figure6_combined.pdf"), plot = combined, width = 17, height = 13.5, device = cairo_pdf)
ggsave(file.path(fig_dir, "Figure6_combined.png"), plot = combined, width = 17, height = 13.5, dpi = 300)
cat("Saved: Figure6_combined.pdf/png\n")

for (nm in c("panelA", "panelB", "panelC", "panelD")) {
  p <- get(nm)
  tag <- toupper(substr(nm, 6, 6))
  fp_pdf <- file.path(fig_dir, paste0("Figure6", tag, ".pdf"))
  fp_png <- file.path(fig_dir, paste0("Figure6", tag, ".png"))
  w <- if (tag %in% c("B", "D")) 8.5 else if (tag == "C") 11 else 10
  h <- if (tag %in% c("B", "D")) 7 else 7
  ggsave(fp_pdf, plot = p, width = w, height = h, device = cairo_pdf)
  ggsave(fp_png, plot = p, width = w, height = h, dpi = 300)
  cat("Saved:", fp_pdf, "/", fp_png, "\n")
}

writeLines(check, file.path(log_dir, "31_Figure6_consistency_check.txt"))
cat(paste(check, collapse = "\n"), "\n")

# ============================================================================
# caption + package versions
# ============================================================================
caption <- paste0(
  "Figure 6. Canonical DESeq2/fgsea/ReactomePA re-analysis of Day-5 post-adriamycin substrain differences ",
  "(BALB/cByJcl vs BALB/cAJcl, A-ADR1 excluded, n=3 ByJcl vs n=2 AJcl). All panels from the same run: ",
  "19,662 tested genes, positive log2FC/NES = higher in ByJcl. ",
  "(A) Volcano plot (DESeq2 Wald test). Horizontal line, FDR=0.05; vertical lines, |log2FC|=1. ",
  "Podocyte markers Wt1 and Nphs1 and ECM/adhesion transcripts Serpine1, Loxl1, Col4a1, Col4a2 are labeled. ",
  "(B) Preranked GSEA (fgsea v1.24.0, DESeq2 Wald-statistic ranking) of the Karaiskos et al. (2018, J Am Soc ",
  "Nephrol 29:2060-2068, GEO GSE111107) podocyte marker panel (TOP50, n=49 genes), replacing the Tabula Muris ",
  "Senis podocyte-ageing signature used in the original Fig. 6B (whose sign is not robust to A-ADR1 status; ",
  "see Discussion). NES=-2.04, joint FDR=5.25x10-4. Leading-edge genes: Nphs2, Cdkn1c, Clic3, Nupr1, Dpp4, ",
  "Enpep, Tcf21, Nphs1, Gadd45a, Rab3b, Rhpn1, Tmsb4x, Col4a3, Rasl11a, Mafb, Npnt, Arhgap24, Adm, Pak1, Synpo, ",
  "Foxd2os, Golim4, Igfbp7, Vegfa, Cd59a, Sdc4, Sema3g, Tdrd5, Nap1l1, Shisa3, Eif3m, Thsd7a, Pth1r, Wt1. ",
  "(C) Reactome over-representation analysis (ReactomePA::enrichPathway v1.42.0) of the DESeq2-derived Day-5 ",
  "DEG list (padj<0.05, n=546 genes; background = all 19,662 tested genes), replacing the original ",
  "edgeR-derived ORA. 19 of 734 tested Reactome terms reached BH-adjusted p<0.05; the 7 significant ",
  "ECM/collagen-related terms are shown (dot size = gene count, color = BH-adjusted p-value). ",
  "'Integrin cell surface interactions' did not reach significance under ORA (BH p=0.057) and is not shown. ",
  "(D) Preranked GSEA of the Reactome Integrin Signaling gene set (n=27 genes). NES=+2.03, joint FDR=9.1x10-4. ",
  "Leading-edge genes: Bcar1, Akt1, Rapgef4, Shc1, Ptpn1, Syk, Rap1b, Fn1, Itgb3, Rap1a, Src, Crk, Tln1, Fga, ",
  "Ptk2, Apbb1ip, Csk."
)
writeLines(caption, file.path(log_dir, "31_Figure6_caption_draft.txt"))
cat("\n=== Figure 6 caption (draft) ===\n", caption, "\n")

version_log <- c(
  "=== Figure 6 (A-D) package versions and seed ===",
  sprintf("ggplot2: %s", as.character(packageVersion("ggplot2"))),
  sprintf("ggrepel: %s", as.character(packageVersion("ggrepel"))),
  sprintf("patchwork: %s", as.character(packageVersion("patchwork"))),
  "Upstream analysis (reused, not re-run): DESeq2 1.38.3, fgsea 1.24.0, msigdbr 26.1.0, ReactomePA 1.42.0, clusterProfiler 4.6.2, org.Mm.eg.db 3.16.0",
  "Random seed: 20260220 (fgsea/ORA runs that produced the underlying data)"
)
writeLines(version_log, file.path(log_dir, "31_Figure6_versions.txt"))

cat("\n=== 31_Figure6_combined_final.R complete ===\n")
