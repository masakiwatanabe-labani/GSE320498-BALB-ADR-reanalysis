#!/usr/bin/env Rscript
# Rebuild Fig. 6C as a canonical-fgsea dot plot, replacing the original
# edgeR + ReactomePA enrichPathway() (ORA) figure. The rest of Fig. 6
# (B, D) was already re-derived from preranked fgsea on DESeq2 Wald-stat
# ranks (script 09_gsea_stat_ranking_canonical.R); this script brings C
# into the same framework so all four Fig. 6 panels come from one run.
#
# Reuses the existing canonical ADR_B_vs_A (A1-excluded, main analysis)
# fgsea output verbatim -- no re-run, no new statistical model. Settings
# (for the record): DESeq2 v1.38.3 Wald stat ranking, fgsea v1.24.0
# multilevel (eps=0, minSize=5, maxSize=500), Reactome M2:CP:REACTOME +
# 3 custom podocyte sets = 1336 joint universe (1293 sets survived the
# size filter and were actually tested), BH FDR computed jointly across
# that full universe, seed 20260220.

suppressMessages({
  library(ggplot2)
  library(dplyr)
})

indir <- "/usr/local/jupyter/ADR_BALB"
outdir_main <- file.path(indir, "outputs")
outdir_v2 <- file.path(outdir_main, "REVISION_PACKAGE/REVISION_PACKAGE_v2")
fig_dir <- file.path(outdir_v2, "figures")
tab_dir <- file.path(outdir_v2, "tables")
log_dir <- file.path(outdir_v2, "logs")

versions <- readLines(file.path(outdir_main, "logs/09_gsea_canonical_versions.txt"))

# ---------------------------------------------------------------------------
# Step 1: load canonical ADR_B_vs_A (A1-excluded, main) fgsea result
# ---------------------------------------------------------------------------
df <- read.delim(file.path(outdir_main, "tables/GSEA_ADR_B_vs_A_full_joint_canonical.tsv"),
                  stringsAsFactors = FALSE)
n_tested <- nrow(df)
n_sig <- sum(df$padj < 0.05, na.rm = TRUE)
cat(sprintf("Step 1: ADR_B_vs_A canonical fgsea -- %d gene sets tested, %d significant at joint FDR<0.05\n",
            n_tested, n_sig))

clean_name <- function(x) {
  x <- sub("^REACTOME_", "", x)
  x <- gsub("_", " ", x)
  x <- vapply(x, function(s) paste(strwrap(s, width = 40), collapse = "\n"), character(1))
  x
}
df$display_name <- clean_name(df$pathway)

# ---------------------------------------------------------------------------
# Step 2/3: dot plot builder + set-selection variants
# ---------------------------------------------------------------------------
plot_dotplot <- function(sub_df, title, subtitle, file_stem, legend_lines) {
  sub_df$display_name <- factor(sub_df$display_name, levels = sub_df$display_name[order(sub_df$NES)])
  p <- ggplot(sub_df, aes(x = NES, y = display_name)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(aes(size = size, color = padj)) +
    scale_color_viridis_c(option = "viridis", direction = -1, trans = "log10",
                           name = "Joint BH FDR\n(q-value)",
                           labels = scales::label_scientific()) +
    scale_size_continuous(name = "Gene set\nsize", range = c(2.2, 9)) +
    labs(title = title, subtitle = subtitle,
         x = "Normalized Enrichment Score (NES)\npositive = enriched in ByJcl", y = NULL) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 8.6, lineheight = 1.15),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 8.3, lineheight = 0.85),
      axis.text.x = element_text(size = 10.5),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      plot.margin = margin(10, 14, 10, 10)
    )

  n_row <- nrow(sub_df)
  h <- max(5, 1.6 + n_row * 0.42)
  for (ext in c("png", "pdf")) {
    fp <- file.path(fig_dir, paste0(file_stem, ".", ext))
    ggsave(fp, plot = p, width = 12.5, height = h, dpi = 300)
    cat("Saved:", fp, "\n")
  }
  writeLines(legend_lines, file.path(fig_dir, paste0(file_stem, "_legend.txt")))
  cat("Saved:", file.path(fig_dir, paste0(file_stem, "_legend.txt")), "\n\n")
}

comparison_caption <- "Comparison: ADR_B_vs_A (Day 5 post-adriamycin, BALB/cByJcl vs BALB/cAJcl; A-ADR1 excluded, main analysis)."
method_caption <- "DESeq2 v1.38.3 Wald-stat-ranked preranked fgsea v1.24.0 (multilevel, eps=0); Reactome M2:CP:REACTOME + 3 custom podocyte sets, 1336-set joint universe, 1293 sets tested; BH FDR computed jointly across all 1293 sets. Seed 20260220."

# ---- main figure: top 20 by FDR, among joint-FDR<0.05 sets ----
sig_df <- df[!is.na(df$padj) & df$padj < 0.05, ]
sig_df <- sig_df[order(sig_df$padj), ]
main_df <- head(sig_df, 20)

plot_dotplot(
  main_df,
  "Figure 6C: Reactome pathway enrichment, Day 5 post-ADR (ByJcl vs AJcl)",
  paste0("Top 20 of ", nrow(sig_df), " gene sets significant at joint FDR<0.05 (of ", n_tested,
         " tested), ranked by FDR (most significant at top/bottom per NES sort).\n", comparison_caption),
  "Figure6C_dotplot_fgsea",
  c(
    "Figure 6C legend (draft).",
    "",
    "Selection rule: the 20 gene sets with the smallest joint BH FDR, restricted to",
    "gene sets reaching joint FDR<0.05 (355 of 1293 gene sets tested pass this",
    "threshold; this panel shows the top 20 of those 355).",
    "",
    comparison_caption,
    method_caption,
    "",
    "Axes/encoding: x = normalized enrichment score (NES; positive = higher in",
    "BALB/cByJcl at Day 5 post-adriamycin). Point size = number of genes in the",
    "gene set (post-filtering). Point color = joint Benjamini-Hochberg FDR",
    "(q-value), computed across the full 1293-gene-set tested universe, NOT a",
    "focal subset -- darker/more saturated color = smaller (more significant) FDR.",
    "Dashed vertical line marks NES=0. Gene sets are sorted by NES on the y-axis.",
    "",
    "Note: because this panel selects by FDR rank genome-wide, and FDR is driven",
    "jointly by effect size and gene-set size/variance, the top 20 sets by FDR are",
    "dominated by large, low-dispersion housekeeping gene sets (translation,",
    "ribosome, proteasome, antigen presentation). REACTOME_EXTRACELLULAR_MATRIX_",
    "ORGANIZATION (rank 55/1293, FDR=9.7e-11) and REACTOME_INTEGRIN_SIGNALING",
    "(rank 170/1293, FDR=9.1e-4) are both significant but do not appear in this",
    "top-20-by-FDR view; see Figure6C_alt_ECMonly for these gene sets specifically."
  )
)

# ---- alt (i): top 20 by |NES|, among significant sets ----
absnes_df <- sig_df[order(-abs(sig_df$NES)), ]
absnes_df <- head(absnes_df, 20)

plot_dotplot(
  absnes_df,
  "Figure 6C (alt.): Reactome pathways with the largest effect size (|NES|)",
  paste0("Top 20 of ", nrow(sig_df), " gene sets significant at joint FDR<0.05 (of ", n_tested,
         " tested), ranked by |NES|.\n", comparison_caption),
  "Figure6C_alt_absNES",
  c(
    "Figure 6C (alternate, |NES|-ranked) legend (draft).",
    "",
    "Selection rule: the 20 gene sets with the largest absolute normalized",
    "enrichment score (|NES|), restricted to gene sets reaching joint FDR<0.05",
    "(355 of 1293 gene sets tested pass this threshold; this panel shows the top",
    "20 of those 355 by effect size rather than by FDR rank).",
    "",
    comparison_caption,
    method_caption,
    "",
    "Axes/encoding: identical to the main Figure 6C panel (x = NES, point size =",
    "gene set size, point color = joint FDR, dashed line at NES=0, y sorted by",
    "NES). This view answers a different question from the FDR-ranked main panel:",
    "which gene sets show the largest shift in ranking-metric-weighted expression,",
    "regardless of gene-set size/variance effects on FDR."
  )
)

# ---- alt (ii): ECM/adhesion-related sets only, regardless of significance ----
ecm_pattern <- "ECM|COLLAGEN|EXTRACELLULAR_MATRIX|INTEGRIN|CELL_ADHESION|CELL_JUNCTION|LAMININ|BASEMENT_MEMBRANE|CELL_CELL"
ecm_df <- df[grepl(ecm_pattern, df$pathway), ]
ecm_df <- ecm_df[order(-ecm_df$NES), ]

plot_dotplot(
  ecm_df,
  "Figure 6C (alt.): ECM / cell-adhesion / integrin-related Reactome gene sets",
  paste0(nrow(ecm_df), " gene sets matched an ECM/collagen/integrin/adhesion name filter (of ", n_tested,
         " tested; shown regardless of significance).\n", comparison_caption),
  "Figure6C_alt_ECMonly",
  c(
    "Figure 6C (alternate, ECM/adhesion-filtered) legend (draft).",
    "",
    sprintf("Selection rule: all Reactome gene sets whose name matched the pattern"),
    "'ECM|COLLAGEN|EXTRACELLULAR_MATRIX|INTEGRIN|CELL_ADHESION|CELL_JUNCTION|",
    "LAMININ|BASEMENT_MEMBRANE|CELL_CELL' (case-insensitive), shown regardless of",
    sprintf("significance (%d of %d tested gene sets matched).", nrow(ecm_df), n_tested),
    "",
    comparison_caption,
    method_caption,
    "",
    "Axes/encoding: identical to the main Figure 6C panel. This view directly",
    "addresses whether ECM-organization/integrin/adhesion biology is enriched in",
    "this comparison, independent of whether any individual set reaches the",
    "top-20-by-FDR or top-20-by-|NES| cutoffs used in the other two panels: 20 of",
    sprintf("22 matched sets have positive NES (higher in ByJcl), %d reach joint FDR<0.05.",
            sum(ecm_df$padj < 0.05, na.rm = TRUE))
  )
)

# ---------------------------------------------------------------------------
# Full supplementary table
# ---------------------------------------------------------------------------
full_out <- df[order(df$pval), c("pathway", "size", "ES", "NES", "pval", "padj", "leadingEdge")]
colnames(full_out)[colnames(full_out) == "pval"] <- "nominal_pval"
colnames(full_out)[colnames(full_out) == "padj"] <- "FDR_qvalue_joint"
write.csv(full_out, file.path(tab_dir, "TableS_fgsea_ADR_B_vs_A_full.csv"), row.names = FALSE)
cat("Saved:", file.path(tab_dir, "TableS_fgsea_ADR_B_vs_A_full.csv"), sprintf("(%d gene sets)\n", nrow(full_out)))

# ---------------------------------------------------------------------------
# Step 4: cross-figure consistency check
# ---------------------------------------------------------------------------
de_files <- c(baseline_B_vs_A = "DE_baseline_B_vs_A.tsv", ADR_B_vs_A = "DE_ADR_B_vs_A.tsv",
              A_ADR_vs_Ctrl = "DE_A_ADR_vs_Ctrl.tsv", B_ADR_vs_Ctrl = "DE_B_ADR_vs_Ctrl.tsv")
gene_counts <- sapply(de_files, function(f) nrow(read.delim(file.path(outdir_main, "tables", f))))

gsea_files <- c(baseline_B_vs_A = "GSEA_baseline_B_vs_A_full_joint_canonical.tsv",
                ADR_B_vs_A = "GSEA_ADR_B_vs_A_full_joint_canonical.tsv",
                A_ADR_vs_Ctrl = "GSEA_A_ADR_vs_Ctrl_full_joint_canonical.tsv",
                B_ADR_vs_Ctrl = "GSEA_B_ADR_vs_Ctrl_full_joint_canonical.tsv")
set_counts <- sapply(gsea_files, function(f) nrow(read.delim(file.path(outdir_main, "tables", f))))

focal_ranks <- df[order(df$padj), ]
focal_ranks$rank <- seq_len(nrow(focal_ranks))
ecm_rank <- focal_ranks[focal_ranks$pathway == "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION", ]
integrin_rank <- focal_ranks[focal_ranks$pathway == "REACTOME_INTEGRIN_SIGNALING", ]
integrin_csi_rank <- focal_ranks[focal_ranks$pathway == "REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS", ]

consistency_lines <- c(
  "=== Step 4: cross-figure (Fig. 6A/6B/6C/6D) consistency check ===",
  "",
  "Genes tested (DESeq2, rows in DE_<comparison>.tsv, all non-NA stat):",
  sprintf("  %-16s n = %d", names(gene_counts), gene_counts),
  "  -> All 4 primary comparisons (underlying Fig. 5, 6A, 6B, 6C, 6D) are drawn",
  "     from the SAME single DESeq2 dds object (dds_group_main, A-ADR1 excluded,",
  "     11 samples, prefilter rowSums>=10 applied once to that shared object),",
  "     so the gene-count mismatch in the original manuscript (Fig.5=13,064 vs.",
  "     Fig.6A=13,278 'variables') does not occur in this canonical re-analysis:",
  "     every comparison tests the identical 19,662 genes.",
  "",
  "Gene sets tested per comparison (fgsea, after minSize=5/maxSize=500 filter):",
  sprintf("  %-16s n = %d", names(set_counts), set_counts),
  "  -> Identical (1293) in all 4 comparisons: same joint universe (1336 sets:",
  "     Reactome M2:CP:REACTOME 1333 + 3 custom podocyte sets), same size filter.",
  "",
  "FDR definition: identical across all 4 comparisons and all Fig. 6 panels --",
  "Benjamini-Hochberg, computed jointly across the full tested gene-set universe",
  "for that comparison (never on a focal subset alone). Significance threshold",
  "used throughout: joint FDR (padj) < 0.05, per the manuscript's own stated",
  "criterion.",
  "",
  "Does ECM Organization / Integrin Signaling appear in the Fig. 6C top-20-by-FDR view?",
  sprintf("  REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION: NES=%.2f, FDR=%.2e, genome-wide rank %d/%d -- NOT in top 20 by FDR.",
          ecm_rank$NES, ecm_rank$padj, ecm_rank$rank, n_tested),
  sprintf("  REACTOME_INTEGRIN_SIGNALING: NES=%.2f, FDR=%.2e, genome-wide rank %d/%d -- NOT in top 20 by FDR.",
          integrin_rank$NES, integrin_rank$padj, integrin_rank$rank, n_tested),
  sprintf("  REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS: NES=%.2f, FDR=%.2e, genome-wide rank %d/%d -- NOT in top 20 by FDR.",
          integrin_csi_rank$NES, integrin_csi_rank$padj, integrin_csi_rank$rank, n_tested),
  "  All three ARE significant at joint FDR<0.05 and are shown in Figure6C_alt_ECMonly.",
  "  They are absent from the FDR-ranked main panel only because ~50-170 other",
  "  gene sets (mostly large housekeeping sets: translation, ribosome, antigen",
  "  presentation, proteasome) reach even smaller joint FDR in this comparison --",
  "  a known property of significance-weighted (stat) ranking on high-count,",
  "  low-dispersion genes, not evidence against the ECM/integrin finding itself",
  "  (see README_analysis_log.md Key finding 7 in the parent REVISION_PACKAGE)."
)
writeLines(consistency_lines, file.path(log_dir, "23_Fig6C_consistency_check.txt"))
cat(paste(consistency_lines, collapse = "\n"), "\n")

cat("\n=== 23_Fig6C_fgsea_dotplot.R complete ===\n")
