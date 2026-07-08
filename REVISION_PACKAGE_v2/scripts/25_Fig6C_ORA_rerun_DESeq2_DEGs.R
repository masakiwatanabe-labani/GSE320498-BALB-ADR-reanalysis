#!/usr/bin/env Rscript
# Fig. 6C: re-run the manuscript's stated method (Reactome pathway
# over-representation analysis, ORA) with the DEG list swapped from the
# original edgeR-derived list to the canonical DESeq2-derived list. This is
# the ONE genuinely new analysis in this task -- everything else in this
# reviewer-response task reuses existing DESeq2/fgsea output.
#
# Input DEG list: DE_ADR_B_vs_A.tsv (Day 5, B-ADR vs A-ADR, A-ADR1 excluded,
# main analysis), padj < 0.05 (manuscript's stated significance criterion;
# no additional fold-change filter applied -- noted explicitly below).
# Background/universe: all genes tested in that same comparison (19,662).
# Tool: ReactomePA::enrichPathway() (organism = mouse), BH correction.

suppressMessages({
  library(clusterProfiler)
  library(ReactomePA)
  library(org.Mm.eg.db)
  library(ggplot2)
})

set.seed(20260220)
indir <- "/usr/local/jupyter/ADR_BALB"
outdir_main <- file.path(indir, "outputs")
outdir_v2 <- file.path(outdir_main, "REVISION_PACKAGE/REVISION_PACKAGE_v2")
fig_dir <- file.path(outdir_v2, "figures")
tab_dir <- file.path(outdir_v2, "tables")
log_dir <- file.path(outdir_v2, "logs")

d5 <- read.delim(file.path(outdir_main, "tables/DE_ADR_B_vs_A.tsv"), stringsAsFactors = FALSE)

deg <- d5$gene[!is.na(d5$padj) & d5$padj < 0.05]
universe <- d5$gene
cat(sprintf("Input DEG list: padj<0.05, D5 (B-ADR vs A-ADR, A-ADR1 excluded, main analysis), no fold-change filter applied. n DEGs = %d\n", length(deg)))
cat(sprintf("Background/universe: all tested genes in this comparison. n universe = %d\n", length(universe)))

# ---- gene symbol -> Entrez ID mapping (mouse) ----
map_deg <- suppressWarnings(bitr(deg, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db))
map_universe <- suppressWarnings(bitr(universe, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db))

n_deg_mapped <- length(unique(map_deg$ENTREZID))
n_universe_mapped <- length(unique(map_universe$ENTREZID))
cat(sprintf("DEG symbols mapped to Entrez ID: %d / %d (%.1f%%)\n", n_deg_mapped, length(deg), 100 * n_deg_mapped / length(deg)))
cat(sprintf("Universe symbols mapped to Entrez ID: %d / %d (%.1f%%)\n", n_universe_mapped, length(universe), 100 * n_universe_mapped / length(universe)))

# ---- ORA via ReactomePA::enrichPathway ----
ora <- enrichPathway(
  gene = unique(map_deg$ENTREZID),
  universe = unique(map_universe$ENTREZID),
  organism = "mouse",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  pAdjustMethod = "BH",
  minGSSize = 5,
  maxGSSize = 500,
  readable = TRUE
)
ora_df <- as.data.frame(ora)
n_terms_tested <- nrow(ora_df)
n_sig_terms <- sum(ora_df$p.adjust < 0.05, na.rm = TRUE)

cat(sprintf("\nORA (ReactomePA::enrichPathway, organism=mouse): %d Reactome terms tested, %d significant at BH-adjusted p<0.05\n",
            n_terms_tested, n_sig_terms))

ora_df <- ora_df[order(ora_df$p.adjust), ]
write.csv(ora_df, file.path(tab_dir, "TableS_ORA_D5_full.csv"), row.names = FALSE)
cat("Saved:", file.path(tab_dir, "TableS_ORA_D5_full.csv"), "\n")

top20 <- head(ora_df, 20)
cat("\n=== Top 20 ORA terms (BH p.adjust ranked) ===\n")
print(top20[, c("Description","GeneRatio","BgRatio","pvalue","p.adjust","Count")], row.names = FALSE)

# ---- ECM/collagen check ----
ecm_pattern <- "ECM|COLLAGEN|EXTRACELLULAR MATRIX|EXTRACELLULAR_MATRIX|INTEGRIN|CELL ADHESION|LAMININ|BASEMENT MEMBRANE"
ecm_terms <- ora_df[grepl(ecm_pattern, ora_df$Description, ignore.case = TRUE), ]
cat(sprintf("\nECM/collagen/integrin/adhesion-related terms found in ORA output: %d\n", nrow(ecm_terms)))
if (nrow(ecm_terms) > 0) print(ecm_terms[, c("Description","GeneRatio","pvalue","p.adjust","Count")], row.names = FALSE)

# ---- text summary for manuscript ----
if (nrow(ora_df) > 0) {
  top1 <- ora_df[1, ]
  summary_sentence <- sprintf(
    "Top term (ORA, DESeq2-derived DEG list, D5 B-ADR vs A-ADR, A-ADR1 excluded): %s; %s genes (%s); BH-adjusted p = %s.",
    top1$Description, top1$Count, top1$GeneRatio, format(top1$p.adjust, digits = 3)
  )
} else {
  summary_sentence <- "No ORA terms returned -- 該当なし."
}
cat("\n", summary_sentence, "\n", sep = "")

# ---- dot plot: x=GeneRatio, size=Count, colour=p.adjust (same semantics as original Fig 6C) ----
plot_df <- top20
parse_ratio <- function(x) sapply(strsplit(x, "/"), function(v) as.numeric(v[1]) / as.numeric(v[2]))
plot_df$GeneRatioNum <- parse_ratio(plot_df$GeneRatio)
wrap_label <- function(x) vapply(x, function(s) paste(strwrap(s, width = 42), collapse = "\n"), character(1))
plot_df$DescriptionWrapped <- wrap_label(plot_df$Description)
plot_df$DescriptionWrapped <- factor(plot_df$DescriptionWrapped,
                                      levels = rev(plot_df$DescriptionWrapped[order(plot_df$GeneRatioNum)]))

p <- ggplot(plot_df, aes(x = GeneRatioNum, y = DescriptionWrapped)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_viridis_c(option = "viridis", direction = -1, name = "BH-adjusted\np-value") +
  scale_size_continuous(name = "Gene\ncount", range = c(2.2, 9)) +
  labs(
    title = "Figure 6C (ORA): Reactome pathway over-representation,\nDay 5 post-ADR DEGs",
    subtitle = paste0(
      "Input: DESeq2-derived DEGs (padj<0.05, D5 B-ADR vs A-ADR, A-ADR1 excluded, main analysis), n=", length(deg), ".\n",
      "Universe: all ", length(universe), " genes tested in this comparison. ReactomePA::enrichPathway (organism=mouse), BH-adjusted.\n",
      "Top 20 of ", n_terms_tested, " terms tested by adjusted p-value."
    ),
    x = "GeneRatio", y = NULL
  ) +
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

h <- max(5, 1.6 + nrow(plot_df) * 0.42)
for (ext in c("png", "pdf")) {
  fp <- file.path(fig_dir, paste0("Figure6C_ORA_dotplot.", ext))
  ggsave(fp, plot = p, width = 12.5, height = h, dpi = 300)
  cat("Saved:", fp, "\n")
}

# ---- versions / seed log ----
version_log <- c(
  "=== Fig. 6C ORA re-run: package versions and settings ===",
  sprintf("ReactomePA: %s", as.character(packageVersion("ReactomePA"))),
  sprintf("clusterProfiler: %s", as.character(packageVersion("clusterProfiler"))),
  sprintf("org.Mm.eg.db: %s", as.character(packageVersion("org.Mm.eg.db"))),
  sprintf("AnnotationDbi: %s", as.character(packageVersion("AnnotationDbi"))),
  sprintf("reactome.db: %s", as.character(tryCatch(packageVersion("reactome.db"), error = function(e) NA))),
  "Random seed: 20260220",
  "",
  "Input DEG list: DE_ADR_B_vs_A.tsv, padj<0.05, no additional fold-change filter, A-ADR1 excluded (main analysis).",
  sprintf("n DEGs (gene symbols): %d; mapped to Entrez ID: %d (%.1f%%)", length(deg), n_deg_mapped, 100*n_deg_mapped/length(deg)),
  sprintf("Background/universe: all tested genes in DE_ADR_B_vs_A.tsv, n = %d; mapped to Entrez ID: %d (%.1f%%)",
          length(universe), n_universe_mapped, 100*n_universe_mapped/length(universe)),
  "enrichPathway(organism='mouse', pAdjustMethod='BH', minGSSize=5, maxGSSize=500, pvalueCutoff=1, qvalueCutoff=1) -- all terms returned, filtered post hoc.",
  "",
  sprintf("n Reactome terms tested: %d", n_terms_tested),
  sprintf("n significant terms (BH p.adjust<0.05): %d", n_sig_terms),
  sprintf("n ECM/collagen/integrin/adhesion-name-matched terms: %d", nrow(ecm_terms)),
  "",
  summary_sentence
)
writeLines(version_log, file.path(log_dir, "25_Fig6C_ORA_versions_and_summary.txt"))

cat("\n=== 25_Fig6C_ORA_rerun_DESeq2_DEGs.R complete ===\n")
