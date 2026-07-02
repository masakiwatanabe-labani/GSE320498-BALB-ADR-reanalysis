#!/usr/bin/env Rscript
# Fig.6B replacement panel: full data package for KARAISKOS2018_PODOCYTE_TOP50
# (proposed main set) shown together with KARAISKOS2018_PODOCYTE_EXCLUSIVE
# (proposed adjacent/confirmatory set), across the two comparisons requested:
#   - ADR_B_vs_A (A1-excluded, main manuscript setting; Fig.6-equivalent)
#   - B_ADR_vs_Ctrl (within-ByJcl ADR response; strongest single result in the
#     whole re-analysis, joint FDR=2.1e-10 for TOP50)
# All numbers use the CANONICAL ranking (DESeq2 Wald stat), consistent with
# 09_gsea_stat_ranking_canonical.R / 10_export_enrichment_curves.R.
#
# Outputs (peer-review-defensibility package):
#   (a) outputs/genesets/KARAISKOS2018_PODOCYTE.gmt -- the gene sets themselves
#   (b) tables/KARAISKOS_TOP50_leading_edge_genes.tsv -- leading-edge gene list
#       per comparison, flagged for well-known canonical podocyte markers
#   (c) tables/KARAISKOS_TOP50_vs_EXCLUSIVE_comparison.tsv -- NES/p/FDR side by
#       side for both sets x both comparisons, one table
#   (d) tables/Supplementary_ageing_set_all_configs.tsv -- ageing set kept in
#       full across every comparison/A1 configuration tested (finding 2/3)
#   plus curve/hit/ranking-vector TSVs for the 4 (set x comparison) enrichment
#   plots, under tables/enrichment_curves/

suppressMessages({
  library(fgsea)
  library(msigdbr)
})

set.seed(20260220)
indir <- "/usr/local/jupyter/ADR_BALB"
outdir <- file.path(indir, "outputs")
curvedir <- file.path(outdir, "tables/enrichment_curves")
gsdir <- file.path(outdir, "genesets")
dir.create(curvedir, showWarnings = FALSE, recursive = TRUE)
dir.create(gsdir, showWarnings = FALSE, recursive = TRUE)

## ---------- gene sets (identical definitions to 07/08/09/10) ----------
karaiskos_exclusive <- c("Nphs2","Cdkn1c","Tcf21","Enpep","Nphs1","Synpo",
                          "Npnt","Wt1","Pard3b","Ptpro","Iqgap2","Mafb")
karaiskos_top50 <- c("Nphs2","Cdkn1c","Clic3","Nupr1","Dpp4","Enpep","Tcf21",
                      "Nphs1","Gadd45a","Rab3b","Rhpn1","Tmsb4x","Col4a3",
                      "Rasl11a","Mafb","Npnt","Arhgap24","Adm","Pak1","Synpo",
                      "Foxd2os","Golim4","Igfbp7","Vegfa","Cd59a","Sdc4",
                      "Sema3g","Tdrd5","Nap1l1","Shisa3","Eif3m","Thsd7a",
                      "Pth1r","Sept11","Ctsl","Podxl","Cryab","Mertk","Htra1",
                      "Nes","Wt1","Npr3","Ildr2","Robo2","Pard3b","Tmem150c",
                      "Gas1","Hoxc8","Iqgap2","Sema3e")
stopifnot(length(karaiskos_exclusive) == 12, length(unique(karaiskos_exclusive)) == 12)
stopifnot(length(karaiskos_top50) == 50, length(unique(karaiskos_top50)) == 50)

## ---------- (a) write GMT ----------
gmt_lines <- c(
  paste(c("KARAISKOS2018_PODOCYTE_EXCLUSIVE",
          "Karaiskos_et_al_2018_JASN_29:2060-2068_GSE111107_SupplTable2_cell-type-exclusive_podocyte_markers_n12",
          karaiskos_exclusive), collapse = "\t"),
  paste(c("KARAISKOS2018_PODOCYTE_TOP50",
          "Karaiskos_et_al_2018_JASN_29:2060-2068_GSE111107_SupplTable3_FindAllMarkers_podocyte_cluster_top50_by_avg_logFC",
          karaiskos_top50), collapse = "\t")
)
writeLines(gmt_lines, file.path(gsdir, "KARAISKOS2018_PODOCYTE.gmt"))

## ---------- canonical stat ranking, per comparison ----------
rank_from_file <- function(fn) {
  df <- read.delim(file.path(outdir, "tables", fn), stringsAsFactors = FALSE)
  df <- df[!is.na(df$stat), ]
  sort(setNames(df$stat, df$gene), decreasing = TRUE)
}
ranks_list <- list(
  ADR_B_vs_A    = rank_from_file("DE_ADR_B_vs_A.tsv"),
  B_ADR_vs_Ctrl = rank_from_file("DE_B_ADR_vs_Ctrl.tsv")
)

pathway_defs <- list(
  KARAISKOS_TOP50     = unique(karaiskos_top50),
  KARAISKOS_EXCLUSIVE = unique(karaiskos_exclusive)
)

## ---------- reimplementation of fgsea::plotEnrichment's internal math ----------
export_curve <- function(pathway_genes, stats, gseaParam = 1) {
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
  statsAdj <- statsAdj / max(abs(statsAdj))
  pathway_idx <- unname(as.vector(na.omit(match(pathway_genes, names(statsAdj)))))
  pathway_idx <- sort(pathway_idx)
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway_idx, returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway_idx - 1, pathway_idx))
  ys <- as.vector(rbind(bottoms, tops))
  curve <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  list(curve = curve, hit_positions = pathway_idx, gene_order = names(statsAdj))
}

## known canonical podocyte identity/injury marker genes for annotation purposes
known_markers <- c("Nphs1","Nphs2","Wt1","Podxl","Synpo","Ptpro","Tcf21","Mafb",
                    "Cdkn1c","Ctsl","Mertk","Npnt")

comparison_rows <- list()
leading_edge_rows <- list()

for (comp in names(ranks_list)) {
  ranks <- ranks_list[[comp]]

  set.seed(20260220)
  reac <- msigdbr(db_species = "MM", species = "mouse", collection = "M2", subcollection = "CP:REACTOME")
  reactome_list <- split(reac$gene_symbol, reac$gs_name)
  m8 <- msigdbr(db_species = "MM", species = "mouse", collection = "M8")
  podo_ageing <- unique(m8$gene_symbol[m8$gs_name == "TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING"])
  joint_list <- c(reactome_list,
                   list(TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING = podo_ageing,
                        KARAISKOS2018_PODOCYTE_EXCLUSIVE = unique(karaiskos_exclusive),
                        KARAISKOS2018_PODOCYTE_TOP50 = unique(karaiskos_top50)))
  set.seed(20260220)
  res_joint <- fgsea(pathways = joint_list, stats = ranks, eps = 0, minSize = 5, maxSize = 500)

  write.table(data.frame(rank_position = seq_along(ranks), gene = names(ranks), stat = as.numeric(ranks)),
              file.path(curvedir, sprintf("ranking_vector_%s_stat.tsv", comp)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  for (set_key in names(pathway_defs)) {
    genes <- pathway_defs[[set_key]]
    res <- export_curve(genes, ranks)
    write.table(res$curve, file.path(curvedir, sprintf("curve_%s_%s.tsv", set_key, comp)),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(data.frame(hit_position = res$hit_positions),
                file.path(curvedir, sprintf("hits_%s_%s.tsv", set_key, comp)),
                sep = "\t", quote = FALSE, row.names = FALSE)

    pw_name <- paste0("KARAISKOS2018_PODOCYTE_", sub("KARAISKOS_", "", set_key))
    r <- res_joint[res_joint$pathway == pw_name, ]
    comparison_rows[[paste(comp, set_key)]] <- data.frame(
      comparison = comp, gene_set = pw_name, n_genes = r$size,
      NES = r$NES, nominal_p = r$pval, joint_FDR = r$padj,
      leading_edge_n = lengths(r$leadingEdge)
    )

    if (set_key == "KARAISKOS_TOP50") {
      le_genes <- r$leadingEdge[[1]]
      leading_edge_rows[[comp]] <- data.frame(
        comparison = comp, gene_set = pw_name, leading_edge_gene = le_genes,
        known_canonical_podocyte_marker = le_genes %in% known_markers
      )
    }
  }
}

comparison_df <- do.call(rbind, comparison_rows)
rownames(comparison_df) <- NULL
write.table(comparison_df, file.path(outdir, "tables/KARAISKOS_TOP50_vs_EXCLUSIVE_comparison.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("=== (c) TOP50 vs EXCLUSIVE, both comparisons ===\n")
print(comparison_df, row.names = FALSE, digits = 4)
cat("\n")

leading_edge_df <- do.call(rbind, leading_edge_rows)
rownames(leading_edge_df) <- NULL
write.table(leading_edge_df, file.path(outdir, "tables/KARAISKOS_TOP50_leading_edge_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("=== (b) TOP50 leading-edge genes (known markers flagged) ===\n")
print(leading_edge_df, row.names = FALSE)
cat("\n")

## ---------- (d) ageing set, full A1 in/out, all comparisons ----------
judgment <- read.delim(file.path(outdir, "tables/GSEA_canonical_focal_judgment_table.tsv"), stringsAsFactors = FALSE)
ageing_all <- judgment[judgment$pathway == "TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING", ]
write.table(ageing_all, file.path(outdir, "tables/Supplementary_ageing_set_all_configs.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("=== (d) Podocyte-ageing set, all comparisons/A1 configs (supplementary) ===\n")
print(ageing_all, row.names = FALSE, digits = 4)

cat("\nWrote GMT, leading-edge table, comparison table, ageing supplementary table, and curve/hit/ranking TSVs.\n")
