#!/usr/bin/env Rscript
# Export fgsea running-enrichment-score curve data (identical algorithm to
# fgsea::plotEnrichment, reimplemented here to get raw numbers instead of a
# ggplot object), gene-set hit positions, and the ranking vector, for the 3
# Fig.6-style panels requested:
#   1) REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION (positive NES, Fig.6C theme)
#   2) REACTOME_INTEGRIN_SIGNALING (positive NES, Fig.6D)
#   3) KARAISKOS2018_PODOCYTE_TOP50 (negative NES; chosen over EXCLUSIVE for a
#      denser leading-edge barcode -- 49 vs 12 genes -- while still highly
#      significant; see README Key finding 4)
# All three use the CANONICAL ranking (DESeq2 Wald stat) and the comparison
# ADR_B_vs_A, A1-excluded (main manuscript setting) -- same run convention as
# 09_gsea_stat_ranking_canonical.R, so NES/FDR quoted on the plots are from
# the identical run, not recomputed separately.

suppressMessages({
  library(fgsea)
  library(msigdbr)
})

set.seed(20260220)
indir <- "/usr/local/jupyter/ADR_BALB"
outdir <- file.path(indir, "outputs")
curvedir <- file.path(outdir, "tables/enrichment_curves")
dir.create(curvedir, showWarnings = FALSE, recursive = TRUE)

## ---------- build the same joint universe as script 09 ----------
reac <- msigdbr(db_species = "MM", species = "mouse", collection = "M2", subcollection = "CP:REACTOME")
reactome_list <- split(reac$gene_symbol, reac$gs_name)

m8 <- msigdbr(db_species = "MM", species = "mouse", collection = "M8")
podo_ageing <- unique(m8$gene_symbol[m8$gs_name == "TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING"])

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
custom_sets <- list(
  TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING = podo_ageing,
  KARAISKOS2018_PODOCYTE_EXCLUSIVE = unique(karaiskos_exclusive),
  KARAISKOS2018_PODOCYTE_TOP50 = unique(karaiskos_top50)
)
joint_list <- c(reactome_list, custom_sets)

## ---------- canonical ranking: DESeq2 stat, ADR_B_vs_A, A1-excluded main ----------
df <- read.delim(file.path(outdir, "tables/DE_ADR_B_vs_A.tsv"), stringsAsFactors = FALSE)
df <- df[!is.na(df$stat), ]
ranks <- sort(setNames(df$stat, df$gene), decreasing = TRUE)

set.seed(20260220)
res_joint <- fgsea(pathways = joint_list, stats = ranks, eps = 0, minSize = 5, maxSize = 500)

targets <- c(
  ECM_ORGANIZATION = "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
  INTEGRIN_SIGNALING = "REACTOME_INTEGRIN_SIGNALING",
  KARAISKOS_TOP50 = "KARAISKOS2018_PODOCYTE_TOP50"
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

  list(curve = curve, hit_positions = pathway_idx, statsAdj = statsAdj,
       gene_order = names(statsAdj), raw_stat_ordered = as.numeric(stats[ord]))
}

meta_rows <- list()
for (nm in names(targets)) {
  pw <- targets[[nm]]
  genes <- joint_list[[pw]]
  res <- export_curve(genes, ranks)

  write.table(res$curve, file.path(curvedir, sprintf("curve_%s.tsv", nm)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.frame(hit_position = res$hit_positions),
              file.path(curvedir, sprintf("hits_%s.tsv", nm)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  r <- res_joint[res_joint$pathway == pw, ]
  meta_rows[[nm]] <- data.frame(
    panel = nm, pathway = pw, NES = r$NES, pval = r$pval, FDR_joint = r$padj,
    size = r$size, n_genes_total = length(res$gene_order)
  )
}

# ranking vector shared by all 3 panels (same comparison, same ranking)
ranking_df <- data.frame(rank_position = seq_along(ranks), gene = names(ranks), stat = as.numeric(ranks))
write.table(ranking_df, file.path(curvedir, "ranking_vector_ADR_B_vs_A_stat.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

meta_df <- do.call(rbind, meta_rows)
write.table(meta_df, file.path(curvedir, "panel_metadata.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("=== Panel metadata (NES/FDR from the same run as the curves) ===\n")
print(meta_df, row.names = FALSE, digits = 4)
cat("\nWrote curve/hit/ranking data to outputs/tables/enrichment_curves/\n")
