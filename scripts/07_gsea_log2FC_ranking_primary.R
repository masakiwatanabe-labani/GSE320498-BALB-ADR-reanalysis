#!/usr/bin/env Rscript
# CORRECTION run: the manuscript's own Methods states GSEA was run on genes
# "ordered by signed log2 fold change" (see manuscript text, Differential
# Expression and Pathway Analysis paragraph). Every prior script in this
# project (05_gsea.R, 05b, 05d, 05e, 06) instead ranked by the DESeq2 Wald
# "stat" column. This script re-runs preranked GSEA using signed log2FC as
# the ranking statistic, matching the manuscript's stated method, for all
# gene sets examined so far: Reactome (joint, incl. REACTOME_INTEGRIN_SIGNALING),
# TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING (as in the manuscript, joint + single),
# GOBP_GLOMERULAR_EPITHELIAL_CELL_DIFFERENTIATION, and the two Karaiskos et al.
# 2018 (GSE111107) podocyte marker sets. This is now the PRIMARY ranking
# convention for this project; the earlier stat-ranked numbers are retained
# in their original files as a secondary/alternative-ranking sensitivity check
# (they used a variance-aware statistic, which is a defensible alternative,
# just not what the manuscript itself used).

suppressMessages({
  library(fgsea)
  library(msigdbr)
})

set.seed(20260220)
indir <- "/usr/local/jupyter/ADR_BALB"
outdir <- file.path(indir, "outputs")

reac <- msigdbr(db_species = "MM", species = "mouse", collection = "M2", subcollection = "CP:REACTOME")
reactome_list <- split(reac$gene_symbol, reac$gs_name)

m8 <- msigdbr(db_species = "MM", species = "mouse", collection = "M8")
podo_ageing <- unique(m8$gene_symbol[m8$gs_name == "TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING"])

go <- msigdbr(db_species = "MM", species = "mouse", collection = "M5")
podo_identity <- unique(go$gene_symbol[go$gs_name == "GOBP_GLOMERULAR_EPITHELIAL_CELL_DIFFERENTIATION"])

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
  GOBP_GLOMERULAR_EPITHELIAL_CELL_DIFFERENTIATION = unique(podo_identity),
  KARAISKOS2018_PODOCYTE_EXCLUSIVE = unique(karaiskos_exclusive),
  KARAISKOS2018_PODOCYTE_TOP50 = unique(karaiskos_top50)
)
joint_list <- c(reactome_list, custom_sets)

comparisons <- c(
  baseline_B_vs_A = "DE_baseline_B_vs_A.tsv",
  ADR_B_vs_A = "DE_ADR_B_vs_A.tsv",
  A_ADR_vs_Ctrl = "DE_A_ADR_vs_Ctrl.tsv",
  B_ADR_vs_Ctrl = "DE_B_ADR_vs_Ctrl.tsv"
)

focal_sets <- c("TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING", "REACTOME_INTEGRIN_SIGNALING",
                 "GOBP_GLOMERULAR_EPITHELIAL_CELL_DIFFERENTIATION",
                 "KARAISKOS2018_PODOCYTE_EXCLUSIVE", "KARAISKOS2018_PODOCYTE_TOP50")

rows <- list()
for (nm in names(comparisons)) {
  df <- read.delim(file.path(outdir, "tables", comparisons[[nm]]), stringsAsFactors = FALSE)
  df <- df[!is.na(df$log2FC), ]
  ranks <- sort(setNames(df$log2FC, df$gene), decreasing = TRUE)

  set.seed(20260220)
  res_joint <- fgsea(pathways = joint_list, stats = ranks, eps = 0, minSize = 5, maxSize = 500)

  # single-set podocyte-ageing test, mirroring the manuscript's apparent
  # single-set approach (nominal P == FDR q in Fig. 6B)
  set.seed(20260220)
  res_podo_single <- fgsea(pathways = list(TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING = podo_ageing),
                            stats = ranks, eps = 0, minSize = 5, maxSize = 500)

  res_joint_out <- res_joint[order(res_joint$pval), ]
  res_joint_out$leadingEdge <- sapply(res_joint_out$leadingEdge, paste, collapse = ";")
  write.table(res_joint_out,
              file.path(outdir, sprintf("tables/GSEA_%s_full_joint_log2FCrank.tsv", nm)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  for (pw in focal_sets) {
    rj <- res_joint[res_joint$pathway == pw, ]
    rows[[paste(nm, pw)]] <- data.frame(
      comparison = nm, pathway = pw,
      NES_joint = rj$NES, pval_joint = rj$pval, FDR_joint = rj$padj, size = rj$size
    )
  }
  rows[[paste(nm, "PODO_AGEING_SINGLE")]] <- data.frame(
    comparison = nm, pathway = "TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING_single",
    NES_joint = res_podo_single$NES, pval_joint = res_podo_single$pval,
    FDR_joint = res_podo_single$padj, size = res_podo_single$size
  )
}

summary_df <- do.call(rbind, rows)
rownames(summary_df) <- NULL
write.table(summary_df, file.path(outdir, "tables/GSEA_log2FCrank_summary_all_sets.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("=== GSEA re-run with signed log2FC ranking (matches manuscript Methods) ===\n")
print(summary_df, row.names = FALSE, digits = 4)
