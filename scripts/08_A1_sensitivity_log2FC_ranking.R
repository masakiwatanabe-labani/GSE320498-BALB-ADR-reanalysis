#!/usr/bin/env Rscript
# A1 in/out sensitivity, re-run with the corrected ranking metric (signed
# log2FC, matching manuscript Methods) established in
# 07_gsea_log2FC_ranking_primary.R. Only ADR_B_vs_A and A_ADR_vs_Ctrl can be
# affected by A-ADR1 inclusion/exclusion (A1 is a Day-5 AJcl-ADR sample).

suppressMessages({
  library(DESeq2)
  library(fgsea)
  library(msigdbr)
})

set.seed(20260220)
indir <- "/usr/local/jupyter/ADR_BALB"
outdir <- file.path(indir, "outputs")

dds_group_main <- readRDS(file.path(outdir, "tables/dds_group_main.rds"))
dds_group_sens <- readRDS(file.path(outdir, "tables/dds_group_sens.rds"))

get_de <- function(dds, contrast) {
  res <- results(dds, contrast = contrast, alpha = 0.05)
  df <- as.data.frame(res)
  df$gene <- rownames(df)
  df[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
}

sens_comparisons <- list(
  ADR_B_vs_A    = c("group", "B_ADR", "A_ADR"),
  A_ADR_vs_Ctrl = c("group", "A_ADR", "A_Ctrl")
)
de_main <- lapply(sens_comparisons, function(ctr) get_de(dds_group_main, ctr))
de_sens <- lapply(sens_comparisons, function(ctr) get_de(dds_group_sens, ctr))

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
joint_list <- c(reactome_list,
                list(TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING = podo_ageing,
                     GOBP_GLOMERULAR_EPITHELIAL_CELL_DIFFERENTIATION = unique(podo_identity),
                     KARAISKOS2018_PODOCYTE_EXCLUSIVE = unique(karaiskos_exclusive),
                     KARAISKOS2018_PODOCYTE_TOP50 = unique(karaiskos_top50)))

focal_sets <- c("TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING", "REACTOME_INTEGRIN_SIGNALING",
                 "GOBP_GLOMERULAR_EPITHELIAL_CELL_DIFFERENTIATION",
                 "KARAISKOS2018_PODOCYTE_EXCLUSIVE", "KARAISKOS2018_PODOCYTE_TOP50")

run_gsea_log2fc <- function(df) {
  df <- df[!is.na(df$log2FoldChange), ]
  ranks <- sort(setNames(df$log2FoldChange, df$gene), decreasing = TRUE)
  set.seed(20260220)
  fgsea(pathways = joint_list, stats = ranks, eps = 0, minSize = 5, maxSize = 500)
}

rows <- list()
for (nm in names(sens_comparisons)) {
  rm_ <- run_gsea_log2fc(de_main[[nm]])
  rs_ <- run_gsea_log2fc(de_sens[[nm]])
  for (pw in focal_sets) {
    rows[[paste(nm, pw, "excl")]] <- data.frame(comparison = nm, pathway = pw, dataset = "A1_excluded_main",
      NES = rm_$NES[rm_$pathway == pw], pval = rm_$pval[rm_$pathway == pw], FDR = rm_$padj[rm_$pathway == pw])
    rows[[paste(nm, pw, "incl")]] <- data.frame(comparison = nm, pathway = pw, dataset = "A1_included_sens",
      NES = rs_$NES[rs_$pathway == pw], pval = rs_$pval[rs_$pathway == pw], FDR = rs_$padj[rs_$pathway == pw])
  }
}
summary_df <- do.call(rbind, rows)
rownames(summary_df) <- NULL
write.table(summary_df, file.path(outdir, "tables/GSEA_log2FCrank_A1_sensitivity.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("=== A1 sensitivity, log2FC ranking (matches manuscript Methods) ===\n")
print(summary_df, row.names = FALSE, digits = 4)
