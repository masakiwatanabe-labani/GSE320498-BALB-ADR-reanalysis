#!/usr/bin/env Rscript
# Step 4 follow-up: the MSigDB M8 "TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING"
# set showed a sign discrepancy vs the manuscript's reported NES for the
# Day-5 ADR ByJcl-vs-AJcl comparison (see 05_gsea.R / logs/05_gsea_versions.txt
# and outputs/figures/CHECK_enrichment_podocyte_ADR_B_vs_A.png). That set is
# derived from an AGEING comparison within podocytes (leading-edge genes are
# MHC/immune/senescence genes, not canonical podocyte-identity genes), so it
# may not be capturing what the manuscript intended by "podocyte program".
#
# As a follow-up (requested by user), we additionally test a pre-existing,
# independently curated GO Biological Process gene set that IS explicitly a
# canonical podocyte-identity signature:
#   GOBP_GLOMERULAR_EPITHELIAL_CELL_DIFFERENTIATION (MSigDB M5:GO:BP, n=24
#   genes; includes Nphs1, Nphs2, Podxl, Wt1, Cd2ap, Magi2, Ptpro, Foxc2, etc.)
# This gene set was NOT used in the original manuscript and is added here
# purely as a sensitivity/robustness check post-hoc, selected for being a
# standard, pre-existing curated GO term matching "podocyte identity" rather
# than one hand-picked from many candidates to force a particular answer.
# This is disclosed explicitly as an EXPLORATORY addition, not a replacement
# for the original pre-registered analysis.

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
cat("GOBP_GLOMERULAR_EPITHELIAL_CELL_DIFFERENTIATION genes:", paste(sort(podo_identity), collapse=", "), "\n\n")

joint_list <- c(reactome_list,
                 list(TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING = podo_ageing,
                      GOBP_GLOMERULAR_EPITHELIAL_CELL_DIFFERENTIATION = podo_identity))

comparisons <- c(
  baseline_B_vs_A = "DE_baseline_B_vs_A.tsv",
  ADR_B_vs_A = "DE_ADR_B_vs_A.tsv",
  A_ADR_vs_Ctrl = "DE_A_ADR_vs_Ctrl.tsv",
  B_ADR_vs_Ctrl = "DE_B_ADR_vs_Ctrl.tsv"
)

rows <- list()
for (nm in names(comparisons)) {
  df <- read.delim(file.path(outdir, "tables", comparisons[[nm]]), stringsAsFactors = FALSE)
  df <- df[!is.na(df$stat), ]
  ranks <- sort(setNames(df$stat, df$gene), decreasing = TRUE)

  set.seed(20260220)
  res <- fgsea(pathways = joint_list, stats = ranks, eps = 0, minSize = 5, maxSize = 500)
  r_identity <- res[res$pathway == "GOBP_GLOMERULAR_EPITHELIAL_CELL_DIFFERENTIATION", ]
  r_ageing   <- res[res$pathway == "TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING", ]

  rows[[nm]] <- data.frame(
    comparison = nm,
    n_genesets_tested_joint = nrow(res),
    identity_NES = r_identity$NES, identity_pval = r_identity$pval, identity_FDR = r_identity$padj,
    identity_size = r_identity$size,
    ageing_NES = r_ageing$NES, ageing_pval = r_ageing$pval, ageing_FDR = r_ageing$padj
  )
}

summary_df <- do.call(rbind, rows)
write.table(summary_df, file.path(outdir, "tables/GSEA_alt_podocyte_identity_vs_ageing.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n=== GOBP_GLOMERULAR_EPITHELIAL_CELL_DIFFERENTIATION (podocyte IDENTITY, exploratory) vs TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING (original) ===\n")
print(summary_df, row.names = FALSE)
