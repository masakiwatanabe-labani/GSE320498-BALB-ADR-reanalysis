#!/usr/bin/env Rscript
# Follow-up to 05b_alt_podocyte_geneset.R: A-ADR1 in/out sensitivity check for
# the podocyte-IDENTITY gene set (GOBP_GLOMERULAR_EPITHELIAL_CELL_DIFFERENTIATION,
# n=24: Nphs1/Nphs2/Podxl/Wt1/Cd2ap/Magi2/Ptpro/Foxc2/etc.), requested by the
# user as a completion of the A1-sensitivity work already done in
# 06_A1_sensitivity.R (which only covered the TABULA_MURIS ageing set and, at
# the GSEA level only, ADR_B_vs_A for this identity set). This script adds:
#   (1) GSEA A1 in/out for ALL FOUR comparisons (06 only covered ADR_B_vs_A)
#   (2) gene-level log2FC/padj A1 in/out for all 24 identity genes (06 only
#       tracked Wt1/Nphs1 among its "key_genes" list)

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

# main = A1 excluded (matches manuscript's stated n=2/n=3); sens = A1 included
comparisons <- list(
  baseline_B_vs_A = c("group", "B_Ctrl", "A_Ctrl"),
  ADR_B_vs_A      = c("group", "B_ADR", "A_ADR"),
  A_ADR_vs_Ctrl   = c("group", "A_ADR", "A_Ctrl"),
  B_ADR_vs_Ctrl   = c("group", "B_ADR", "B_Ctrl")
)

de_main <- lapply(comparisons, function(ctr) get_de(dds_group_main, ctr))
de_sens <- lapply(comparisons, function(ctr) get_de(dds_group_sens, ctr))

go <- msigdbr(db_species = "MM", species = "mouse", collection = "M5")
podo_identity <- unique(go$gene_symbol[go$gs_name == "GOBP_GLOMERULAR_EPITHELIAL_CELL_DIFFERENTIATION"])
identity_list <- list(GOBP_GLOMERULAR_EPITHELIAL_CELL_DIFFERENTIATION = podo_identity)

run_gsea <- function(df) {
  df <- df[!is.na(df$stat), ]
  ranks <- sort(setNames(df$stat, df$gene), decreasing = TRUE)
  set.seed(20260220)
  fgsea(pathways = identity_list, stats = ranks, eps = 0, minSize = 5, maxSize = 500)
}

gsea_rows <- list()
for (nm in names(comparisons)) {
  rm_ <- run_gsea(de_main[[nm]]); rs_ <- run_gsea(de_sens[[nm]])
  gsea_rows[[paste0(nm, "_A1excluded")]] <- data.frame(comparison = nm, dataset = "A1_excluded_main",
                                                        NES = rm_$NES, pval = rm_$pval, padj = rm_$padj)
  gsea_rows[[paste0(nm, "_A1included")]] <- data.frame(comparison = nm, dataset = "A1_included_sens",
                                                        NES = rs_$NES, pval = rs_$pval, padj = rs_$padj)
}
gsea_summary <- do.call(rbind, gsea_rows)
rownames(gsea_summary) <- NULL

cat("=== Podocyte-IDENTITY set (GOBP_GLOMERULAR_EPITHELIAL_CELL_DIFFERENTIATION): A1 sensitivity across all 4 comparisons ===\n")
print(gsea_summary, row.names = FALSE)

# gene-level table for all 24 identity genes, A1 in vs out, for the two
# comparisons where A1 membership can actually matter (A-ADR1 is only in the
# Day-5 AJcl ADR group)
gene_level <- function(nm) {
  m <- merge(de_main[[nm]][de_main[[nm]]$gene %in% podo_identity, c("gene", "log2FoldChange", "padj")],
             de_sens[[nm]][de_sens[[nm]]$gene %in% podo_identity, c("gene", "log2FoldChange", "padj")],
             by = "gene", suffixes = c("_A1excl", "_A1incl"))
  m$comparison <- nm
  m[order(m$padj_A1excl), c("comparison", "gene", "log2FoldChange_A1excl", "padj_A1excl",
                             "log2FoldChange_A1incl", "padj_A1incl")]
}
gene_level_table <- rbind(gene_level("ADR_B_vs_A"), gene_level("A_ADR_vs_Ctrl"))

cat("\n=== Gene-level (24 identity genes), A1 excluded vs included ===\n")
print(gene_level_table, row.names = FALSE, digits = 3)

write.table(gsea_summary, file.path(outdir, "tables/GSEA_alt_podocyte_identity_A1_sensitivity.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gene_level_table, file.path(outdir, "tables/DE_alt_podocyte_identity_genes_A1_sensitivity.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nSaved outputs/tables/GSEA_alt_podocyte_identity_A1_sensitivity.tsv and DE_alt_podocyte_identity_genes_A1_sensitivity.tsv\n")
