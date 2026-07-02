#!/usr/bin/env Rscript
# Step 5 (R1-7, R2-4): A1 (=A-ADR1) inclusion/exclusion sensitivity analysis.
#
# A-ADR1 only appears in the Day-5 ADR AJcl group. It therefore CANNOT affect
# baseline_B_vs_A (no overlapping samples) -- that comparison is included
# below purely as an explicit negative control demonstrating no change, since
# the task instructions asked for baseline B-vs-A to be part of this
# comparison. The comparisons where A1 inclusion/exclusion can plausibly
# matter are: ADR_B_vs_A (Day-5 substrain comparison, = original Fig.6A),
# A_ADR_vs_Ctrl (within-AJcl ADR response), and the substrain:treatment
# interaction term -- all three are run both ways below, together with the
# corresponding GSEA for ADR_B_vs_A (podocyte ageing set, podocyte-identity
# GO set, and Reactome Integrin signaling).

suppressMessages({
  library(DESeq2)
  library(fgsea)
  library(msigdbr)
  library(openxlsx)
})

set.seed(20260220)
indir <- "/usr/local/jupyter/ADR_BALB"
outdir <- file.path(indir, "outputs")

dds_group_main <- readRDS(file.path(outdir, "tables/dds_group_main.rds"))
dds_group_sens <- readRDS(file.path(outdir, "tables/dds_group_sens.rds"))
dds_int_main   <- readRDS(file.path(outdir, "tables/dds_int_main.rds"))
dds_int_sens   <- readRDS(file.path(outdir, "tables/dds_int_sens.rds"))

get_de <- function(dds, contrast) {
  res <- results(dds, contrast = contrast, alpha = 0.05)
  df <- as.data.frame(res)
  df$gene <- rownames(df)
  df[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
}

# ---- (1) baseline B vs A: negative control, A1 not in this dataset ----
base_main <- get_de(dds_group_main, c("group", "B_Ctrl", "A_Ctrl"))
base_sens <- get_de(dds_group_sens, c("group", "B_Ctrl", "A_Ctrl"))

# ---- (2) ADR_B_vs_A: Day-5 substrain comparison (Fig.6A-equivalent) ----
adr_main <- get_de(dds_group_main, c("group", "B_ADR", "A_ADR"))
adr_sens <- get_de(dds_group_sens, c("group", "B_ADR", "A_ADR"))

# ---- (3) A_ADR_vs_Ctrl: within-AJcl ADR response ----
a_main <- get_de(dds_group_main, c("group", "A_ADR", "A_Ctrl"))
a_sens <- get_de(dds_group_sens, c("group", "A_ADR", "A_Ctrl"))

# ---- (4) interaction term ----
int_main <- as.data.frame(results(dds_int_main, name = "substrainB.treatmentADR"))
int_main$gene <- rownames(int_main)
int_sens <- as.data.frame(results(dds_int_sens, name = "substrainB.treatmentADR"))
int_sens$gene <- rownames(int_sens)

n_sig <- function(df) sum(df$padj < 0.05, na.rm = TRUE)

key_genes <- c("Wt1", "Nphs1", "Tmem215", "Nlrp1b", "Glp1r", "Serpine1", "Col4a1", "Col4a2", "Loxl1")
key_compare <- function(main_df, sens_df, label) {
  m <- merge(main_df[main_df$gene %in% key_genes, c("gene", "log2FoldChange", "padj")],
             sens_df[sens_df$gene %in% key_genes, c("gene", "log2FoldChange", "padj")],
             by = "gene", suffixes = c("_main_A1excl", "_sens_A1incl"))
  m$comparison <- label
  m
}

key_tables <- rbind(
  key_compare(base_main, base_sens, "baseline_B_vs_A"),
  key_compare(adr_main, adr_sens, "ADR_B_vs_A"),
  key_compare(a_main, a_sens, "A_ADR_vs_Ctrl"),
  key_compare(int_main, int_sens, "interaction_substrainxtreatment")
)
key_tables <- key_tables[, c("comparison", "gene", "log2FoldChange_main_A1excl", "padj_main_A1excl",
                              "log2FoldChange_sens_A1incl", "padj_sens_A1incl")]

overview <- data.frame(
  comparison = c("baseline_B_vs_A", "ADR_B_vs_A", "A_ADR_vs_Ctrl", "interaction_substrainxtreatment"),
  n_sig_A1_excluded_main = c(n_sig(base_main), n_sig(adr_main), n_sig(a_main), n_sig(int_main)),
  n_sig_A1_included_sens = c(n_sig(base_sens), n_sig(adr_sens), n_sig(a_sens), n_sig(int_sens))
)

cat("=== Step 5: A1 sensitivity -- significant gene counts (padj<0.05) ===\n")
print(overview, row.names = FALSE)
cat("\n=== Key gene log2FC/padj, A1-excluded (main) vs A1-included (sens) ===\n")
print(key_tables, row.names = FALSE)

# ---- GSEA sensitivity for ADR_B_vs_A (the comparison most affected by A1) ----
reac <- msigdbr(db_species = "MM", species = "mouse", collection = "M2", subcollection = "CP:REACTOME")
reactome_list <- split(reac$gene_symbol, reac$gs_name)
m8 <- msigdbr(db_species = "MM", species = "mouse", collection = "M8")
podo_ageing <- unique(m8$gene_symbol[m8$gs_name == "TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING"])
go <- msigdbr(db_species = "MM", species = "mouse", collection = "M5")
podo_identity <- unique(go$gene_symbol[go$gs_name == "GOBP_GLOMERULAR_EPITHELIAL_CELL_DIFFERENTIATION"])
joint_list <- c(reactome_list,
                 list(TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING = podo_ageing,
                      GOBP_GLOMERULAR_EPITHELIAL_CELL_DIFFERENTIATION = podo_identity))

run_gsea <- function(df) {
  df <- df[!is.na(df$stat), ]
  ranks <- sort(setNames(df$stat, df$gene), decreasing = TRUE)
  set.seed(20260220)
  fgsea(pathways = joint_list, stats = ranks, eps = 0, minSize = 5, maxSize = 500)
}

res_adr_main <- run_gsea(adr_main)
res_adr_sens <- run_gsea(adr_sens)

pull <- function(res, pw) res[res$pathway == pw, c("pathway", "NES", "pval", "padj")]

gsea_sens_compare <- rbind(
  cbind(dataset = "A1_excluded_main", pull(res_adr_main, "TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING")),
  cbind(dataset = "A1_included_sens", pull(res_adr_sens, "TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING")),
  cbind(dataset = "A1_excluded_main", pull(res_adr_main, "GOBP_GLOMERULAR_EPITHELIAL_CELL_DIFFERENTIATION")),
  cbind(dataset = "A1_included_sens", pull(res_adr_sens, "GOBP_GLOMERULAR_EPITHELIAL_CELL_DIFFERENTIATION")),
  cbind(dataset = "A1_excluded_main", pull(res_adr_main, "REACTOME_INTEGRIN_SIGNALING")),
  cbind(dataset = "A1_included_sens", pull(res_adr_sens, "REACTOME_INTEGRIN_SIGNALING"))
)
cat("\n=== GSEA sensitivity (ADR_B_vs_A, Fig.6-equivalent gene sets) ===\n")
print(gsea_sens_compare, row.names = FALSE)

# ---- write outputs ----
write.table(overview, file.path(outdir, "tables/Step5_A1_sensitivity_sig_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(key_tables, file.path(outdir, "tables/Step5_A1_sensitivity_key_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gsea_sens_compare, file.path(outdir, "tables/Step5_A1_sensitivity_GSEA.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

write.table(adr_main, file.path(outdir, "tables/DE_ADR_B_vs_A_A1excluded_main.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(adr_sens, file.path(outdir, "tables/DE_ADR_B_vs_A_A1included_sens.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

wb <- createWorkbook()
addWorksheet(wb, "overview_sig_counts"); writeData(wb, "overview_sig_counts", overview)
addWorksheet(wb, "key_genes_compare"); writeData(wb, "key_genes_compare", key_tables)
addWorksheet(wb, "GSEA_compare"); writeData(wb, "GSEA_compare", gsea_sens_compare)
addWorksheet(wb, "DE_ADR_B_vs_A_A1excluded"); writeDataTable(wb, "DE_ADR_B_vs_A_A1excluded", adr_main[order(adr_main$pvalue), ], tableStyle = "TableStyleLight9")
addWorksheet(wb, "DE_ADR_B_vs_A_A1included"); writeDataTable(wb, "DE_ADR_B_vs_A_A1included", adr_sens[order(adr_sens$pvalue), ], tableStyle = "TableStyleLight9")
saveWorkbook(wb, file.path(outdir, "tables/Supplementary_Step5_A1_sensitivity.xlsx"), overwrite = TRUE)
cat("\nSaved outputs/tables/Supplementary_Step5_A1_sensitivity.xlsx\n")
