#!/usr/bin/env Rscript
# Control-side sensitivity analysis, analogous to the A-ADR1 (A1) sensitivity
# check but for the baseline/Ctrl arm.
#
# Motivation: 18_A1_contamination_QC_figure.R's purity check (podocyte:tubular
# marker-sum CPM ratio) was extended to all 12 samples, not just the 6 ADR
# samples the manuscript's "A1" flag covers. That extension shows A-Ctrl3 has
# the lowest ratio among the 6 Ctrl (baseline) samples (4.50, vs 8.52-50.58 for
# the other 5) and the 2nd-highest summed tubular-marker CPM of ALL 12 samples
# (3483.9, second only to A-ADR1's 3739.5). This was NOT previously flagged by
# the manuscript authors (their "data not shown" sensitivity language only
# covers A1/A-ADR1) -- it is a new QC observation from this re-analysis.
#
# This script builds the equivalent main-vs-sensitivity DE/GSEA comparison for
# baseline_B_vs_A (Fig. 5) and A_ADR_vs_Ctrl (both use the A-Ctrl arm):
#   "main"              : A-Ctrl3 INCLUDED (= the current main analysis, all
#                          6 Ctrl samples; A-ADR1 still excluded, unrelated)
#   "actrl3excluded_sens": A-Ctrl3 EXCLUDED (this sensitivity check)

suppressMessages({
  library(DESeq2)
  library(fgsea)
  library(msigdbr)
  library(openxlsx)
})

set.seed(20260220)
indir <- "/usr/local/jupyter/ADR_BALB"
outdir <- file.path(indir, "outputs")

counts <- read.delim(file.path(outdir, "tables/00_merged_counts.tsv"),
                      row.names = 1, check.names = FALSE)

build_coldata <- function(samples) {
  data.frame(
    sample = samples,
    substrain = factor(ifelse(grepl("^A-", samples), "A", "B"), levels = c("A", "B")),
    treatment = factor(ifelse(grepl("Ctrl", samples), "Ctrl", "ADR"), levels = c("Ctrl", "ADR")),
    row.names = samples
  )
}

make_dds_group <- function(samples) {
  coldata <- build_coldata(samples)
  coldata$group <- factor(paste(coldata$substrain, coldata$treatment, sep = "_"),
                           levels = c("A_Ctrl", "B_Ctrl", "A_ADR", "B_ADR"))
  cts <- as.matrix(counts[, samples])
  dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~group)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  DESeq(dds)
}

all_samples <- colnames(counts)
main_samples <- setdiff(all_samples, "A-ADR1")                    # current main analysis
sens_samples <- setdiff(all_samples, c("A-ADR1", "A-Ctrl3"))      # + A-Ctrl3 excluded

cat("=== Building dds: main (A-Ctrl3 included, n=", length(main_samples),
    ") vs actrl3excluded_sens (A-Ctrl3 excluded, n=", length(sens_samples), ") ===\n", sep = "")

dds_main <- make_dds_group(main_samples)
dds_sens <- make_dds_group(sens_samples)

saveRDS(dds_main, file.path(outdir, "tables/dds_group_main.rds"))   # unchanged, re-saved for parity
saveRDS(dds_sens, file.path(outdir, "tables/dds_group_actrl3excluded_sens.rds"))

get_de <- function(dds, contrast) {
  res <- results(dds, contrast = contrast, alpha = 0.05)
  df <- as.data.frame(res)
  df$gene <- rownames(df)
  df[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
}
colnames_fix <- function(df) { colnames(df)[colnames(df) == "log2FoldChange"] <- "log2FC"; df }

# ---- (1) baseline_B_vs_A: the comparison A-Ctrl3 directly affects (Fig. 5) ----
base_main <- get_de(dds_main, c("group", "B_Ctrl", "A_Ctrl"))
base_sens <- get_de(dds_sens, c("group", "B_Ctrl", "A_Ctrl"))

# ---- (2) A_ADR_vs_Ctrl: within-AJcl ADR response, also uses the A-Ctrl arm ----
a_main <- get_de(dds_main, c("group", "A_ADR", "A_Ctrl"))
a_sens <- get_de(dds_sens, c("group", "A_ADR", "A_Ctrl"))

n_sig <- function(df) sum(df$padj < 0.05, na.rm = TRUE)
overview <- data.frame(
  comparison = c("baseline_B_vs_A", "A_ADR_vs_Ctrl"),
  n_sig_ACtrl3_included_main = c(n_sig(base_main), n_sig(a_main)),
  n_sig_ACtrl3_excluded_sens = c(n_sig(base_sens), n_sig(a_sens))
)
cat("=== A-Ctrl3 sensitivity: significant gene counts (padj<0.05) ===\n")
print(overview, row.names = FALSE)

key_genes <- c("Wt1", "Nphs1", "Tmem215", "Nlrp1b", "Glp1r", "Serpine1", "Col4a1", "Col4a2", "Loxl1")
key_compare <- function(main_df, sens_df, label) {
  m <- merge(main_df[main_df$gene %in% key_genes, c("gene", "log2FoldChange", "padj")],
             sens_df[sens_df$gene %in% key_genes, c("gene", "log2FoldChange", "padj")],
             by = "gene", suffixes = c("_main_ACtrl3incl", "_sens_ACtrl3excl"))
  m$comparison <- label
  m
}
key_tables <- rbind(
  key_compare(base_main, base_sens, "baseline_B_vs_A"),
  key_compare(a_main, a_sens, "A_ADR_vs_Ctrl")
)
key_tables <- key_tables[, c("comparison", "gene", "log2FoldChange_main_ACtrl3incl", "padj_main_ACtrl3incl",
                              "log2FoldChange_sens_ACtrl3excl", "padj_sens_ACtrl3excl")]
cat("\n=== Key gene log2FC/padj, A-Ctrl3-included (main) vs A-Ctrl3-excluded (sens) ===\n")
print(key_tables, row.names = FALSE)

# genome-wide log2FC concordance, baseline_B_vs_A
merged_base <- merge(base_main[, c("gene", "log2FoldChange", "padj")],
                      base_sens[, c("gene", "log2FoldChange", "padj")],
                      by = "gene", suffixes = c("_main", "_sens"))
pearson_r  <- cor(merged_base$log2FoldChange_main, merged_base$log2FoldChange_sens, method = "pearson")
spearman_r <- cor(merged_base$log2FoldChange_main, merged_base$log2FoldChange_sens, method = "spearman")
sig_either <- merged_base[(!is.na(merged_base$padj_main) & merged_base$padj_main < 0.05) |
                           (!is.na(merged_base$padj_sens) & merged_base$padj_sens < 0.05), ]
dir_concord_all <- mean(sign(merged_base$log2FoldChange_main) == sign(merged_base$log2FoldChange_sens), na.rm = TRUE)
dir_concord_sigeither <- mean(sign(sig_either$log2FoldChange_main) == sign(sig_either$log2FoldChange_sens), na.rm = TRUE)

concordance_summary <- c(
  sprintf("baseline_B_vs_A genome-wide log2FC concordance, A-Ctrl3 included (main) vs excluded (sens):"),
  sprintf("  n genes compared: %d", nrow(merged_base)),
  sprintf("  Pearson r = %.4f, Spearman rho = %.4f", pearson_r, spearman_r),
  sprintf("  Direction concordance: %.1f%% genome-wide, %.1f%% among genes significant in either configuration (n=%d)",
          100 * dir_concord_all, 100 * dir_concord_sigeither, nrow(sig_either))
)
cat("\n", paste(concordance_summary, collapse = "\n"), "\n", sep = "")

# ---- GSEA canonical (stat-ranked), baseline_B_vs_A, both configurations ----
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
focal_sets <- c("TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING", "REACTOME_INTEGRIN_SIGNALING",
                "REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS", "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
                "KARAISKOS2018_PODOCYTE_EXCLUSIVE", "KARAISKOS2018_PODOCYTE_TOP50")

run_gsea_stat <- function(df) {
  df <- df[!is.na(df$stat), ]
  ranks <- sort(setNames(df$stat, df$gene), decreasing = TRUE)
  set.seed(20260220)
  fgsea(pathways = joint_list, stats = ranks, eps = 0, minSize = 5, maxSize = 500)
}
extract_focal <- function(res, comparison, actrl3_status) {
  do.call(rbind, lapply(focal_sets, function(pw) {
    r <- res[res$pathway == pw, ]
    data.frame(comparison = comparison, actrl3_status = actrl3_status, pathway = pw,
               NES = if (nrow(r)) r$NES else NA, pval = if (nrow(r)) r$pval else NA,
               FDR_joint = if (nrow(r)) r$padj else NA, size = if (nrow(r)) r$size else NA,
               pass_FDR0.05 = if (nrow(r)) (!is.na(r$padj) && r$padj < 0.05) else NA)
  }))
}

res_base_main <- run_gsea_stat(colnames_fix(base_main))
res_base_sens <- run_gsea_stat(colnames_fix(base_sens))
judgment_table <- rbind(
  extract_focal(res_base_main, "baseline_B_vs_A", "ACtrl3_included_main"),
  extract_focal(res_base_sens, "baseline_B_vs_A", "ACtrl3_excluded_sens")
)
cat("\n=== Canonical (stat-ranked) GSEA, baseline_B_vs_A, A-Ctrl3 included vs excluded ===\n")
print(judgment_table, row.names = FALSE, digits = 4)

# ---- write outputs ----
write.table(colnames_fix(base_main), file.path(outdir, "tables/DE_baseline_B_vs_A_ACtrl3included_main.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(colnames_fix(base_sens), file.path(outdir, "tables/DE_baseline_B_vs_A_ACtrl3excluded_sens.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(overview, file.path(outdir, "tables/Step19_ACtrl3_sensitivity_sig_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(key_tables, file.path(outdir, "tables/Step19_ACtrl3_sensitivity_key_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(judgment_table, file.path(outdir, "tables/GSEA_ACtrl3_sensitivity_judgment_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(concordance_summary, file.path(outdir, "logs/19_ACtrl3_sensitivity_concordance.txt"))

wb <- createWorkbook()
addWorksheet(wb, "overview_sig_counts"); writeData(wb, "overview_sig_counts", overview)
addWorksheet(wb, "key_genes_compare"); writeData(wb, "key_genes_compare", key_tables)
addWorksheet(wb, "GSEA_baseline_compare"); writeData(wb, "GSEA_baseline_compare", judgment_table)
addWorksheet(wb, "DE_baseline_ACtrl3included"); writeDataTable(wb, "DE_baseline_ACtrl3included", colnames_fix(base_main)[order(base_main$pvalue), ], tableStyle = "TableStyleLight9")
addWorksheet(wb, "DE_baseline_ACtrl3excluded"); writeDataTable(wb, "DE_baseline_ACtrl3excluded", colnames_fix(base_sens)[order(base_sens$pvalue), ], tableStyle = "TableStyleLight9")
saveWorkbook(wb, file.path(outdir, "tables/Supplementary_Step19_ACtrl3_sensitivity.xlsx"), overwrite = TRUE)

cat("\n=== 19_ACtrl3_sensitivity_dds_and_DE.R complete ===\n")
