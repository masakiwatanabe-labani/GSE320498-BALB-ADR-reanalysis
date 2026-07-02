#!/usr/bin/env Rscript
# DE-tool consistency check: edgeR classic exactTest (the platform's actual
# original method, confirmed from Rscript/08_EnhancedVolcano.R on the shared
# "Amelieff QSP" server -- see Methods_rewrite_and_reviewer_response_notes.md)
# vs. DESeq2 v1.38.3 (this revision's adopted method), on the same deposited
# GSE320498 count matrices.
#
# edgeR call sequence below is copied AS-IS from the discovered template
# script, not from a generic "best practice" edgeR recipe, so this is a
# faithful reproduction of the original method, not a stand-in:
#   keep <- rowSums(cpm(merge) > 2) >= 1          (low-expression filter,
#                                                   computed PER COMPARISON on
#                                                   only that comparison's samples)
#   calcNormFactors(expressed)                     (TMM, edgeR default)
#   estimateCommonDisp(norm) -> estimateTagwiseDisp(com)   (classic 2-step
#                                                   dispersion, NOT the newer
#                                                   unified estimateDisp() --
#                                                   this is what the template
#                                                   actually calls, so it is
#                                                   used here for fidelity)
#   exactTest(mod) -> topTags(out, n=nrow(out$table))
# Significance criterion matches both the manuscript's stated threshold and
# the template's own EnhancedVolcano call: |log2FC| > 1 (fold change > 2) AND
# FDR < 0.05.

suppressMessages({
  library(edgeR)
  library(fgsea)
  library(msigdbr)
})

set.seed(20260220)
indir <- "/usr/local/jupyter/ADR_BALB"
outdir <- file.path(indir, "outputs")

## ---------- load raw counts ----------
counts <- read.delim(file.path(outdir, "tables/00_merged_counts.tsv"), row.names = 1, check.names = FALSE)

## ---------- exact edgeR classic exactTest, as in the discovered template ----------
run_edgeR_exacttest <- function(control_samples, case_samples) {
  control <- counts[, control_samples, drop = FALSE]
  case    <- counts[, case_samples, drop = FALSE]
  merge   <- cbind(control, case)
  group   <- factor(c(rep("A", ncol(control)), rep("B", ncol(case))))  # A=CONTROL, B=CASE (template's own placeholder labels)

  gene <- DGEList(counts = merge, group = group)
  keep <- rowSums(cpm(merge) > 2) >= 1
  expressed <- gene[keep, ]

  norm <- calcNormFactors(expressed)
  if (all(c(length(control_samples), length(case_samples)) == 1)) {
    com <- estimateGLMCommonDisp(norm)
  } else {
    com <- estimateCommonDisp(norm)
  }
  mod <- estimateTagwiseDisp(com)
  out <- exactTest(mod)  # logFC = log2(CASE/CONTROL), matches project sign convention when CASE is chosen correctly below
  final <- topTags(out, n = nrow(out$table))$table
  final$gene <- rownames(final)
  final[, c("gene", "logFC", "logCPM", "PValue", "FDR")]
}

## ---------- comparisons, CONTROL/CASE chosen to match this project's fixed sign convention ----------
comparisons <- list(
  baseline_B_vs_A            = list(control = c("A-Ctrl1","A-Ctrl2","A-Ctrl3"), case = c("B-Ctrl1","B-Ctrl2","B-Ctrl3")),
  ADR_B_vs_A_A1excluded_main = list(control = c("A-ADR2","A-ADR3"),             case = c("B-ADR1","B-ADR2","B-ADR3")),
  ADR_B_vs_A_A1included_sens = list(control = c("A-ADR1","A-ADR2","A-ADR3"),    case = c("B-ADR1","B-ADR2","B-ADR3")),
  A_ADR_vs_Ctrl_A1excluded_main = list(control = c("A-Ctrl1","A-Ctrl2","A-Ctrl3"), case = c("A-ADR2","A-ADR3")),
  A_ADR_vs_Ctrl_A1included_sens = list(control = c("A-Ctrl1","A-Ctrl2","A-Ctrl3"), case = c("A-ADR1","A-ADR2","A-ADR3")),
  B_ADR_vs_Ctrl               = list(control = c("B-Ctrl1","B-Ctrl2","B-Ctrl3"), case = c("B-ADR1","B-ADR2","B-ADR3"))
)

edgeR_results <- list()
for (nm in names(comparisons)) {
  cat("=== edgeR exactTest:", nm, "===\n")
  cfg <- comparisons[[nm]]
  res <- run_edgeR_exacttest(cfg$control, cfg$case)
  edgeR_results[[nm]] <- res
  write.table(res, file.path(outdir, sprintf("tables/edgeR_%s.tsv", nm)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  %d genes tested, %d significant (|log2FC|>1 & FDR<0.05)\n",
              nrow(res), sum(abs(res$logFC) > 1 & res$FDR < 0.05, na.rm = TRUE)))
}
saveRDS(edgeR_results, file.path(outdir, "tables/edgeR_all_results.rds"))

## ---------- Step 2: edgeR vs DESeq2 concordance ----------
deseq2_files <- list(
  baseline_B_vs_A                = "DE_baseline_B_vs_A.tsv",
  ADR_B_vs_A_A1excluded_main     = "DE_ADR_B_vs_A_A1excluded_main.tsv",
  ADR_B_vs_A_A1included_sens     = "DE_ADR_B_vs_A_A1included_sens.tsv",
  A_ADR_vs_Ctrl_A1excluded_main  = "DE_A_ADR_vs_Ctrl.tsv",
  B_ADR_vs_Ctrl                  = "DE_B_ADR_vs_Ctrl.tsv"
)
# A_ADR_vs_Ctrl_A1included_sens has no pre-existing DESeq2 tsv on disk; derive it here for a fair comparison
dds_group_sens <- readRDS(file.path(outdir, "tables/dds_group_sens.rds"))
suppressMessages(library(DESeq2))
r <- results(dds_group_sens, contrast = c("group", "A_ADR", "A_Ctrl"), alpha = 0.05)
d <- as.data.frame(r); d$gene <- rownames(d)
d <- d[, c("gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
colnames(d)[colnames(d) == "log2FoldChange"] <- "log2FC"
write.table(d, file.path(outdir, "tables/DE_A_ADR_vs_Ctrl_A1included_sens.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

deseq2_results <- list()
for (nm in names(deseq2_files)) {
  df <- read.delim(file.path(outdir, "tables", deseq2_files[[nm]]), stringsAsFactors = FALSE)
  # some earlier scripts (06_A1_sensitivity.R) kept DESeq2's native "log2FoldChange"
  # column name instead of renaming to "log2FC" as 02_DE_tables.R does; normalize here
  if ("log2FoldChange" %in% colnames(df) && !"log2FC" %in% colnames(df)) {
    colnames(df)[colnames(df) == "log2FoldChange"] <- "log2FC"
  }
  deseq2_results[[nm]] <- df
}
deseq2_results[["A_ADR_vs_Ctrl_A1included_sens"]] <- d

named_genes <- c("Wt1","Nphs1","Tmem215","Nlrp1b","Glp1r","Serpine1","Loxl1","Col4a1","Col4a2")

concordance_rows <- list()
named_gene_rows <- list()
for (nm in names(edgeR_results)) {
  if (!nm %in% names(deseq2_results)) next
  e <- edgeR_results[[nm]]
  d2 <- deseq2_results[[nm]]
  m <- merge(e, d2, by = "gene", suffixes = c("_edgeR", "_DESeq2"))
  m <- m[!is.na(m$logFC) & !is.na(m$log2FC), ]

  rho_p <- cor(m$logFC, m$log2FC, method = "pearson")
  rho_s <- cor(m$logFC, m$log2FC, method = "spearman")

  sig_edgeR  <- m$gene[abs(m$logFC) > 1 & m$FDR < 0.05 & !is.na(m$FDR)]
  sig_deseq2 <- m$gene[abs(m$log2FC) > 1 & m$padj < 0.05 & !is.na(m$padj)]
  inter <- intersect(sig_edgeR, sig_deseq2)
  uni   <- union(sig_edgeR, sig_deseq2)
  jaccard <- if (length(uni) > 0) length(inter) / length(uni) else NA

  common_sig <- m[m$gene %in% inter, ]
  dir_concord <- if (nrow(common_sig) > 0) mean(sign(common_sig$logFC) == sign(common_sig$log2FC)) else NA

  concordance_rows[[nm]] <- data.frame(
    comparison = nm,
    n_genes_common = nrow(m),
    pearson_r_log2FC = rho_p,
    spearman_rho_log2FC = rho_s,
    n_sig_edgeR = length(sig_edgeR),
    n_sig_DESeq2 = length(sig_deseq2),
    n_sig_common = length(inter),
    jaccard = jaccard,
    direction_concordance_common_sig = dir_concord
  )

  ng <- m[m$gene %in% named_genes, c("gene","logFC","PValue","FDR","log2FC","pvalue","padj")]
  if (nrow(ng) > 0) {
    ng$comparison <- nm
    named_gene_rows[[nm]] <- ng
  }
}
concordance_df <- do.call(rbind, concordance_rows)
rownames(concordance_df) <- NULL
write.table(concordance_df, file.path(outdir, "tables/edgeR_vs_DESeq2_concordance.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n=== Step 2: edgeR vs DESeq2 concordance ===\n")
print(concordance_df, row.names = FALSE, digits = 4)

named_gene_df <- do.call(rbind, named_gene_rows)
rownames(named_gene_df) <- NULL
named_gene_df <- named_gene_df[, c("comparison","gene","logFC","PValue","FDR","log2FC","pvalue","padj")]
write.table(named_gene_df, file.path(outdir, "tables/edgeR_vs_DESeq2_named_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n=== Named genes: edgeR vs DESeq2 ===\n")
print(named_gene_df, row.names = FALSE, digits = 4)

## ---------- Step 3: GSEA impact -- rank by edgeR signed -log10(PValue) vs DESeq2 canonical stat ----------
reac <- msigdbr(db_species = "MM", species = "mouse", collection = "M2", subcollection = "CP:REACTOME")
reactome_list <- split(reac$gene_symbol, reac$gs_name)
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
                list(KARAISKOS2018_PODOCYTE_EXCLUSIVE = unique(karaiskos_exclusive),
                     KARAISKOS2018_PODOCYTE_TOP50 = unique(karaiskos_top50)))
focal_sets <- c("REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION", "REACTOME_INTEGRIN_SIGNALING",
                "KARAISKOS2018_PODOCYTE_EXCLUSIVE", "KARAISKOS2018_PODOCYTE_TOP50")

run_gsea <- function(ranks) {
  set.seed(20260220)
  fgsea(pathways = joint_list, stats = ranks, eps = 0, minSize = 5, maxSize = 500)
}

gsea_compare_rows <- list()
for (nm in c("ADR_B_vs_A_A1excluded_main", "B_ADR_vs_Ctrl")) {
  e <- edgeR_results[[nm]]
  e <- e[!is.na(e$PValue) & !is.na(e$logFC) & e$PValue > 0, ]
  edger_ranks <- sort(setNames(-log10(e$PValue) * sign(e$logFC), e$gene), decreasing = TRUE)
  res_edgeR <- run_gsea(edger_ranks)

  d2 <- deseq2_results[[nm]]
  d2 <- d2[!is.na(d2$stat), ]
  deseq2_ranks <- sort(setNames(d2$stat, d2$gene), decreasing = TRUE)
  res_deseq2 <- run_gsea(deseq2_ranks)

  for (pw in focal_sets) {
    re <- res_edgeR[res_edgeR$pathway == pw, ]
    rd <- res_deseq2[res_deseq2$pathway == pw, ]
    gsea_compare_rows[[paste(nm, pw)]] <- data.frame(
      comparison = nm, pathway = pw,
      NES_edgeR = if (nrow(re)) re$NES else NA, FDR_edgeR = if (nrow(re)) re$padj else NA,
      NES_DESeq2 = if (nrow(rd)) rd$NES else NA, FDR_DESeq2 = if (nrow(rd)) rd$padj else NA,
      sign_concordant = if (nrow(re) && nrow(rd)) sign(re$NES) == sign(rd$NES) else NA
    )
  }
}
gsea_compare_df <- do.call(rbind, gsea_compare_rows)
rownames(gsea_compare_df) <- NULL
write.table(gsea_compare_df, file.path(outdir, "tables/edgeR_vs_DESeq2_GSEA_impact.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n=== Step 3: GSEA impact, edgeR-ranked vs DESeq2-ranked (canonical) ===\n")
print(gsea_compare_df, row.names = FALSE, digits = 4)

cat("\nDone. Wrote edgeR_*.tsv (per-comparison DE), edgeR_vs_DESeq2_concordance.tsv, edgeR_vs_DESeq2_named_genes.tsv, edgeR_vs_DESeq2_GSEA_impact.tsv\n")
