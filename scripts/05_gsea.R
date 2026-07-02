#!/usr/bin/env Rscript
# Step 4 (R1-6, R2-2, R2-4): preranked GSEA for all four DE comparisons.
#
# Ranking metric: DESeq2 Wald "stat" column (signed test statistic; captures
# both direction and evidence strength without additional shrinkage
# assumptions). This choice is recorded here and in the manuscript-facing log.
#
# Gene-set databases (versions recorded for methods/citation):
#   - msigdbr package version: recorded below (programmatically)
#   - MSigDB db_version (mouse-native, db_species="MM"): recorded below
#   - Reactome pathways: MSigDB M2:CP:REACTOME collection (mouse-native symbols)
#   - Podocyte signature: MSigDB M8 collection, gene set
#     "TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING" (same GMT identity as used
#     in the original manuscript; NOTE current MSigDB release may have a
#     slightly different gene count than the version used in the original
#     manuscript's analysis -- reported explicitly, not silently reconciled)
#   - fgsea package version: recorded below
#
# Two ways of testing the podocyte gene set are BOTH reported because they
# answer different questions and give different FDR/q values:
#   (1) "joint" -- podocyte set tested together with all 1333 Reactome sets
#       in a single BH correction (1334 tests). This is the statistically
#       appropriate way to report a family-wise FDR/q value and is the
#       PRIMARY result used to answer "does this pass FDR<0.05".
#   (2) "single-set" -- podocyte set tested alone (as apparently done in the
#       original manuscript, which reported a "single set estimate" q-value).
#       Reported for direct comparability with the original manuscript text,
#       but this is not a proper multiple-testing correction across the
#       tested universe of gene sets and should not be the basis for claims
#       of genome-wide-corrected significance.

suppressMessages({
  library(fgsea)
  library(msigdbr)
  library(dplyr)
})

set.seed(20260220)

indir <- "/usr/local/jupyter/ADR_BALB"
outdir <- file.path(indir, "outputs")

# ---- gene sets ----
reac <- msigdbr(db_species = "MM", species = "mouse", collection = "M2", subcollection = "CP:REACTOME")
reactome_list <- split(reac$gene_symbol, reac$gs_name)

m8 <- msigdbr(db_species = "MM", species = "mouse", collection = "M8")
podo_genes <- unique(m8$gene_symbol[m8$gs_name == "TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING"])
podo_list <- list(TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING = podo_genes)

joint_list <- c(reactome_list, podo_list)

db_version <- unique(reac$db_version)

version_log <- c(
  sprintf("fgsea package version: %s", as.character(packageVersion("fgsea"))),
  sprintf("msigdbr package version: %s", as.character(packageVersion("msigdbr"))),
  sprintf("MSigDB db_version (mouse-native, db_species=MM): %s", db_version),
  sprintf("Reactome gene sets tested (M2:CP:REACTOME): %d", length(reactome_list)),
  sprintf("Podocyte gene set (M8 TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING) size in current MSigDB release: %d genes", length(podo_genes)),
  "Original manuscript reported gene set size = 192 for this same named gene set (Fig. 6B legend); current release differs (see note above) -- MSigDB periodically revises gene set membership between releases.",
  sprintf("Total gene sets in joint preranked test: %d", length(joint_list)),
  "Ranking metric: DESeq2 Wald 'stat' (signed test statistic) from each comparison's main-dataset DE table.",
  "Random seed: 20260220 (fixed for fgsea's permutation/multilevel estimation, for reproducibility)."
)
writeLines(version_log, file.path(outdir, "logs/05_gsea_versions.txt"))
cat(paste(version_log, collapse = "\n"), "\n\n")

run_gsea_for_comparison <- function(de_file, comparison_name) {
  df <- read.delim(file.path(outdir, "tables", de_file), stringsAsFactors = FALSE)
  df <- df[!is.na(df$stat), ]
  # collapse potential duplicate gene symbols by keeping the max |stat| entry (none expected; Geneid was unique in Step 0)
  ranks <- setNames(df$stat, df$gene)
  ranks <- sort(ranks, decreasing = TRUE)

  set.seed(20260220)
  res_joint <- fgsea(pathways = joint_list, stats = ranks, eps = 0, minSize = 5, maxSize = 500)
  res_joint <- res_joint[order(res_joint$pval), ]
  res_joint$leadingEdge <- sapply(res_joint$leadingEdge, paste, collapse = ";")

  set.seed(20260220)
  res_podo_single <- fgsea(pathways = podo_list, stats = ranks, eps = 0, minSize = 5, maxSize = 500)
  res_podo_single$leadingEdge <- sapply(res_podo_single$leadingEdge, paste, collapse = ";")

  write.table(res_joint, file.path(outdir, paste0("tables/GSEA_", comparison_name, "_full_joint.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  podo_joint_row <- res_joint[res_joint$pathway == "TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING", ]
  integrin_joint_row <- res_joint[res_joint$pathway == "REACTOME_INTEGRIN_SIGNALING", ]

  list(
    comparison = comparison_name,
    n_genesets_tested = nrow(res_joint),
    podo_joint = podo_joint_row,
    podo_single = res_podo_single,
    integrin_joint = integrin_joint_row,
    res_joint = res_joint
  )
}

comparisons <- c(
  baseline_B_vs_A = "DE_baseline_B_vs_A.tsv",
  ADR_B_vs_A = "DE_ADR_B_vs_A.tsv",
  A_ADR_vs_Ctrl = "DE_A_ADR_vs_Ctrl.tsv",
  B_ADR_vs_Ctrl = "DE_B_ADR_vs_Ctrl.tsv"
)

all_results <- list()
summary_rows <- list()

for (nm in names(comparisons)) {
  cat("=== Running GSEA for:", nm, "===\n")
  r <- run_gsea_for_comparison(comparisons[[nm]], nm)
  all_results[[nm]] <- r

  pj <- r$podo_joint
  ps <- r$podo_single
  ij <- r$integrin_joint

  summary_rows[[nm]] <- data.frame(
    comparison = nm,
    n_genesets_tested_joint = r$n_genesets_tested,
    podocyte_NES_joint = if (nrow(pj) > 0) pj$NES else NA,
    podocyte_pval_joint = if (nrow(pj) > 0) pj$pval else NA,
    podocyte_FDR_joint = if (nrow(pj) > 0) pj$padj else NA,
    podocyte_NES_single = if (nrow(ps) > 0) ps$NES else NA,
    podocyte_pval_single = if (nrow(ps) > 0) ps$pval else NA,
    integrin_NES_joint = if (nrow(ij) > 0) ij$NES else NA,
    integrin_pval_joint = if (nrow(ij) > 0) ij$pval else NA,
    integrin_FDR_joint = if (nrow(ij) > 0) ij$padj else NA
  )
}

summary_df <- do.call(rbind, summary_rows)
write.table(summary_df, file.path(outdir, "tables/GSEA_podocyte_integrin_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n=== Summary: podocyte + integrin gene sets across all comparisons ===\n")
print(summary_df, row.names = FALSE)

saveRDS(all_results, file.path(outdir, "tables/gsea_all_results.rds"))
cat("\nSaved full GSEA result tables (outputs/tables/GSEA_<comparison>_full_joint.tsv) and summary (GSEA_podocyte_integrin_summary.tsv)\n")
