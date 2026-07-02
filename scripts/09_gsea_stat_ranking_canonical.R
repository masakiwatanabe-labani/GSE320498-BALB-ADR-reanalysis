#!/usr/bin/env Rscript
# CANONICAL re-run, superseding 07/08's log2FC-as-primary framing.
#
# Rationale (recorded here per reviewer-response instructions): the manuscript
# Methods text says genes were ranked by "signed log2 fold change", but the
# reported Fig. 6B/6D statistics (nominal P alongside a DIFFERENT FDR q, with
# q < p in one case) are only consistent with a GSEA run on a significance-
# weighted ranking, not a magnitude-only ranking. Step 1 below tests this
# directly: signed(-log10 pvalue) * sign(log2FC) is compared to the DESeq2
# Wald "stat" column. stat IS (up to numerical noise) signed(-log10 p)*sign(FC)
# for a Wald test, so a near-1.0 Spearman correlation confirms that ranking by
# "stat" is the correct reconstruction of whatever the authors actually ran,
# regardless of what the Methods paragraph literally says. This makes stat-
# ranked GSEA the PRIMARY/canonical convention for this project. The log2FC-
# ranked results from 07_gsea_log2FC_ranking_primary.R / 08 are retained as a
# named sensitivity check (matches Methods text verbatim, not the inferred
# actual procedure) -- not deleted, just demoted.
#
# Software: fgsea (established, peer-reviewed package) replaces the project's
# earlier ad hoc preranked implementation. All FDR here is Benjamini-Hochberg
# computed JOINTLY across the full tested gene-set universe per comparison
# (Reactome M2:CP:REACTOME + custom podocyte sets), never focal-subset-only.

suppressMessages({
  library(DESeq2)
  library(fgsea)
  library(msigdbr)
})

set.seed(20260220)
indir <- "/usr/local/jupyter/ADR_BALB"
outdir <- file.path(indir, "outputs")

## ---------- Step 1: ranking-metric identity check ----------
de_files <- c(
  baseline_B_vs_A = "DE_baseline_B_vs_A.tsv",
  ADR_B_vs_A      = "DE_ADR_B_vs_A.tsv",
  A_ADR_vs_Ctrl   = "DE_A_ADR_vs_Ctrl.tsv",
  B_ADR_vs_Ctrl   = "DE_B_ADR_vs_Ctrl.tsv"
)

rank_check_rows <- list()
for (nm in names(de_files)) {
  df <- read.delim(file.path(outdir, "tables", de_files[[nm]]), stringsAsFactors = FALSE)
  df <- df[!is.na(df$stat) & !is.na(df$pvalue) & !is.na(df$log2FC), ]
  signed_sig <- -log10(df$pvalue) * sign(df$log2FC)
  # ties at sign(log2FC)==0 (impossible in practice with continuous FC, but guard anyway)
  rho <- suppressWarnings(cor(signed_sig, df$stat, method = "spearman"))
  rank_check_rows[[nm]] <- data.frame(comparison = nm, n_genes = nrow(df), spearman_rho = rho)
}
rank_check_df <- do.call(rbind, rank_check_rows)
write.table(rank_check_df, file.path(outdir, "tables/Step1_ranking_metric_identity_check.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("=== Step 1: Spearman rho, signed(-log10 p)*sign(log2FC) vs DESeq2 Wald stat ===\n")
print(rank_check_df, row.names = FALSE, digits = 6)
cat("\n")

## ---------- Step 2: gene sets ----------
reac <- msigdbr(db_species = "MM", species = "mouse", collection = "M2", subcollection = "CP:REACTOME")
reactome_list <- split(reac$gene_symbol, reac$gs_name)
db_version <- unique(reac$db_version)

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

focal_sets <- c(
  "TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING",
  "REACTOME_INTEGRIN_SIGNALING",
  "REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS",
  "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
  "KARAISKOS2018_PODOCYTE_EXCLUSIVE",
  "KARAISKOS2018_PODOCYTE_TOP50"
)

version_log <- c(
  sprintf("fgsea package version: %s", as.character(packageVersion("fgsea"))),
  sprintf("msigdbr package version: %s", as.character(packageVersion("msigdbr"))),
  sprintf("DESeq2 package version: %s", as.character(packageVersion("DESeq2"))),
  sprintf("MSigDB db_version (mouse-native, db_species=MM): %s", db_version),
  sprintf("Reactome gene sets (M2:CP:REACTOME): %d", length(reactome_list)),
  sprintf("Custom podocyte sets added: %s", paste(names(custom_sets), collapse = ", ")),
  sprintf("Total joint universe per comparison: %d gene sets", length(joint_list)),
  "Ranking metric (CANONICAL): DESeq2 Wald 'stat' column.",
  "Random seed: 20260220 (reset immediately before every fgsea() call)."
)
writeLines(version_log, file.path(outdir, "logs/09_gsea_canonical_versions.txt"))
cat(paste(version_log, collapse = "\n"), "\n\n")

run_gsea_stat <- function(df) {
  df <- df[!is.na(df$stat), ]
  ranks <- sort(setNames(df$stat, df$gene), decreasing = TRUE)
  set.seed(20260220)
  fgsea(pathways = joint_list, stats = ranks, eps = 0, minSize = 5, maxSize = 500)
}

extract_focal <- function(res, comparison, a1_status) {
  do.call(rbind, lapply(focal_sets, function(pw) {
    r <- res[res$pathway == pw, ]
    data.frame(
      comparison = comparison, a1_status = a1_status, pathway = pw,
      NES = if (nrow(r)) r$NES else NA,
      pval = if (nrow(r)) r$pval else NA,
      FDR_joint = if (nrow(r)) r$padj else NA,
      size = if (nrow(r)) r$size else NA,
      pass_FDR0.05 = if (nrow(r)) (!is.na(r$padj) && r$padj < 0.05) else NA
    )
  }))
}

## ---------- Step 3/4: main (A1-excluded) run, all 4 comparisons ----------
focal_rows <- list()
full_res <- list()
for (nm in names(de_files)) {
  df <- read.delim(file.path(outdir, "tables", de_files[[nm]]), stringsAsFactors = FALSE)
  res <- run_gsea_stat(df)
  full_res[[nm]] <- res
  res_out <- res[order(res$pval), ]
  res_out$leadingEdge <- sapply(res_out$leadingEdge, paste, collapse = ";")
  write.table(res_out, file.path(outdir, sprintf("tables/GSEA_%s_full_joint_canonical.tsv", nm)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  focal_rows[[nm]] <- extract_focal(res, nm, "A1_excluded_main")
}

## ---------- A1 sensitivity: ADR_B_vs_A and A_ADR_vs_Ctrl, A1 included ----------
dds_group_main <- readRDS(file.path(outdir, "tables/dds_group_main.rds"))
dds_group_sens <- readRDS(file.path(outdir, "tables/dds_group_sens.rds"))
get_de <- function(dds, contrast) {
  r <- results(dds, contrast = contrast, alpha = 0.05)
  d <- as.data.frame(r)
  d$gene <- rownames(d)
  d[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
}

sens_contrasts <- list(
  ADR_B_vs_A    = c("group", "B_ADR", "A_ADR"),
  A_ADR_vs_Ctrl = c("group", "A_ADR", "A_Ctrl")
)
for (nm in names(sens_contrasts)) {
  df_sens <- get_de(dds_group_sens, sens_contrasts[[nm]])
  res_sens <- run_gsea_stat(df_sens)
  res_sens_out <- res_sens[order(res_sens$pval), ]
  res_sens_out$leadingEdge <- sapply(res_sens_out$leadingEdge, paste, collapse = ";")
  write.table(res_sens_out, file.path(outdir, sprintf("tables/GSEA_%s_full_joint_canonical_A1included.tsv", nm)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  focal_rows[[paste0(nm, "_A1incl")]] <- extract_focal(res_sens, nm, "A1_included_sens")
}

judgment_table <- do.call(rbind, focal_rows)
rownames(judgment_table) <- NULL
write.table(judgment_table, file.path(outdir, "tables/GSEA_canonical_focal_judgment_table.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("=== Canonical (stat-ranked) focal gene-set judgment table ===\n")
print(judgment_table, row.names = FALSE, digits = 4)
cat("\n")

## ---------- Step 5: Integrin verdict (ADR_B_vs_A, A1-excluded main) ----------
main_adr <- judgment_table[judgment_table$comparison == "ADR_B_vs_A" & judgment_table$a1_status == "A1_excluded_main", ]

cat("=== Step 5: Integrin verdict, ADR_B_vs_A, A1-excluded (main manuscript setting) ===\n")
print(main_adr[main_adr$pathway %in% c("REACTOME_INTEGRIN_SIGNALING",
                                        "REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS",
                                        "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION"), ],
      row.names = FALSE, digits = 4)
cat("\n")

## ---------- Step 6: podocyte set recommendation ----------
cat("=== Step 6: ageing vs Karaiskos, ADR_B_vs_A and B_ADR_vs_Ctrl (A1-excluded main) ===\n")
podo_compare <- judgment_table[judgment_table$comparison %in% c("ADR_B_vs_A", "B_ADR_vs_Ctrl") &
                                  judgment_table$a1_status == "A1_excluded_main" &
                                  judgment_table$pathway %in% c("TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING",
                                                                 "KARAISKOS2018_PODOCYTE_EXCLUSIVE",
                                                                 "KARAISKOS2018_PODOCYTE_TOP50"), ]
print(podo_compare, row.names = FALSE, digits = 4)
write.table(podo_compare, file.path(outdir, "tables/GSEA_canonical_ageing_vs_karaiskos.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# leading-edge gene overlap, ADR_B_vs_A
le <- function(res, pw) unlist(res[res$pathway == pw, ]$leadingEdge)
le_ageing_adr    <- le(full_res[["ADR_B_vs_A"]], "TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING")
le_top50_adr     <- le(full_res[["ADR_B_vs_A"]], "KARAISKOS2018_PODOCYTE_TOP50")
le_ageing_b      <- le(full_res[["B_ADR_vs_Ctrl"]], "TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING")
le_top50_b       <- le(full_res[["B_ADR_vs_Ctrl"]], "KARAISKOS2018_PODOCYTE_TOP50")

cat("\nLeading-edge overlap (ADR_B_vs_A): ageing n=", length(le_ageing_adr),
    " / top50 n=", length(le_top50_adr),
    " / shared=", length(intersect(le_ageing_adr, le_top50_adr)), "\n", sep = "")
cat("Leading-edge overlap (B_ADR_vs_Ctrl): ageing n=", length(le_ageing_b),
    " / top50 n=", length(le_top50_b),
    " / shared=", length(intersect(le_ageing_b, le_top50_b)), "\n", sep = "")

writeLines(c(
  sprintf("ADR_B_vs_A leading-edge (ageing set, n=%d): %s", length(le_ageing_adr), paste(le_ageing_adr, collapse = ",")),
  sprintf("ADR_B_vs_A leading-edge (Karaiskos TOP50, n=%d): %s", length(le_top50_adr), paste(le_top50_adr, collapse = ",")),
  sprintf("ADR_B_vs_A shared leading-edge genes (n=%d): %s", length(intersect(le_ageing_adr, le_top50_adr)), paste(intersect(le_ageing_adr, le_top50_adr), collapse = ",")),
  sprintf("B_ADR_vs_Ctrl leading-edge (ageing set, n=%d): %s", length(le_ageing_b), paste(le_ageing_b, collapse = ",")),
  sprintf("B_ADR_vs_Ctrl leading-edge (Karaiskos TOP50, n=%d): %s", length(le_top50_b), paste(le_top50_b, collapse = ",")),
  sprintf("B_ADR_vs_Ctrl shared leading-edge genes (n=%d): %s", length(intersect(le_ageing_b, le_top50_b)), paste(intersect(le_ageing_b, le_top50_b), collapse = ","))
), file.path(outdir, "logs/09_leading_edge_overlap.txt"))

## ---------- sessionInfo ----------
writeLines(capture.output(sessionInfo()), file.path(outdir, "logs/09_sessionInfo.txt"))
cat("\nWrote canonical (stat-ranked) joint tables, judgment table, ageing-vs-Karaiskos comparison, leading-edge overlap log, sessionInfo.\n")
