#!/usr/bin/env Rscript
# Follow-up requested by user: test a podocyte gene set derived from a real,
# peer-reviewed, MOUSE-native single-cell dataset (not an annotation-based GO
# term, not an "ageing" signature) -- Karaiskos N et al. 2018, J Am Soc Nephrol
# 29:2060-2068, "A Single-Cell Transcriptome Atlas of the Mouse Glomerulus"
# (GEO GSE111107). Marker genes were extracted directly from the paper's own
# Supplemental Data 1 PDF (Seurat FindAllMarkers output, cluster = podocytes
# vs all other glomerular cells; genes below are exactly as tabulated by the
# authors, sorted by avg_logFC as in their table):
#
#   Table 2 (cell-type EXCLUSIVE markers, used by the authors themselves to
#     identify doublets) -- the strictest, most specific list, n=12
#   Table 3 "Full dataset" (top FindAllMarkers hits for the podocyte cluster,
#     unfiltered beyond the table's own cutoff) -- broader empirical list, n=50
#
# Both are tested here exactly as done for the prior gene sets (05_gsea.R,
# 05b_alt_podocyte_geneset.R): joint fgsea with the 1333 Reactome sets, same
# 4 comparisons, same fixed seed.

suppressMessages({
  library(fgsea)
  library(msigdbr)
})

set.seed(20260220)
indir <- "/usr/local/jupyter/ADR_BALB"
outdir <- file.path(indir, "outputs")

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
cat(sprintf("Karaiskos exclusive markers: n=%d (%d unique)\n",
            length(karaiskos_exclusive), length(unique(karaiskos_exclusive))))
cat(sprintf("Karaiskos top-50 FindAllMarkers list: n=%d (%d unique)\n",
            length(karaiskos_top50), length(unique(karaiskos_top50))))

reac <- msigdbr(db_species = "MM", species = "mouse", collection = "M2", subcollection = "CP:REACTOME")
reactome_list <- split(reac$gene_symbol, reac$gs_name)

joint_list <- c(reactome_list,
                 list(KARAISKOS2018_PODOCYTE_EXCLUSIVE = unique(karaiskos_exclusive),
                      KARAISKOS2018_PODOCYTE_TOP50 = unique(karaiskos_top50)))

comparisons <- c(
  baseline_B_vs_A = "DE_baseline_B_vs_A.tsv",
  ADR_B_vs_A = "DE_ADR_B_vs_A.tsv",
  A_ADR_vs_Ctrl = "DE_A_ADR_vs_Ctrl.tsv",
  B_ADR_vs_Ctrl = "DE_B_ADR_vs_Ctrl.tsv"
)

rows <- list()
gene_in_matrix_report <- NULL
for (nm in names(comparisons)) {
  df <- read.delim(file.path(outdir, "tables", comparisons[[nm]]), stringsAsFactors = FALSE)
  df <- df[!is.na(df$stat), ]
  ranks <- sort(setNames(df$stat, df$gene), decreasing = TRUE)

  if (is.null(gene_in_matrix_report)) {
    missing_excl <- setdiff(unique(karaiskos_exclusive), df$gene)
    missing_top50 <- setdiff(unique(karaiskos_top50), df$gene)
    cat(sprintf("\nGenes from exclusive list NOT found in count matrix: %s\n",
                if (length(missing_excl)) paste(missing_excl, collapse=", ") else "(none)"))
    cat(sprintf("Genes from top-50 list NOT found in count matrix: %s\n",
                if (length(missing_top50)) paste(missing_top50, collapse=", ") else "(none)"))
    gene_in_matrix_report <- TRUE
  }

  set.seed(20260220)
  res <- fgsea(pathways = joint_list, stats = ranks, eps = 0, minSize = 5, maxSize = 500)
  r_excl <- res[res$pathway == "KARAISKOS2018_PODOCYTE_EXCLUSIVE", ]
  r_top50 <- res[res$pathway == "KARAISKOS2018_PODOCYTE_TOP50", ]

  rows[[nm]] <- data.frame(
    comparison = nm,
    n_genesets_tested_joint = nrow(res),
    exclusive_NES = r_excl$NES, exclusive_pval = r_excl$pval, exclusive_FDR = r_excl$padj, exclusive_size = r_excl$size,
    top50_NES = r_top50$NES, top50_pval = r_top50$pval, top50_FDR = r_top50$padj, top50_size = r_top50$size
  )
}

summary_df <- do.call(rbind, rows)
write.table(summary_df, file.path(outdir, "tables/GSEA_karaiskos_podocyte_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n=== Karaiskos et al. 2018 (GSE111107) podocyte marker gene sets: GSEA across all 4 comparisons ===\n")
print(summary_df, row.names = FALSE)

# ---- A1 sensitivity (mirrors 05d_alt_podocyte_A1_sensitivity.R) ----
suppressMessages(library(DESeq2))
dds_group_main <- readRDS(file.path(outdir, "tables/dds_group_main.rds"))
dds_group_sens <- readRDS(file.path(outdir, "tables/dds_group_sens.rds"))

get_de <- function(dds, contrast) {
  res <- results(dds, contrast = contrast, alpha = 0.05)
  df <- as.data.frame(res)
  df$gene <- rownames(df)
  df[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
}
sens_comparisons <- list(
  baseline_B_vs_A = c("group", "B_Ctrl", "A_Ctrl"),
  ADR_B_vs_A      = c("group", "B_ADR", "A_ADR"),
  A_ADR_vs_Ctrl   = c("group", "A_ADR", "A_Ctrl"),
  B_ADR_vs_Ctrl   = c("group", "B_ADR", "B_Ctrl")
)
de_main <- lapply(sens_comparisons, function(ctr) get_de(dds_group_main, ctr))
de_sens <- lapply(sens_comparisons, function(ctr) get_de(dds_group_sens, ctr))

run_gsea_2sets <- function(df) {
  df <- df[!is.na(df$stat), ]
  ranks <- sort(setNames(df$stat, df$gene), decreasing = TRUE)
  set.seed(20260220)
  fgsea(pathways = list(KARAISKOS2018_PODOCYTE_EXCLUSIVE = unique(karaiskos_exclusive),
                         KARAISKOS2018_PODOCYTE_TOP50 = unique(karaiskos_top50)),
        stats = ranks, eps = 0, minSize = 5, maxSize = 500)
}

sens_rows <- list()
for (nm in names(sens_comparisons)) {
  rm_ <- run_gsea_2sets(de_main[[nm]]); rs_ <- run_gsea_2sets(de_sens[[nm]])
  sens_rows[[paste0(nm, "_excl")]] <- data.frame(comparison = nm, pathway = "EXCLUSIVE",
    NES_A1excl = rm_$NES[rm_$pathway=="KARAISKOS2018_PODOCYTE_EXCLUSIVE"],
    pval_A1excl = rm_$pval[rm_$pathway=="KARAISKOS2018_PODOCYTE_EXCLUSIVE"],
    NES_A1incl = rs_$NES[rs_$pathway=="KARAISKOS2018_PODOCYTE_EXCLUSIVE"],
    pval_A1incl = rs_$pval[rs_$pathway=="KARAISKOS2018_PODOCYTE_EXCLUSIVE"])
  sens_rows[[paste0(nm, "_top50")]] <- data.frame(comparison = nm, pathway = "TOP50",
    NES_A1excl = rm_$NES[rm_$pathway=="KARAISKOS2018_PODOCYTE_TOP50"],
    pval_A1excl = rm_$pval[rm_$pathway=="KARAISKOS2018_PODOCYTE_TOP50"],
    NES_A1incl = rs_$NES[rs_$pathway=="KARAISKOS2018_PODOCYTE_TOP50"],
    pval_A1incl = rs_$pval[rs_$pathway=="KARAISKOS2018_PODOCYTE_TOP50"])
}
sens_summary <- do.call(rbind, sens_rows)
rownames(sens_summary) <- NULL
write.table(sens_summary, file.path(outdir, "tables/GSEA_karaiskos_podocyte_A1_sensitivity.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n=== A1 in/out sensitivity ===\n")
print(sens_summary, row.names = FALSE)

# ---- gene-level DE for both marker lists, ADR_B_vs_A (Fig.6A-equivalent) ----
adr_de <- read.delim(file.path(outdir, "tables/DE_ADR_B_vs_A.tsv"), stringsAsFactors = FALSE)
gene_level <- adr_de[adr_de$gene %in% unique(karaiskos_top50), c("gene","log2FC","padj")]
gene_level <- gene_level[order(gene_level$padj), ]
write.table(gene_level, file.path(outdir, "tables/DE_karaiskos_podocyte_genes_ADR_B_vs_A.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n=== Gene-level log2FC/padj for Karaiskos top-50 podocyte genes, ADR_B_vs_A ===\n")
print(gene_level, row.names = FALSE)
