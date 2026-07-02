#!/usr/bin/env Rscript
# Assemble Step 4 GSEA supplementary workbook: full results per comparison +
# summary sheets, with captions documenting methods/versions and the
# podocyte gene-set sign discrepancy investigation.

suppressMessages(library(openxlsx))

indir <- "/usr/local/jupyter/ADR_BALB"
outdir <- file.path(indir, "outputs")

version_log <- readLines(file.path(outdir, "logs/05_gsea_versions.txt"))
podo_integrin_summary <- read.delim(file.path(outdir, "tables/GSEA_podocyte_integrin_summary.tsv"), stringsAsFactors = FALSE)
alt_summary <- read.delim(file.path(outdir, "tables/GSEA_alt_podocyte_identity_vs_ageing.tsv"), stringsAsFactors = FALSE)

wb <- createWorkbook()

addWorksheet(wb, "README_methods")
readme <- c(
  "Preranked GSEA methods and versions:",
  version_log,
  "",
  "IMPORTANT FINDING (podocyte gene-set sign check, Day-5 ADR ByJcl-vs-AJcl / Fig.6B-equivalent contrast):",
  "The MSigDB M8 set TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING (same gene set used in the original manuscript) reproduces with the OPPOSITE sign from the manuscript's reported NES=-1.42:",
  "  This reanalysis: NES=+1.93, nominal p=8.1e-7, joint FDR=9.5e-6 (statistically significant, opposite direction).",
  "  Verified not a coding artifact: DESeq2 stat sign always matches log2FC sign; enrichment plot (see figures/CHECK_enrichment_podocyte_ADR_B_vs_A.png) shows a clean, genuine positive enrichment.",
  "  Likely explanation: this MSigDB set is derived from an AGEING comparison within podocytes (Tabula Muris Senis), not a podocyte-identity signature; its leading-edge genes here are MHC/immune/senescence genes (H2-K1, H2-D1, B2m, Bst2, Ly6a), not canonical podocyte markers (Wt1/Nphs1/Nphs2/Podxl are not leading-edge).",
  "  As an EXPLORATORY follow-up (not part of the original manuscript methodology), we additionally tested a pre-existing curated GO Biological Process gene set representing canonical podocyte identity: GOBP_GLOMERULAR_EPITHELIAL_CELL_DIFFERENTIATION (MSigDB M5:GO:BP, n=24 genes: Nphs1, Nphs2, Podxl, Wt1, Cd2ap, Magi2, Ptpro, Foxc2, etc.).",
  "  This identity-focused gene set DOES show the direction consistent with the manuscript's narrative at Day-5 post-ADR (NES=-1.86, FDR=0.009) while showing NO significant difference at baseline (NES=+1.17, FDR=0.44) -- see sheet 'alt_podocyte_identity_vs_ageing'.",
  "  Recommendation for the authors: (1) re-check the original GSEA pipeline/sign convention used for Fig.6B against the deposited count matrices, since the originally reported sign does not reproduce with the same named gene set from current MSigDB; (2) consider reporting/cross-validating with a canonical podocyte-identity gene set (e.g. GOBP_GLOMERULAR_EPITHELIAL_CELL_DIFFERENTIATION) alongside or instead of the ageing-derived set, since it gives a cleaner and more directly interpretable result.",
  "",
  "All full per-comparison GSEA result tables (every tested gene set, NES/pval/FDR/leadingEdge) are in the 'GSEA_<comparison>' sheets."
)
writeData(wb, "README_methods", data.frame(README = readme))
setColWidths(wb, "README_methods", cols = 1, widths = 160)

addWorksheet(wb, "podocyte_integrin_summary")
writeData(wb, "podocyte_integrin_summary", podo_integrin_summary)

addWorksheet(wb, "alt_podocyte_identity_vs_ageing")
writeData(wb, "alt_podocyte_identity_vs_ageing", alt_summary)

for (nm in c("baseline_B_vs_A", "ADR_B_vs_A", "A_ADR_vs_Ctrl", "B_ADR_vs_Ctrl")) {
  df <- read.delim(file.path(outdir, paste0("tables/GSEA_", nm, "_full_joint.tsv")), stringsAsFactors = FALSE)
  addWorksheet(wb, paste0("GSEA_", nm))
  writeDataTable(wb, paste0("GSEA_", nm), df, tableStyle = "TableStyleLight9")
}

saveWorkbook(wb, file.path(outdir, "tables/Supplementary_GSEA_tables.xlsx"), overwrite = TRUE)
cat("Saved outputs/tables/Supplementary_GSEA_tables.xlsx\n")
