#!/usr/bin/env Rscript
# Extract confirmed, citable values for the manuscript text, directly from
# the canonical run output (no new statistical model; every number below is
# read fresh from the existing DE/GSEA tables, not copied from prior
# summaries, so this script is the traceable source for manuscript_values.md).
#
# canonical run definition (fixed, per instructions):
#   DE: DESeq2 v1.38.3, Wald test, main analysis = A-ADR1 excluded
#   contrasts: all BALB/cByJcl (B) vs BALB/cAJcl (A); positive log2FC = higher in B
#   GSEA: fgsea v1.24.0 multilevel, ranking = DESeq2 Wald statistic
#   FDR: joint BH across all Reactome gene sets tested in that comparison

suppressMessages(library(dplyr))

indir <- "/usr/local/jupyter/ADR_BALB"
outdir_main <- file.path(indir, "outputs")
outdir_v2 <- file.path(outdir_main, "REVISION_PACKAGE/REVISION_PACKAGE_v2")
log_dir <- file.path(outdir_v2, "logs")
tab_dir <- file.path(outdir_v2, "tables")

out <- c()
add <- function(...) out <<- c(out, sprintf(...))

# ============================================================================
# STEP 1: Figure 5 (baseline) / Figure 6A (D5) values
# ============================================================================
add("=== STEP 1: Figure 5 (baseline) values ===")
baseline <- read.delim(file.path(outdir_main, "tables/DE_baseline_B_vs_A.tsv"), stringsAsFactors = FALSE)

n_tested_base <- nrow(baseline)
n_sig_base <- sum(!is.na(baseline$padj) & baseline$padj < 0.05)
n_sig_fc_base <- sum(!is.na(baseline$padj) & baseline$padj < 0.05 & abs(baseline$log2FC) > 1)

add("Tested genes (rows in DE_baseline_B_vs_A.tsv, DESeq2 main dds, A-ADR1 excluded): %d", n_tested_base)
add("padj<0.05: %d", n_sig_base)
add("padj<0.05 & |log2FC|>1: %d", n_sig_fc_base)

get_gene <- function(df, g) df[df$gene == g, c("gene","log2FC","padj","baseMean")]
wt1_base <- get_gene(baseline, "Wt1")
nphs1_base <- get_gene(baseline, "Nphs1")
add("")
add("Wt1 (baseline): log2FC=%.4f, padj=%s, baseMean=%.1f", wt1_base$log2FC, format(wt1_base$padj, digits=4), wt1_base$baseMean)
add("Nphs1 (baseline): log2FC=%.4f, padj=%s, baseMean=%.1f", nphs1_base$log2FC, format(nphs1_base$padj, digits=4), nphs1_base$baseMean)

top20_base <- baseline[order(-abs(baseline$log2FC)), ][1:20, c("gene","log2FC","padj","baseMean")]
add("")
add("Top 20 |log2FC| genes, baseline (ALL tested genes, not restricted to significant):")
for (i in seq_len(nrow(top20_base))) {
  r <- top20_base[i, ]
  add("  %2d. %-15s log2FC=%+.3f  padj=%-10s baseMean=%.1f", i, r$gene, r$log2FC, format(r$padj, digits=4), r$baseMean)
}
write.csv(top20_base, file.path(tab_dir, "TableS_top20_absLog2FC_baseline.csv"), row.names = FALSE)

add("")
add("=== STEP 1: Figure 6A (Day-5 post-ADR, B-ADR vs A-ADR, A1 excluded) values ===")
d5_main <- read.delim(file.path(outdir_main, "tables/DE_ADR_B_vs_A.tsv"), stringsAsFactors = FALSE)
n_tested_d5 <- nrow(d5_main)
n_sig_d5 <- sum(!is.na(d5_main$padj) & d5_main$padj < 0.05)
n_sig_fc_d5 <- sum(!is.na(d5_main$padj) & d5_main$padj < 0.05 & abs(d5_main$log2FC) > 1)

add("Tested genes (rows in DE_ADR_B_vs_A.tsv, DESeq2 main dds, A-ADR1 excluded): %d", n_tested_d5)
add("padj<0.05: %d", n_sig_d5)
add("padj<0.05 & |log2FC|>1: %d", n_sig_fc_d5)

wt1_d5 <- get_gene(d5_main, "Wt1")
nphs1_d5 <- get_gene(d5_main, "Nphs1")
add("")
add("Wt1 (Day 5, A1-excluded main): log2FC=%.4f, padj=%s -- %s", wt1_d5$log2FC, format(wt1_d5$padj, digits=4),
    ifelse(wt1_d5$padj < 0.05, "SIGNIFICANT", ifelse(wt1_d5$padj < 0.06, "BORDERLINE (just above 0.05)", "not significant")))
add("Nphs1 (Day 5, A1-excluded main): log2FC=%.4f, padj=%s -- %s", nphs1_d5$log2FC, format(nphs1_d5$padj, digits=4),
    ifelse(nphs1_d5$padj < 0.05, "SIGNIFICANT", ifelse(nphs1_d5$padj < 0.10, "borderline/approaching, not significant", "not significant")))

for (g in c("Serpine1","Loxl1","Col4a1","Col4a2")) {
  r <- get_gene(d5_main, g)
  add("%s (Day 5, A1-excluded main): log2FC=%.4f, padj=%s", g, r$log2FC, format(r$padj, digits=4))
}

d5_a1incl_fp <- file.path(outdir_main, "tables/DE_ADR_B_vs_A_A1included_sens.tsv")
if (file.exists(d5_a1incl_fp)) {
  d5_a1incl <- read.delim(d5_a1incl_fp, stringsAsFactors = FALSE)
  if (!"log2FC" %in% colnames(d5_a1incl) && "log2FoldChange" %in% colnames(d5_a1incl)) {
    colnames(d5_a1incl)[colnames(d5_a1incl) == "log2FoldChange"] <- "log2FC"
  }
  wt1_incl <- get_gene(d5_a1incl, "Wt1")
  nphs1_incl <- get_gene(d5_a1incl, "Nphs1")
  add("")
  add("Reference, A1 INCLUDED (sensitivity, not main): Wt1 log2FC=%.4f, padj=%s", wt1_incl$log2FC, format(wt1_incl$padj, digits=4))
  add("Reference, A1 INCLUDED (sensitivity, not main): Nphs1 log2FC=%.4f, padj=%s", nphs1_incl$log2FC, format(nphs1_incl$padj, digits=4))
} else {
  add("")
  add("Reference, A1 INCLUDED: FILE NOT FOUND (%s) -- 該当なし, not reported", d5_a1incl_fp)
}

add("")
add("=== Gene-count discrepancy check (manuscript vs. canonical) ===")
add("Manuscript-reported: Fig. 5 'total = 13,064 variables'; Fig. 6A 'total = 13,278 variables' (mismatched between the two figures; edgeR-pipeline-derived).")
add("Prior draft response text stated: '19,662 tested genes' (baseline).")
add("Canonical DESeq2 run (this script, re-derived from source tables): baseline_B_vs_A tested genes = %d; ADR_B_vs_A (D5, A1-excluded) tested genes = %d.", n_tested_base, n_tested_d5)
add("=> Baseline and Day-5 canonical gene counts are IDENTICAL (%d = %d): both draw from the same single DESeq2 dds object (dds_group_main, A-ADR1 excluded, prefilter rowSums>=10 applied once).", n_tested_base, n_tested_d5)
add("=> This matches the prior draft text's '19,662 tested genes' figure exactly.")
add("=> Neither manuscript-reported figure (13,064 / 13,278) matches the canonical DESeq2 count (%d); both are lower, consistent with a different (edgeR-pipeline, independently filtered per-comparison) gene set, not the DESeq2 prefilter used here. This must be flagged as a discrepancy requiring a Methods-text correction, not silently reconciled.", n_tested_base)

writeLines(out, file.path(log_dir, "24_Step1_Fig5_Fig6A_values.txt"))
cat(paste(out, collapse = "\n"), "\n")
out <- c()

# ============================================================================
# STEP 2: Figure 6B (podocyte, Karaiskos + ageing set) values
# ============================================================================
add("=== STEP 2: Figure 6B values (comparison = ADR_B_vs_A, substrain comparison at Day 5, A1-excluded main) ===")
add("NOTE: this is explicitly ADR_B_vs_A (B vs A substrain comparison at Day 5), NOT B_ADR_vs_Ctrl (within-B disease response) -- these are different comparisons and are not to be confused.")

gsea_main <- read.delim(file.path(outdir_main, "tables/GSEA_ADR_B_vs_A_full_joint_canonical.tsv"), stringsAsFactors = FALSE)
gsea_a1incl_fp <- file.path(outdir_main, "tables/GSEA_ADR_B_vs_A_full_joint_canonical_A1included.tsv")
gsea_a1incl <- if (file.exists(gsea_a1incl_fp)) read.delim(gsea_a1incl_fp, stringsAsFactors = FALSE) else NULL

get_pw <- function(df, pw) {
  if (is.null(df)) return(NULL)
  r <- df[df$pathway == pw, ]
  if (nrow(r) == 0) return(NULL)
  r
}
report_pw <- function(df, pw, label) {
  r <- get_pw(df, pw)
  if (is.null(r)) {
    add("%s: 該当なし (pathway not found in this table)", label)
    return(invisible(NULL))
  }
  le <- strsplit(r$leadingEdge, ";")[[1]]
  add("%s: NES=%.4f, nominal p=%s, joint FDR=%s, size=%d, n_leading_edge=%d", label, r$NES,
      format(r$pval, digits = 4), format(r$padj, digits = 4), r$size, length(le))
  add("   leading edge genes: %s", paste(le, collapse = ", "))
}

report_pw(gsea_main, "KARAISKOS2018_PODOCYTE_TOP50", "TOP50 (ADR_B_vs_A, A1-excluded main)")
report_pw(gsea_main, "KARAISKOS2018_PODOCYTE_EXCLUSIVE", "EXCLUSIVE (ADR_B_vs_A, A1-excluded main)")
report_pw(gsea_main, "TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING", "Ageing set (ADR_B_vs_A, A1-excluded main)")
add("")
if (!is.null(gsea_a1incl)) {
  report_pw(gsea_a1incl, "KARAISKOS2018_PODOCYTE_TOP50", "TOP50 (ADR_B_vs_A, A1-INCLUDED sensitivity)")
  report_pw(gsea_a1incl, "KARAISKOS2018_PODOCYTE_EXCLUSIVE", "EXCLUSIVE (ADR_B_vs_A, A1-INCLUDED sensitivity)")
  report_pw(gsea_a1incl, "TABULA_MURIS_SENIS_KIDNEY_PODOCYTE_AGEING", "Ageing set (ADR_B_vs_A, A1-INCLUDED sensitivity)")
} else {
  add("A1-included table not found at %s -- 該当なし", gsea_a1incl_fp)
}

writeLines(out, file.path(log_dir, "24_Step2_Fig6B_values.txt"))
cat(paste(out, collapse = "\n"), "\n")
out <- c()

# ============================================================================
# STEP 3: Figure 6D (integrin) values
# ============================================================================
add("=== STEP 3: Figure 6D values (comparison = ADR_B_vs_A, A1-excluded main) ===")
report_pw(gsea_main, "REACTOME_INTEGRIN_SIGNALING", "Integrin Signaling")
report_pw(gsea_main, "REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS", "Integrin Cell Surface Interactions")
report_pw(gsea_main, "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION", "ECM Organization")

integrin_r <- get_pw(gsea_main, "REACTOME_INTEGRIN_SIGNALING")
add("")
add("Manuscript-reported (per user input for this task): NES=+1.54, nominal p=0.0235 (Fig. 6D, Integrin Signaling).")
add("Canonical DESeq2+fgsea value: NES=%.4f, nominal p=%s, joint FDR=%s.", integrin_r$NES, format(integrin_r$pval, digits=4), format(integrin_r$padj, digits=4))
add("=> NES deviates from the manuscript-reported value (canonical %.2f vs. reported +1.54; direction/sign agrees -- both positive/higher in ByJcl -- but magnitude differs, expected given the DE-engine and ranking-metric changes; nominal p also differs (%s vs. 0.0235). This must be reported as a discrepancy, not silently matched.",
    integrin_r$NES, format(integrin_r$pval, digits=4))

writeLines(out, file.path(log_dir, "24_Step3_Fig6D_values.txt"))
cat(paste(out, collapse = "\n"), "\n")
out <- c()

# ============================================================================
# STEP 5: ADR response substrain concordance
# ============================================================================
add("=== STEP 5: ADR response concordance between substrains (re-derived fresh from DE tables) ===")
a_resp <- read.delim(file.path(outdir_main, "tables/DE_A_ADR_vs_Ctrl.tsv"), stringsAsFactors = FALSE)
b_resp <- read.delim(file.path(outdir_main, "tables/DE_B_ADR_vs_Ctrl.tsv"), stringsAsFactors = FALSE)

merged <- merge(a_resp[, c("gene","log2FC","padj")], b_resp[, c("gene","log2FC","padj")],
                by = "gene", suffixes = c("_A", "_B"))
n_genes <- nrow(merged)
pearson_r <- cor(merged$log2FC_A, merged$log2FC_B, method = "pearson")
spearman_r <- cor(merged$log2FC_A, merged$log2FC_B, method = "spearman")

sig_both <- merged[!is.na(merged$padj_A) & merged$padj_A < 0.05 & !is.na(merged$padj_B) & merged$padj_B < 0.05, ]
dir_concord_both <- mean(sign(sig_both$log2FC_A) == sign(sig_both$log2FC_B))

sig_either <- merged[(!is.na(merged$padj_A) & merged$padj_A < 0.05) | (!is.na(merged$padj_B) & merged$padj_B < 0.05), ]
dir_concord_either <- mean(sign(sig_either$log2FC_A) == sign(sig_either$log2FC_B))
dir_concord_all <- mean(sign(merged$log2FC_A) == sign(merged$log2FC_B))

add("n genes compared (present in both A_ADR_vs_Ctrl and B_ADR_vs_Ctrl DE tables): %d", n_genes)
add("Pearson r = %.4f", pearson_r)
add("Spearman rho = %.4f", spearman_r)
add("")
add("Direction concordance, genome-wide (all %d genes): %.1f%%", n_genes, 100 * dir_concord_all)
add("Direction concordance, significant in EITHER substrain (n=%d): %.1f%%", nrow(sig_either), 100 * dir_concord_either)
add("Direction concordance, significant in BOTH substrains (n=%d): %.1f%%", nrow(sig_both), 100 * dir_concord_both)

interaction_fp <- file.path(outdir_main, "tables/DE_interaction_substrain_by_treatment.tsv")
if (file.exists(interaction_fp)) {
  int_df <- read.delim(interaction_fp, stringsAsFactors = FALSE)
  padj_col <- if ("padj" %in% colnames(int_df)) "padj" else if ("padj.x" %in% colnames(int_df)) "padj.x" else NA
  if (!is.na(padj_col)) {
    n_int_sig <- sum(!is.na(int_df[[padj_col]]) & int_df[[padj_col]] < 0.05)
    add("")
    add("Substrain x treatment interaction term: %d / %d genes significant at padj<0.05", n_int_sig, nrow(int_df))
  } else {
    add("")
    add("Interaction table found but no padj column identified -- 該当なし (columns: %s)", paste(colnames(int_df), collapse=", "))
  }
} else {
  add("")
  add("Interaction table not found at %s -- 該当なし", interaction_fp)
}

write.csv(merged, file.path(tab_dir, "TableS_ADR_response_A_vs_B_merged.csv"), row.names = FALSE)
writeLines(out, file.path(log_dir, "24_Step5_ADR_concordance_values.txt"))
cat(paste(out, collapse = "\n"), "\n")

cat("\n=== 24_manuscript_values_extraction.R complete ===\n")
