#!/usr/bin/env Rscript
# Formal statistics table for Figure 7 (Serpine1, Loxl1, Col4a1, Col4a2),
# baseline and Day-5 ByJcl-vs-AJcl contrasts, extracted verbatim from the
# canonical DESeq2 v1.38.3 run -- no new statistical model. Produces:
#   - tables/TableS_Figure7_stats.csv   (4 genes x 2 contrasts, full DESeq2 columns)
#   - logs/37_Figure7_legend_text.txt   (legend-ready summary text)
#   - logs/37_Figure7_verification.txt (manuscript-text match check +
#     baseline-significance check + CPM/DESeq2 sample-count cross-check)
#
# Day-5 source: DE_ADR_B_vs_A_A1excluded_main.tsv (A-ADR1 excluded, AJcl n=2,
# ByJcl n=3) -- explicit filename for auditability.
# Baseline source: DE_baseline_B_vs_A.tsv (all n=3/3; A-ADR1 is a Day-5
# sample and was never part of this contrast).

suppressMessages(library(dplyr))

indir <- "/usr/local/jupyter/ADR_BALB"
outdir_main <- file.path(indir, "outputs")
outdir_v2 <- file.path(outdir_main, "REVISION_PACKAGE/REVISION_PACKAGE_v2")
tab_dir <- file.path(outdir_v2, "tables")
log_dir <- file.path(outdir_v2, "logs")

genes <- c("Serpine1", "Loxl1", "Col4a1", "Col4a2")
cols <- c("gene", "baseMean", "log2FC", "lfcSE", "stat", "pvalue", "padj")

# ============================================================================
# Step 1: extract the 4 genes from both canonical contrasts
# ============================================================================
baseline <- read.delim(file.path(outdir_main, "tables/DE_baseline_B_vs_A.tsv"), stringsAsFactors = FALSE)
d5 <- read.delim(file.path(outdir_main, "tables/DE_ADR_B_vs_A_A1excluded_main.tsv"), stringsAsFactors = FALSE)
names(d5)[names(d5) == "log2FoldChange"] <- "log2FC"

# cross-check against the generically-named canonical file used elsewhere
# (Figure 6 Panel A, Figure 7 CPM plot) -- must be identical
d5_generic <- read.delim(file.path(outdir_main, "tables/DE_ADR_B_vs_A.tsv"), stringsAsFactors = FALSE)
stopifnot(isTRUE(all.equal(
  d5[d5$gene %in% genes, cols][order(d5$gene[d5$gene %in% genes]), ],
  d5_generic[d5_generic$gene %in% genes, cols][order(d5_generic$gene[d5_generic$gene %in% genes]), ],
  check.attributes = FALSE
)))

base_sub <- baseline[baseline$gene %in% genes, cols]
base_sub$contrast <- "Baseline: ByJcl vs AJcl (n=3 vs 3)"
base_sub$timepoint <- "Baseline"
base_sub$n_AJcl <- 3; base_sub$n_ByJcl <- 3

d5_sub <- d5[d5$gene %in% genes, cols]
d5_sub$contrast <- "Day 5: ByJcl vs AJcl (A-ADR1 excluded; n=2 vs 3)"
d5_sub$timepoint <- "Day 5"
d5_sub$n_AJcl <- 2; d5_sub$n_ByJcl <- 3

out <- rbind(base_sub, d5_sub)
out <- out[order(match(out$gene, genes), match(out$timepoint, c("Baseline", "Day 5"))), ]
out$sign_convention <- "positive log2FC = higher in ByJcl"
out$FDR_method <- "Benjamini-Hochberg, DESeq2::results(), all tested genes in that contrast"
out$significant_FDR0.05 <- ifelse(!is.na(out$padj) & out$padj < 0.05, "YES", "no")

out_full <- out  # keep full precision for the CSV
out_full$baseMean <- round(out_full$baseMean, 2)
out_full$log2FC <- round(out_full$log2FC, 4)
out_full$lfcSE <- round(out_full$lfcSE, 4)
out_full$stat <- round(out_full$stat, 4)
# pvalue/padj kept at full double precision (not rounded) in the CSV, per
# "数値を捏造しない" -- readers can round themselves; sci-notation columns
# added separately for convenience
out_full$pvalue_sci <- formatC(out_full$pvalue, format = "e", digits = 3)
out_full$padj_sci <- formatC(out_full$padj, format = "e", digits = 3)

out_cols <- c("gene", "contrast", "timepoint", "n_AJcl", "n_ByJcl", "baseMean", "log2FC", "lfcSE",
              "stat", "pvalue", "pvalue_sci", "padj", "padj_sci", "significant_FDR0.05",
              "sign_convention", "FDR_method")
out_fp <- file.path(tab_dir, "TableS_Figure7_stats.csv")
write.csv(out_full[, out_cols], out_fp, row.names = FALSE)
cat("Saved:", out_fp, "\n\n")
print(out_full[, c("gene", "timepoint", "n_AJcl", "n_ByJcl", "log2FC", "padj_sci", "significant_FDR0.05")], row.names = FALSE)

# ============================================================================
# Step 2: verify against manuscript text values (Day 5) and check the
# "baseline all n.s." assumption (Step 2 of the task)
# ============================================================================
verif <- c("=== Figure 7 formal statistics -- verification log ===", "")

verif <- c(verif, "--- Day 5 vs manuscript text [ref 85] ---")
text_vals <- data.frame(
  gene = genes,
  log2FC_text = c(1.84, 1.38, 0.97, 1.22),
  padj_text = c(1.7e-9, 4.9e-12, 4.5e-3, 1.5e-3)
)
d5_out <- out_full[out_full$timepoint == "Day 5", ]
all_d5_match <- TRUE
for (g in genes) {
  r <- d5_out[d5_out$gene == g, ]
  t <- text_vals[text_vals$gene == g, ]
  lfc_ok <- abs(r$log2FC - t$log2FC_text) < 0.02
  # padj compared on log10 scale at ~2 significant figures (matches how the
  # manuscript text itself rounds these very small numbers)
  padj_ok <- abs(log10(r$padj) - log10(t$padj_text)) < 0.15
  if (!lfc_ok || !padj_ok) all_d5_match <- FALSE
  verif <- c(verif, sprintf("  %-10s log2FC: canonical=%.4f  text=%.2f  %s  |  padj: canonical=%s  text=%.1e  %s",
                             g, r$log2FC, t$log2FC_text, ifelse(lfc_ok, "MATCH", "*** MISMATCH ***"),
                             r$padj_sci, t$padj_text, ifelse(padj_ok, "MATCH", "*** MISMATCH ***")))
}
verif <- c(verif, "",
  ifelse(all_d5_match,
         "RESULT: all 4 genes' Day-5 log2FC and padj MATCH the manuscript text within rounding tolerance. Continuing.",
         "RESULT: *** at least one gene's Day-5 value MISMATCHES the manuscript text. STOPPING per instructions -- see lines above. ***"))

if (!all_d5_match) {
  writeLines(verif, file.path(log_dir, "37_Figure7_verification.txt"))
  cat(paste(verif, collapse = "\n"), "\n")
  stop("Day-5 values do not match manuscript text -- halted per task instructions ('一致しなければ止めて報告').")
}

verif <- c(verif, "", "--- Baseline significance check (Step 2 assumption: 'all 4 genes padj>0.05 at baseline') ---")
base_out <- out_full[out_full$timepoint == "Baseline", ]
n_sig_baseline <- sum(base_out$padj < 0.05)
for (g in genes) {
  r <- base_out[base_out$gene == g, ]
  verif <- c(verif, sprintf("  %-10s log2FC=%+.4f  padj=%.6f (%s)  -- %s",
                             g, r$log2FC, r$padj, r$padj_sci,
                             ifelse(r$padj < 0.05, "*** SIGNIFICANT (padj<0.05) ***", "n.s. (padj>=0.05)")))
}
verif <- c(verif, "")
if (n_sig_baseline > 0) {
  verif <- c(verif,
    "*** TASK ASSUMPTION MISMATCH: the task states 'baseline は4遺伝子とも padj > 0.05（非有意）であることを確認' --",
    sprintf("this is FALSE for %d of 4 genes. Loxl1 is already significant at baseline", n_sig_baseline),
    "(padj=0.005156 < 0.05), with a smaller effect size than at Day 5 (baseline log2FC=+0.68 vs",
    "Day-5 log2FC=+1.37). This was already flagged in this same project in",
    "29_Figure7_value_check.txt and TableS_ECM_genes_statistics.csv -- not a new finding, but",
    "reported again here explicitly because this task's Step 2 assumed otherwise.",
    "Serpine1 (padj=0.4521), Col4a1 (padj=0.6548), and Col4a2 (padj=0.6385) ARE all non-significant",
    "at baseline as assumed; only Loxl1 is the exception. Actual padj values are reported above and",
    "in TableS_Figure7_stats.csv regardless -- not silently forced to match the 'all n.s.' assumption."
  )
} else {
  verif <- c(verif, "RESULT: all 4 genes are non-significant (padj>=0.05) at baseline, confirming the task's Step 2 assumption.")
}

# ============================================================================
# Step 4: CPM table / DESeq2 sample-count cross-check
# ============================================================================
verif <- c(verif, "", "--- Step 4: TableS_ECM_genes_group_summary.csv vs DESeq2 sample composition ---")
grp <- read.csv(file.path(tab_dir, "TableS_ECM_genes_group_summary.csv"), stringsAsFactors = FALSE)
grp_check <- grp %>% distinct(substrain, timepoint, n) %>% arrange(timepoint, substrain)
expected <- data.frame(
  substrain = c("AJcl", "ByJcl", "AJcl", "ByJcl"),
  timepoint = c("Baseline", "Baseline", "Day 5", "Day 5"),
  n_expected = c(3, 3, 2, 3)
)
merged_check <- merge(grp_check, expected, by = c("substrain", "timepoint"))
merged_check$match <- merged_check$n == merged_check$n_expected
for (i in seq_len(nrow(merged_check))) {
  r <- merged_check[i, ]
  verif <- c(verif, sprintf("  %-6s %-10s CPM-table n=%d  DESeq2-implied n=%d  -- %s",
                             r$substrain, r$timepoint, r$n, r$n_expected,
                             ifelse(r$match, "MATCH", "*** MISMATCH ***")))
}
all_n_match <- all(merged_check$match)
verif <- c(verif, "",
  ifelse(all_n_match,
    "RESULT: CPM summary table sample counts match the DESeq2 A1-excluded run exactly (AJcl Day 5 n=2, all others n=3). Figure and statistics are on the same A1-excluded sample basis.",
    "*** RESULT: sample-count mismatch between CPM table and DESeq2 run -- see above. ***"))

writeLines(verif, file.path(log_dir, "37_Figure7_verification.txt"))
cat("\n", paste(verif, collapse = "\n"), "\n", sep = "")

# ============================================================================
# Step 3: legend-ready text
# ============================================================================
fmt_sci <- function(x) formatC(x, format = "e", digits = 2)
legend_text <- paste0(
  "Baseline: Serpine1 padj=", sprintf("%.3f", base_out$padj[base_out$gene == "Serpine1"]),
  ", Col4a1 padj=", sprintf("%.3f", base_out$padj[base_out$gene == "Col4a1"]),
  ", Col4a2 padj=", sprintf("%.3f", base_out$padj[base_out$gene == "Col4a2"]),
  " (all n.s.); Loxl1 padj=", sprintf("%.4f", base_out$padj[base_out$gene == "Loxl1"]),
  " (already significant at baseline, unlike the other three genes; effect size smaller than at Day 5). ",
  "Day 5 (A-ADR1 excluded): Serpine1 padj=", fmt_sci(d5_out$padj[d5_out$gene == "Serpine1"]),
  " (log2FC=+", sprintf("%.2f", d5_out$log2FC[d5_out$gene == "Serpine1"]), ")",
  ", Loxl1 padj=", fmt_sci(d5_out$padj[d5_out$gene == "Loxl1"]),
  " (log2FC=+", sprintf("%.2f", d5_out$log2FC[d5_out$gene == "Loxl1"]), ")",
  ", Col4a1 padj=", fmt_sci(d5_out$padj[d5_out$gene == "Col4a1"]),
  " (log2FC=+", sprintf("%.2f", d5_out$log2FC[d5_out$gene == "Col4a1"]), ")",
  ", Col4a2 padj=", fmt_sci(d5_out$padj[d5_out$gene == "Col4a2"]),
  " (log2FC=+", sprintf("%.2f", d5_out$log2FC[d5_out$gene == "Col4a2"]), ")",
  " (all significant; DESeq2 Wald test, ByJcl vs AJcl)."
)
writeLines(c(
  "=== Figure 7 legend-ready text (verbatim, transcribable) ===",
  "",
  "NOTE: unlike the task's draft template, this cannot say 'Baseline: ... (all n.s.)' for all",
  "four genes, because Loxl1 is significant at baseline (padj=0.0052). The wording below reflects",
  "the actual canonical values.",
  "",
  legend_text,
  "",
  "Baseline padj range across the 3 non-significant genes: ",
  sprintf("%.3f - %.3f", min(base_out$padj[base_out$gene != "Loxl1"]), max(base_out$padj[base_out$gene != "Loxl1"])),
  "(Serpine1, Col4a1, Col4a2). Loxl1 baseline padj=0.0052 falls outside this range (significant)."
), file.path(log_dir, "37_Figure7_legend_text.txt"))

cat("\n=== Legend text ===\n", legend_text, "\n")
cat("\n=== 37_Figure7_formal_statistics.R complete ===\n")
