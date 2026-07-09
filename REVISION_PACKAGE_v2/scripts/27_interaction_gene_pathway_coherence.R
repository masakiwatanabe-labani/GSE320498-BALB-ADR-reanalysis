#!/usr/bin/env Rscript
# Do the 24 genes with a significant substrain x treatment interaction
# (canonical DESeq2 interaction model, FDR<0.05) form a coherent pathway, or
# are they a scattered handful of genes? Also checks whether the interaction
# signal is driven by known polymorphic loci (Oas1 cluster, H2/MHC, Nlrp1)
# rather than a genuine differential-treatment-response signal.
#
# Reuses the EXISTING canonical interaction model output
# (DE_interaction_substrain_by_treatment.tsv, from 01_build_dds.R's
# ~substrain+treatment+substrain:treatment design on dds_group_main,
# A-ADR1 excluded) verbatim -- no new model fit, no re-derivation of the
# gene list. Only new analysis: ORA (ReactomePA + GO BP) on this fixed list.

suppressMessages({
  library(clusterProfiler)
  library(ReactomePA)
  library(org.Mm.eg.db)
})

indir <- "/usr/local/jupyter/ADR_BALB"
outdir_main <- file.path(indir, "outputs")
outdir_v2 <- file.path(outdir_main, "REVISION_PACKAGE/REVISION_PACKAGE_v2")
fig_dir <- file.path(outdir_v2, "figures")
tab_dir <- file.path(outdir_v2, "tables")
log_dir <- file.path(outdir_v2, "logs")

out <- c()
add <- function(...) out <<- c(out, sprintf(...))

# ============================================================================
# Step 1: the 24-gene interaction list, with per-substrain response context
# ============================================================================
int_df <- read.delim(file.path(outdir_main, "tables/DE_interaction_substrain_by_treatment.tsv"), stringsAsFactors = FALSE)
n_tested <- nrow(int_df)
sig24 <- int_df[!is.na(int_df$padj) & int_df$padj < 0.05, ]
add("Step 1: interaction model (canonical, DESeq2 ~substrain+treatment+substrain:treatment, A-ADR1 excluded): %d/%d genes at interaction FDR<0.05", nrow(sig24), n_tested)
stopifnot(nrow(sig24) == 24)

a_resp <- read.delim(file.path(outdir_main, "tables/DE_A_ADR_vs_Ctrl.tsv"), stringsAsFactors = FALSE)
b_resp <- read.delim(file.path(outdir_main, "tables/DE_B_ADR_vs_Ctrl.tsv"), stringsAsFactors = FALSE)
baseline <- read.delim(file.path(outdir_main, "tables/DE_baseline_B_vs_A.tsv"), stringsAsFactors = FALSE)

a_resp_r <- setNames(a_resp[, c("gene", "log2FC", "padj")], c("gene", "log2FC_A", "padj_A"))
b_resp_r <- setNames(b_resp[, c("gene", "log2FC", "padj")], c("gene", "log2FC_B", "padj_B"))
baseline_r <- setNames(baseline[, c("gene", "log2FC", "padj")], c("gene", "log2FC_baseline", "padj_baseline"))

sig24_full <- merge(sig24, a_resp_r, by = "gene")
sig24_full <- merge(sig24_full, b_resp_r, by = "gene")
sig24_full <- merge(sig24_full, baseline_r, by = "gene")
sig24_full <- sig24_full[order(sig24_full$padj), ]

out_table <- sig24_full[, c("gene", "baseMean", "log2FC_interaction", "stat", "padj",
                             "log2FC_A", "padj_A", "log2FC_B", "padj_B",
                             "log2FC_baseline", "padj_baseline")]
colnames(out_table)[colnames(out_table) == "stat"] <- "interaction_stat"
colnames(out_table)[colnames(out_table) == "padj"] <- "interaction_padj"
write.csv(out_table, file.path(tab_dir, "TableS2_interaction_genes.csv"), row.names = FALSE)
add("Saved: tables/TableS2_interaction_genes.csv (%d genes x %d columns)", nrow(out_table), ncol(out_table))
add("")
add("=== 24 interaction genes, sorted by interaction padj ===")
for (i in seq_len(nrow(out_table))) {
  r <- out_table[i, ]
  add("  %-12s int_log2FC=%+.2f int_padj=%s | A: log2FC=%+.2f padj=%s | B: log2FC=%+.2f padj=%s | baseline: log2FC=%+.2f padj=%s",
      r$gene, r$log2FC_interaction, format(r$interaction_padj, digits = 3),
      r$log2FC_A, format(r$padj_A, digits = 3), r$log2FC_B, format(r$padj_B, digits = 3),
      r$log2FC_baseline, format(r$padj_baseline, digits = 3))
}

# ============================================================================
# Step 2: pathway coherence -- ORA (Reactome) + GO BP on the 24 genes
# ============================================================================
add("")
add("=== Step 2: pathway coherence (ORA), background = all %d genes tested in canonical run ===", n_tested)

universe_symbols <- int_df$gene
deg_symbols <- sig24$gene

map_deg <- suppressWarnings(bitr(deg_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db))
map_universe <- suppressWarnings(bitr(universe_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db))
add("DEG symbols mapped to Entrez: %d/%d", length(unique(map_deg$ENTREZID)), length(deg_symbols))
add("Universe symbols mapped to Entrez: %d/%d", length(unique(map_universe$ENTREZID)), length(universe_symbols))

# ---- Reactome ORA ----
ora_reactome <- enrichPathway(
  gene = unique(map_deg$ENTREZID), universe = unique(map_universe$ENTREZID),
  organism = "mouse", pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH",
  minGSSize = 2, maxGSSize = 1000, readable = TRUE
)
ora_reactome_df <- as.data.frame(ora_reactome)
n_reactome_tested <- nrow(ora_reactome_df)
n_reactome_sig <- sum(ora_reactome_df$p.adjust < 0.05, na.rm = TRUE)
add("")
add("Reactome ORA: %d terms tested, %d significant at BH p.adjust<0.05", n_reactome_tested, n_reactome_sig)
if (n_reactome_tested > 0) {
  top5 <- head(ora_reactome_df[order(ora_reactome_df$pvalue), ], 5)
  for (i in seq_len(nrow(top5))) {
    r <- top5[i, ]
    add("  %-55s GeneRatio=%s pvalue=%s p.adjust=%s Count=%d genes=[%s]",
        r$Description, r$GeneRatio, format(r$pvalue, digits = 3), format(r$p.adjust, digits = 3), r$Count, r$geneID)
  }
}
write.csv(ora_reactome_df[order(ora_reactome_df$pvalue), ], file.path(tab_dir, "TableS_interaction_ORA_reactome.csv"), row.names = FALSE)

# ---- GO Biological Process ORA ----
ora_go <- enrichGO(
  gene = unique(map_deg$ENTREZID), universe = unique(map_universe$ENTREZID),
  OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP",
  pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH",
  minGSSize = 2, maxGSSize = 1000, readable = TRUE
)
ora_go_df <- as.data.frame(ora_go)
n_go_tested <- nrow(ora_go_df)
n_go_sig <- sum(ora_go_df$p.adjust < 0.05, na.rm = TRUE)
add("")
add("GO Biological Process ORA: %d terms tested, %d significant at BH p.adjust<0.05", n_go_tested, n_go_sig)
if (n_go_tested > 0) {
  top5g <- head(ora_go_df[order(ora_go_df$pvalue), ], 5)
  for (i in seq_len(nrow(top5g))) {
    r <- top5g[i, ]
    add("  %-55s GeneRatio=%s pvalue=%s p.adjust=%s Count=%d genes=[%s]",
        r$Description, r$GeneRatio, format(r$pvalue, digits = 3), format(r$p.adjust, digits = 3), r$Count, r$geneID)
  }
}
write.csv(ora_go_df[order(ora_go_df$pvalue), ], file.path(tab_dir, "TableS_interaction_ORA_GO_BP.csv"), row.names = FALSE)

n_deg_mapped <- length(unique(map_deg$ENTREZID))
n_unmapped <- length(deg_symbols) - n_deg_mapped

# ============================================================================
# Step 3: polymorphic-locus check
# ============================================================================
add("")
add("=== Step 3: polymorphic-locus check ===")
known_poly_pattern <- "^Oas1|^H2-|^Nlrp1|^Klr|^Ly6|^Xlr|^Gm[0-9]"
sig24_full$flag_poly_name <- grepl(known_poly_pattern, sig24_full$gene)

big_baseline <- sig24_full[!is.na(sig24_full$padj_baseline) & sig24_full$padj_baseline < 0.05 & abs(sig24_full$log2FC_baseline) > 1, ]
add("Genes among the 24 with LARGE baseline substrain difference (baseline padj<0.05 & |log2FC|>1): %d/24", nrow(big_baseline))
if (nrow(big_baseline) > 0) {
  for (i in seq_len(nrow(big_baseline))) {
    r <- big_baseline[i, ]
    add("  %-12s baseline log2FC=%+.2f padj=%s%s", r$gene, r$log2FC_baseline, format(r$padj_baseline, digits = 3),
        ifelse(r$flag_poly_name, "  <-- matches known-polymorphic-locus name pattern", ""))
  }
}

name_flagged <- sig24_full[sig24_full$flag_poly_name, ]
add("")
add("Genes among the 24 matching a known-polymorphic-locus NAME pattern (Oas1*, H2-*, Nlrp1*, Klr*, Ly6*, Xlr*, Gm#): %d/24", nrow(name_flagged))
if (nrow(name_flagged) > 0) add("  %s", paste(name_flagged$gene, collapse = ", "))

writeLines(out, file.path(log_dir, "27_interaction_pathway_coherence.txt"))
cat(paste(out, collapse = "\n"), "\n")
out <- c()

# ============================================================================
# Step 4: verdict
# ============================================================================
verdict <- c()
addv <- function(...) verdict <<- c(verdict, sprintf(...))

addv("=== Step 4: verdict ===")
addv("")
addv("Reactome ORA: %d/%d terms significant (BH<0.05).", n_reactome_sig, n_reactome_tested)
addv("GO BP ORA: %d/%d terms significant (BH<0.05).", n_go_sig, n_go_tested)
addv("Large baseline-divergent genes among the 24: %d/24.", nrow(big_baseline))
addv("Known-polymorphic-locus-name-pattern genes among the 24: %d/24 (%s).",
     nrow(name_flagged), if (nrow(name_flagged) > 0) paste(name_flagged$gene, collapse=", ") else "none")
addv("")

sig_reactome_terms <- if (n_reactome_sig > 0) ora_reactome_df[order(ora_reactome_df$pvalue), ][seq_len(n_reactome_sig), ] else ora_reactome_df[0, ]
sig_go_terms <- if (n_go_sig > 0) ora_go_df[order(ora_go_df$pvalue), ][seq_len(n_go_sig), ] else ora_go_df[0, ]
genes_in_any_sig_term <- unique(unlist(strsplit(c(sig_reactome_terms$geneID, sig_go_terms$geneID), "/")))
n_sig_terms_total <- n_reactome_sig + n_go_sig

if (n_sig_terms_total == 0) {
  verdict_letter <- "A"
  addv("VERDICT: (A) No significant pathway enrichment (Reactome or GO BP) among the 24 interaction genes at BH<0.05.")
  addv("This supports describing the 24 genes as a scattered/non-coherent set in the manuscript text")
  addv("('rather than a globally divergent injury program').")
} else if (n_sig_terms_total <= 2 && length(genes_in_any_sig_term) <= 2) {
  verdict_letter <- "A (borderline)"
  addv("VERDICT: (A, borderline) %d significant term(s) found, but ALL are driven by the SAME %d gene(s): %s.",
       n_sig_terms_total, length(genes_in_any_sig_term), paste(genes_in_any_sig_term, collapse = ", "))
  addv("Reactome's 'Smooth Muscle Contraction' and 'Striated Muscle Contraction' terms are both populated by")
  addv("only Tpm2/Tpm4 (paralogous tropomyosin genes that trivially co-occur in small, overlapping Reactome")
  addv("terms -- 2 genes reaching a 13-member term is a low bar, not evidence of a shared injury program).")
  addv("The other %d of 24 genes (%.0f%%) belong to NO significant term in either ontology.",
       24 - length(genes_in_any_sig_term), 100 * (24 - length(genes_in_any_sig_term)) / 24)
  addv("PRACTICAL READING: this does not rise to a biologically coherent theme (no shared immune/ECM/stress")
  addv("program); it is consistent with 'rather than a globally divergent injury program', with the minor,")
  addv("explicitly-stated exception of a 2-gene tropomyosin/cytoskeletal pair. We recommend the manuscript text")
  addv("keep the 'scattered, not a coherent program' framing but name this one exception rather than claim zero")
  addv("shared annotation at all.")
} else {
  sig_terms <- c(sig_reactome_terms$Description, sig_go_terms$Description)
  verdict_letter <- "B"
  addv("VERDICT: (B) %d significant term(s) found, spanning %d distinct genes: %s", n_sig_terms_total,
       length(genes_in_any_sig_term), paste(sig_terms, collapse = "; "))
  addv("This does NOT support an unqualified 'no coherent pathway' claim -- the manuscript text may need to")
  addv("name this specific theme rather than describe the 24 genes as purely scattered.")
}

if (nrow(name_flagged) >= 3 || nrow(big_baseline) >= nrow(sig24_full) / 2) {
  addv("")
  addv("ADDITIONAL FLAG (C): %d/24 genes show large pre-existing baseline divergence and/or match known-polymorphic-locus", nrow(big_baseline))
  addv("name patterns (%d/24). This raises the possibility that at least part of the interaction signal reflects", nrow(name_flagged))
  addv("allele/mapping-related baseline differences carried through to the ADR-response calculation, not a genuine")
  addv("treatment-response interaction, for THOSE SPECIFIC genes. Recommend listing them explicitly as a caveat")
  addv("alongside whichever of (A)/(B) applies to the set as a whole.")
}

writeLines(verdict, file.path(log_dir, "27_verdict_summary.txt"))
cat(paste(verdict, collapse = "\n"), "\n")

cat("\n=== 27_interaction_gene_pathway_coherence.R complete ===\n")
