#!/usr/bin/env Rscript
# Build DESeqDataSet objects used by all downstream steps.
#
# Sample "A1" was confirmed by the user to be A-ADR1 (Step 0: also corroborated
# by glomerular-purity marker CPM check in 00_qc_report.txt: A-ADR1 has the
# lowest podocyte:tubular marker ratio among AJcl Day-5 ADR samples).
#
# Two sample sets are built for every design:
#   "main"      -> A-ADR1 EXCLUDED (11 samples) -- matches the manuscript's
#                  primary reported analysis (Fig. 6A: AJcl n=2, ByJcl n=3)
#   "sens"      -> A-ADR1 INCLUDED (12 samples) -- used only for the Step 5
#                  sensitivity analysis
#
# Two designs are built for each sample set:
#   group design: ~ group   (group = A_Ctrl/A_ADR/B_Ctrl/B_ADR)
#                 used for the three pairwise contrasts (baseline B vs A,
#                 A: ADR vs Ctrl, B: ADR vs Ctrl) via contrast=c("group",...)
#   interaction design: ~ substrain + treatment + substrain:treatment
#                 used for the substrain:treatment interaction term (Step 3)
#
# Both designs fit the same full-rank model to the same data (same residual
# df), so results for equivalent contrasts are numerically identical; the
# group design simply gives unambiguous, easy-to-audit contrast= calls for
# the pairwise comparisons DESeq2 recommends for 2x2 factorial data.

suppressMessages(library(DESeq2))
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

make_dds_pair <- function(samples, label) {
  coldata <- build_coldata(samples)
  coldata$group <- factor(paste(coldata$substrain, coldata$treatment, sep = "_"),
                           levels = c("A_Ctrl", "B_Ctrl", "A_ADR", "B_ADR"))
  cts <- as.matrix(counts[, samples])

  dds_group <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~group)
  dds_group <- dds_group[rowSums(counts(dds_group)) >= 10, ]
  dds_group <- DESeq(dds_group)

  dds_int <- DESeqDataSetFromMatrix(countData = cts, colData = coldata,
                                     design = ~ substrain + treatment + substrain:treatment)
  dds_int <- dds_int[rowSums(counts(dds_int)) >= 10, ]
  dds_int <- DESeq(dds_int)

  saveRDS(dds_group, file.path(outdir, paste0("tables/dds_group_", label, ".rds")))
  saveRDS(dds_int,   file.path(outdir, paste0("tables/dds_int_", label, ".rds")))

  cat("=== ", label, " (n=", length(samples), ") ===\n", sep = "")
  cat("Samples: ", paste(samples, collapse = ", "), "\n")
  cat("Genes after prefilter (rowSums>=10), group design: ", nrow(dds_group), "\n")
  cat("Genes after prefilter (rowSums>=10), interaction design: ", nrow(dds_int), "\n")
  cat("resultsNames (group): ", paste(resultsNames(dds_group), collapse=", "), "\n")
  cat("resultsNames (interaction): ", paste(resultsNames(dds_int), collapse=", "), "\n\n")
}

all_samples <- colnames(counts)
main_samples <- setdiff(all_samples, "A-ADR1")   # primary analysis: excludes A1
sens_samples <- all_samples                       # sensitivity analysis: includes A1

log_file <- file.path(outdir, "logs/01_build_dds_log.txt")
sink(log_file, split = TRUE)
make_dds_pair(main_samples, "main")
make_dds_pair(sens_samples, "sens")
sink()

cat(readLines(log_file), sep = "\n")
