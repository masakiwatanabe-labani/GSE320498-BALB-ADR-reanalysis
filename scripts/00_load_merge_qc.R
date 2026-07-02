#!/usr/bin/env Rscript
# Step 0: load, merge, QC the GSE320498 gene count matrices
# BALB/cAJcl (A) vs BALB/cByJcl (B) adriamycin nephropathy RNA-seq re-analysis
#
# Inputs (read-only, never modified):
#   GSE320498_BL_gene_counts (1).tsv   -> A-Ctrl1/2/3, B-Ctrl1/2/3
#   GSE320498_ADR_gene_counts (1).tsv  -> A-ADR1/2/3, B-ADR1/2/3
#
# Output: outputs/tables/00_merged_counts.tsv, outputs/logs/00_qc_report.txt

set.seed(20260220)  # matches seed already reported in the manuscript methods (image analysis); kept for consistency across all scripts in this project

sessionInfo_path <- "/usr/local/jupyter/ADR_BALB/outputs/logs/00_sessionInfo.txt"

indir <- "/usr/local/jupyter/ADR_BALB"
bl_file  <- file.path(indir, "GSE320498_BL_gene_counts (1).tsv")
adr_file <- file.path(indir, "GSE320498_ADR_gene_counts (1).tsv")

bl  <- read.delim(bl_file,  check.names = FALSE, stringsAsFactors = FALSE)
adr <- read.delim(adr_file, check.names = FALSE, stringsAsFactors = FALSE)

log_lines <- c()
add_log <- function(...) log_lines <<- c(log_lines, paste0(...))

add_log("=== Step 0 QC report ===")
add_log("BL file dim: ", nrow(bl), " genes x ", ncol(bl), " cols; colnames: ", paste(colnames(bl), collapse=", "))
add_log("ADR file dim: ", nrow(adr), " genes x ", ncol(adr), " cols; colnames: ", paste(colnames(adr), collapse=", "))

# Duplicate Geneid check within each file
dup_bl  <- bl$Geneid[duplicated(bl$Geneid)]
dup_adr <- adr$Geneid[duplicated(adr$Geneid)]
add_log("Duplicate Geneid count in BL file: ", length(dup_bl))
if (length(dup_bl) > 0) add_log("  duplicated BL genes: ", paste(unique(dup_bl), collapse=", "))
add_log("Duplicate Geneid count in ADR file: ", length(dup_adr))
if (length(dup_adr) > 0) add_log("  duplicated ADR genes: ", paste(unique(dup_adr), collapse=", "))

# Geneid set comparison
only_bl  <- setdiff(bl$Geneid, adr$Geneid)
only_adr <- setdiff(adr$Geneid, bl$Geneid)
add_log("Geneids only in BL file (not in ADR): ", length(only_bl))
add_log("Geneids only in ADR file (not in BL): ", length(only_adr))

# Merge by Geneid (inner join on identical gene sets; report if any dropped)
merged <- merge(bl, adr, by = "Geneid", all = FALSE)
add_log("Merged dim (inner join): ", nrow(merged), " genes x ", ncol(merged), " cols")

# NA check
na_count <- sum(is.na(merged))
add_log("Total NA cells in merged matrix: ", na_count)

# Reorder columns explicitly
sample_order <- c("A-Ctrl1","A-Ctrl2","A-Ctrl3","B-Ctrl1","B-Ctrl2","B-Ctrl3",
                   "A-ADR1","A-ADR2","A-ADR3","B-ADR1","B-ADR2","B-ADR3")
stopifnot(all(sample_order %in% colnames(merged)))
merged <- merged[, c("Geneid", sample_order)]
add_log("Final column order: ", paste(colnames(merged), collapse=", "))
add_log("Final dim: ", nrow(merged), " genes x ", ncol(merged)-1, " samples")

# Basic library size / non-numeric check
counts_only <- merged[, sample_order]
all_numeric <- all(sapply(counts_only, is.numeric))
add_log("All count columns numeric: ", all_numeric)
add_log("Any negative counts: ", any(counts_only < 0))
libsizes <- colSums(counts_only)
add_log("Library sizes (colSums):")
for (s in names(libsizes)) add_log("  ", s, ": ", format(libsizes[s], big.mark=","))

# Zero-count genes across all 12 samples
all_zero <- rowSums(counts_only) == 0
add_log("Genes with 0 counts across all 12 samples: ", sum(all_zero))

# Row names = Geneid, write merged matrix
rownames(merged) <- merged$Geneid
write.table(merged, file.path(indir, "outputs/tables/00_merged_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Purity check to help identify sample "A1" flagged in the manuscript
# (one AJcl Day-5 ADR sample reportedly showed reduced glomerular purity)
podo_markers   <- c("Nphs1","Nphs2","Wt1","Podxl","Synpo","Nes")
tubule_markers <- c("Lrp2","Slc34a1","Aqp1","Slc12a1","Umod","Aqp2")

cpm <- sweep(counts_only, 2, colSums(counts_only), FUN = "/") * 1e6
rownames(cpm) <- merged$Geneid

get_rows <- function(genes) cpm[rownames(cpm) %in% genes, , drop = FALSE]

add_log("")
add_log("=== Purity marker CPM check (for identifying sample A1) ===")
add_log("Podocyte marker CPM (A-ADR1/2/3, B-ADR1/2/3):")
pm <- get_rows(podo_markers)
for (g in rownames(pm)) add_log("  ", g, ": ", paste(names(pm), round(as.numeric(pm[g,]),1), sep="=", collapse=", "))
add_log("Tubular marker CPM (A-ADR1/2/3, B-ADR1/2/3):")
tm <- get_rows(tubule_markers)
for (g in rownames(tm)) add_log("  ", g, ": ", paste(names(tm), round(as.numeric(tm[g,]),1), sep="=", collapse=", "))

# summed z-scores across ADR samples only, to flag outlier
adr_samples <- c("A-ADR1","A-ADR2","A-ADR3","B-ADR1","B-ADR2","B-ADR3")
podo_sum <- colSums(get_rows(podo_markers)[, adr_samples])
tubule_sum <- colSums(get_rows(tubule_markers)[, adr_samples])
add_log("")
add_log("Summed podocyte-marker CPM per ADR sample: ", paste(adr_samples, round(podo_sum,1), sep="=", collapse=", "))
add_log("Summed tubular-marker CPM per ADR sample: ", paste(adr_samples, round(tubule_sum,1), sep="=", collapse=", "))
add_log("Podocyte:Tubular marker ratio per ADR sample: ", paste(adr_samples, round(podo_sum/tubule_sum,2), sep="=", collapse=", "))

writeLines(log_lines, file.path(indir, "outputs/logs/00_qc_report.txt"))
cat(paste(log_lines, collapse = "\n"), "\n")

writeLines(capture.output(sessionInfo()), sessionInfo_path)
