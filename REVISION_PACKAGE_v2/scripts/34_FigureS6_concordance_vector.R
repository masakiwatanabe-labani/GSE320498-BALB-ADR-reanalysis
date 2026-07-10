#!/usr/bin/env Rscript
# Vector-PDF version of Figure S6 (ADR-response concordance scatter,
# AJcl vs ByJcl). Only a PNG existed for this plot anywhere in the project
# (flagged explicitly in submission_supplementary/SUPPLEMENTARY_MANIFEST.txt).
# Reuses the already-computed merged log2FC table (TableS_ADR_response_A_vs_
# B_merged.csv, from 24_manuscript_values_extraction.R, itself a re-derivation
# of 04_interaction_and_crosscomparison.R's original scatter) -- no new
# statistics, same numbers (Pearson r=0.535, Spearman rho=0.425, etc.),
# just re-plotted with a vector PDF device and submission-quality fonts.

suppressMessages(library(ggplot2))

indir <- "/usr/local/jupyter/ADR_BALB"
outdir_main <- file.path(indir, "outputs")
outdir_v2 <- file.path(outdir_main, "REVISION_PACKAGE/REVISION_PACKAGE_v2")
fig_dir <- file.path(outdir_v2, "figures")

pt_mm <- function(pt) pt / 2.845276

m <- read.csv(file.path(outdir_v2, "tables/TableS_ADR_response_A_vs_B_merged.csv"), stringsAsFactors = FALSE)

pear <- cor.test(m$log2FC_A, m$log2FC_B, method = "pearson")
spear <- cor.test(m$log2FC_A, m$log2FC_B, method = "spearman", exact = FALSE)
fit <- lm(log2FC_B ~ log2FC_A, data = m)
fit_sum <- summary(fit)

sig_either <- (!is.na(m$padj_A) & m$padj_A < 0.05) | (!is.na(m$padj_B) & m$padj_B < 0.05)

check <- c(
  sprintf("n genes = %d", nrow(m)),
  sprintf("Pearson r = %.4f (expected 0.535): %s", pear$estimate, ifelse(round(pear$estimate,3)==0.535,"MATCH","*** CHECK ***")),
  sprintf("Spearman rho = %.4f (expected 0.425): %s", spear$estimate, ifelse(round(spear$estimate,3)==0.425,"MATCH","*** CHECK ***")),
  sprintf("OLS slope = %.4f, R2 = %.4f", coef(fit)[2], fit_sum$r.squared)
)
cat(paste(check, collapse = "\n"), "\n\n")

p <- ggplot(m, aes(x = log2FC_A, y = log2FC_B)) +
  geom_point(aes(color = sig_either), size = 1.5, alpha = 0.55) +
  scale_color_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "grey70"),
                      labels = c(`TRUE` = "padj<0.05 in A or B response", `FALSE` = "n.s. in both"),
                      name = NULL) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "steelblue") +
  geom_hline(yintercept = 0, color = "grey85") + geom_vline(xintercept = 0, color = "grey85") +
  annotate("label", x = -Inf, y = Inf, hjust = -0.03, vjust = 1.1,
           label = sprintf("Pearson r = %.2f, p = %.1e\nSpearman rho = %.2f\nslope = %.2f, R2 = %.2f",
                            pear$estimate, pear$p.value, spear$estimate, coef(fit)[2], fit_sum$r.squared),
           size = pt_mm(9.5), fontface = "bold", label.size = 0.3) +
  labs(title = "Figure S6. ADR response concordance between substrains",
       subtitle = "x = log2FC (A: ADR vs Ctrl), y = log2FC (B: ADR vs Ctrl). Dashed line = identity (slope=1); solid line = OLS fit.",
       x = "log2FC, AJcl ADR vs Ctrl", y = "log2FC, ByJcl ADR vs Ctrl") +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(size = 13.5, face = "bold"),
    plot.subtitle = element_text(size = 9.5),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 10.5),
    legend.text = element_text(size = 9.5),
    legend.position = "bottom"
  )

for (ext in c("pdf", "png")) {
  fp <- file.path(fig_dir, paste0("FigureS6_concordance.", ext))
  ggsave(fp, plot = p, width = 9, height = 7.5, device = if (ext == "pdf") cairo_pdf else ext, dpi = 300)
  cat("Saved:", fp, "\n")
}

cat("\n=== 34_FigureS6_concordance_vector.R complete ===\n")
