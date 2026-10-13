#!/usr/bin/env Rscript
# =============================================================================
# RNA-seq differential expression analysis — NTH-exposed crucian carp liver
# (Carassius auratus, aligned to the Carassius gibelio reference carGib1.2-hapl.c)
#
# Complete, executable DEG workflow reported in the manuscript:
#   1. Count import           (featureCounts raw count matrix)
#   2. Low-expression filtering
#   3. TMM normalization
#   4. Design matrix
#   5. Dispersion estimation
#   6. Contrast specification
#   7. Differential-expression testing
#
# DEGs are declared at FDR < 0.05 and |log2FC| > 1, where the log2 fold change
# is computed from FPKM values as the log2 ratio of group means, consistent with
# the Methods ("FPKM values were not used as input for edgeR").
#
# Environment: R 3.6.3, edgeR (>= 3.28), run on the CFViSA analysis platform.
# =============================================================================

suppressPackageStartupMessages(library(edgeR))

## ---- File paths (adjust to your layout) ------------------------------------
count_file <- "RNAseq_supplementary/05_Expression/genes.read.txt"   # raw counts
fpkm_file  <- "RNAseq_supplementary/05_Expression/genes.fpkm.txt"   # FPKM matrix
out_dir    <- "RNAseq_supplementary/06_DEG/details"

## ---- 1. Count import -------------------------------------------------------
# genes.read.txt is a gene x sample raw count matrix (featureCounts output):
# first column = gene ID; remaining 9 columns = samples C-1..C-3, H-1..H-3,
# L-1..L-3. Genes retain their reference IDs (LOC numbers or symbols).
counts <- read.delim(count_file, header = TRUE, row.names = 1,
                     check.names = FALSE, stringsAsFactors = FALSE)
fpkm   <- read.delim(fpkm_file,  header = TRUE, row.names = 1,
                     check.names = FALSE, stringsAsFactors = FALSE)

# Sample groups derived from the column-name prefix (C- / L- / H-)
group <- factor(sub("-.*", "", colnames(counts)), levels = c("C", "L", "H"))

stopifnot(ncol(counts) == 9,
          nlevels(group) == 3,
          identical(rownames(counts), rownames(fpkm)),
          identical(colnames(counts), colnames(fpkm)))
cat("Imported", nrow(counts), "genes x", ncol(counts), "samples.\n")

## ---- 2. Low-expression filtering -------------------------------------------
dge  <- DGEList(counts = counts, group = group)
keep <- filterByExpr(dge, group = group)
dge  <- dge[keep, , keep.lib.sizes = FALSE]
cat("Low-expression filtering:", nrow(dge), "genes retained of",
    nrow(counts), ".\n")

## ---- 3. TMM normalization --------------------------------------------------
dge <- calcNormFactors(dge, method = "TMM")
cat("TMM normalization factors:\n")
print(round(dge$samples$norm.factors, 4))

## ---- 4. Design matrix ------------------------------------------------------
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)          # C, L, H
cat("Design matrix:\n")
print(design)

## ---- 5. Dispersion estimation ----------------------------------------------
dge <- estimateDisp(dge, design, robust = TRUE)
cat("Common BCV:", round(sqrt(dge$common.dispersion), 4),
    "| common dispersion:", round(dge$common.dispersion, 4), "\n")

## ---- 6. Contrast specification ---------------------------------------------
contrasts <- makeContrasts(
  L_vs_C = L - C,
  H_vs_C = H - C,
  L_vs_H = L - H,
  levels = design
)
cat("Contrast matrix:\n")
print(contrasts)

## ---- log2 fold change from FPKM group means (per Methods) -------------------
# The statistical test (FDR) comes from edgeR on TMM-normalized counts; the
# log2FC used for the DEG threshold is the log2 ratio of FPKM group means.
grp_mean_fpkm <- sapply(levels(group), function(g) {
  rowMeans(fpkm[rownames(dge), group == g, drop = FALSE])
})
fpkm_log2fc <- function(num, den) {
  log2((grp_mean_fpkm[, num] + 1e-3) / (grp_mean_fpkm[, den] + 1e-3))
}

## ---- 7. Differential-expression testing (quasi-likelihood F-test) -----------
fit <- glmQLFit(dge, design, robust = TRUE)

contrast_map <- list(
  L_vs_C = c(num = "L", den = "C"),
  H_vs_C = c(num = "H", den = "C"),
  L_vs_H = c(num = "L", den = "H")
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

run_contrast <- function(name) {
  qlf <- glmQLFTest(fit, contrast = contrasts[, name])
  res <- topTags(qlf, n = Inf, sort.by = "none")$table

  nm <- contrast_map[[name]]
  l2fc <- fpkm_log2fc(nm["num"], nm["den"])

  out <- data.frame(
    gene       = rownames(res),
    log2FC_fpkm = round(l2fc, 6),
    logCPM     = res$logCPM,
    F          = res$F,
    PValue     = res$PValue,
    FDR        = res$FDR,
    stringsAsFactors = FALSE
  )
  out$direction <- ifelse(out$log2FC_fpkm > 0, "up", "down")

  # DEG: FDR < 0.05 AND |log2FC| > 1 (FPKM-based)
  deg <- out[out$FDR < 0.05 & abs(out$log2FC_fpkm) > 1, , drop = FALSE]
  deg <- deg[order(deg$FDR), ]

  cat(sprintf("%s: %d DEGs (up %d / down %d)\n",
              name, nrow(deg), sum(deg$direction == "up"), sum(deg$direction == "down")))

  write.table(out[order(out$FDR), ],
              file.path(out_dir, paste0(name, ".fpkm.DEG.xls")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(deg,
              file.path(out_dir, paste0(name, ".fpkm.DEG.sig.xls")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(deg)
}

results <- lapply(names(contrast_map), run_contrast)
names(results) <- names(contrast_map)

cat("DEG tables written to:", out_dir, "\n")
