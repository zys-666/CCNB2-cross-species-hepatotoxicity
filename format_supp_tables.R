# format_supp_tables.R : S3 / S4 / S6 formal supplementary tables with header notes
suppressMessages({library(data.table)})
setwd("D:/MLroute2")
out <- "outputs/paper/appendices"; dir.create(out, showWarnings = FALSE)

## ---------- S3 : classifier CV ----------
cv  <- fread("outputs/v2/modelPerf_CV.txt")
ncv <- fread("outputs/v2/nestedCV_perf.txt")
ord <- c("RF", "SVM", "XGB", "GLM", "GBM", "KNN", "NNET", "LASSO", "DT")
s3 <- merge(cv, ncv, by = "Model", all = TRUE)
s3 <- s3[match(ord, Model)]
setnames(s3, c("Model", "CV_AUC_mean", "CV_AUC_sd", "CV_Sens", "CV_Spec",
               "NestedCV_AUC_mean", "NestedCV_AUC_sd", "NestedCV_min", "NestedCV_max"))
f3 <- file.path(out, "Supplementary_Table_S3_classifier_CV.txt")
notes3 <- c(
  "# Supplementary Table S3. Screening-classifier performance (45 network-derived candidate genes, rat cohort",
  "# GSE57815, GEO-verified labels: 1,939 exposed vs 279 controls).",
  "# CV_* = stratified 5-fold x 3-repeated cross-validation (mean +/- SD); Sensitivity/Specificity at the 0.5",
  "# probability cut-off averaged over resamples (native XGBoost was evaluated by a manual stratified 5-fold x 3",
  "# CV; its Sens/Spec come from pooled out-of-fold predictions at the 0.5 cut-off). NestedCV_* = fully nested CV (outer stratified 5-fold x 2 repeats; hyperparameters",
  "# re-tuned by inner 5-fold CV inside each outer training partition), AUC per outer fold.",
  "# Hyperparameters fixed and documented in training_v2.R (seed 123); R version records in sessionInfo.txt.", "")
writeLines(notes3, f3)
s3 <- s3[, lapply(.SD, function(x) if (is.numeric(x)) sprintf("%.3f", x) else x)]
write.table(s3, f3, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", append = TRUE)

## ---------- S4 : 45-gene x 9-model ranking matrix ----------
s4 <- fread("outputs/v2/ranks_all_metrics.txt")
setnames(s4, c("gene", "AUC_treatVsCtrl", "discAUC", "rank_RF", "rank_SVM", "rank_GLM",
               "rank_GBM", "rank_KNN", "rank_NNET", "rank_LASSO", "rank_DT", "rank_XGB",
               "rank_uniAUC", "meanRank"))
setnames(s4, c("gene", "AUC_TreatHigher", "AUC_Discriminative",
               paste0("Rank_", c("RF", "SVM", "GLM", "GBM", "KNN", "NNET", "LASSO", "DT", "XGB")),
               "Rank_UnivariateAUC", "MeanRank"))
s4[, AUC_TreatHigher := round(AUC_TreatHigher, 3)][, AUC_Discriminative := round(AUC_Discriminative, 3)]
s4[, MeanRank := round(MeanRank, 2)]
f4 <- file.path(out, "Supplementary_Table_S4_ranking_matrix.txt")
notes4 <- c(
  "# Supplementary Table S4. Consensus ranking of the 45 candidate genes across nine classifiers (rat cohort,",
  "# GSE57815, corrected labels). Rank_* = position of the gene by permutation-based variable importance",
  "# (DALEX, training-side sample n = 500, mean over permutations) within that classifier (1 = most important).",
  "# Rank_UnivariateAUC = rank by single-gene discrimination; AUC_TreatHigher = univariate ROC AUC under the",
  "# convention 'exposed higher' (direction '<'); AUC_Discriminative = max(AUC, 1-AUC) (direction-agnostic).",
  "# MeanRank = mean of the ten ranks (lower = more reproducibly important). CCNB2 ranked first in all nine",
  "# classifiers and in univariate AUC (MeanRank = 1.0). Sorted by MeanRank.", "")
writeLines(notes4, f4)
write.table(s4, f4, sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE)

## ---------- S6 : LOCO + Youden ----------
loco <- fread("outputs/paper/GSE44783_LOCO.txt")
setnames(loco, c("compound", "n_trt", "LOCO_AUC_CCNB2", "LOCO_AUC_3gene"))
locoM <- rbindlist(list(loco, data.table(compound = "mean", n_trt = NA_integer_,
  LOCO_AUC_CCNB2 = round(mean(loco$LOCO_AUC_CCNB2), 3),
  LOCO_AUC_3gene = round(mean(loco$LOCO_AUC_3gene), 3))))
yd <- fread("outputs/paper/CCNB2_youden.txt")
f6 <- file.path(out, "Supplementary_Table_S6_LOCO_Youden.txt")
notes6 <- c(
  "# Supplementary Table S6. Leave-one-compound-out (LOCO) and Youden-threshold discrimination, mouse cohort GSE44783.",
  "# LOCO: for each of 13 compounds, logistic model trained on all other compounds plus a random half of the 120 vehicle",
  "# controls (5 seeded repeats); AUC evaluated on the held-out compound together with held-out controls. CCNB2 = model on",
  "# CCNB2 alone; g3 = model on Ccnb2/Plk1/Ube2c. Youden block: single-gene CCNB2 at Youden-optimal threshold (AUC with",
  "# DeLong 95% CI, exposed-higher convention); threshold units = log2 RMA expression.", "")
con <- file(f6, "w"); writeLines(notes6, con); close(con)
write.table(locoM, f6, sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE)
con <- file(f6, "a"); writeLines(c("", "--- Youden-optimal discrimination (CCNB2 single gene) ---"), con); close(con)
write.table(yd, f6, sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE)
cat("S6:", f6, "\n")
cat("rows S3:", nrow(s3), "| S4:", nrow(s4), "| LOCO:", nrow(loco), "+mean | Youden:", nrow(yd), "\n")
