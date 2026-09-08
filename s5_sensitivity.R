# s5_sensitivity.R : Supplementary Table S5 — mouse GSE44783 CCNB2 sensitivity strata (time x sex) + per compound
suppressMessages({library(data.table); library(pROC)})
setwd("D:/MLroute2")
st  <- fread("inputs/mouse_samples.tsv")
grp <- fread("inputs/mouse_group.csv"); setnames(grp, c("SampleID", "group"))
mg  <- fread("inputs/mouse_expr.csv", sep = ",", header = TRUE, check.names = FALSE)
M   <- as.matrix(mg[, -1, with = FALSE]); rownames(M) <- mg[[1]]
v   <- M[toupper(rownames(M)) == "CCNB2", ]            # log2 CCNB2 per sample
sid <- colnames(M)
y   <- ifelse(grp$group[match(sid, grp$SampleID)] == "Treat", 1, 0)
meta <- st[match(sid, st$acc)]
comp <- meta$treatment; tm <- meta$time; sex <- tolower(meta$gender)
cat("unique treatment:\n"); print(sort(unique(comp[y == 1])))

rowf <- function(label, ii) {
  tt <- y[ii] == 1; cc <- !tt
  nT <- sum(tt); nC <- sum(cc)
  if (nT == 0 || nC == 0) return(NULL)
  mT <- mean(v[ii][tt]); mC <- mean(v[ii][cc]); d <- mT - mC
  tt1 <- t.test(v[ii][tt], v[ii][cc])
  ci <- tt1$conf.int
  r <- roc(y[ii], v[ii], quiet = TRUE, direction = "<", ci = TRUE)
  data.frame(Stratum = label, n_Treat = nT, n_Control = nC,
             Treat_mean_log2 = round(mT, 3), Ctrl_mean_log2 = round(mC, 3),
             Delta_log2 = round(d, 3), CI_lo = round(ci[1], 3), CI_hi = round(ci[2], 3),
             p_Welch = format.pval(tt1$p.value, digits = 3, eps = 1e-300),
             AUC = round(as.numeric(auc(r)), 3),
             AUC_CI_lo = round(ci.auc(r)[1], 3), AUC_CI_hi = round(ci.auc(r)[3], 3))
}
rows <- list(
  rowf("All (330 vs 120)", seq_along(y)),
  rowf("Time: 4 days", which(tm == "4 days")),
  rowf("Time: 15 days", which(tm == "15 days")),
  rowf("Sex: male", which(sex == "male")),
  rowf("Sex: female", which(sex == "female")),
  rowf("4 days, male", which(tm == "4 days" & sex == "male")),
  rowf("4 days, female", which(tm == "4 days" & sex == "female")),
  rowf("15 days, male", which(tm == "15 days" & sex == "male")),
  rowf("15 days, female", which(tm == "15 days" & sex == "female")))
for (c in sort(unique(comp[y == 1]))) rows[[length(rows) + 1]] <-
  rowf(paste0("Compound: ", c, " (vs all 120 controls)"),
       which((!is.na(comp) & comp == c) | y == 0))
tab <- do.call(rbind, rows)
names(tab)[names(tab) %in% c("AUC", "AUC_CI_lo", "AUC_CI_hi")] <-
  c("AUC_TreatHigher", "AUC_TreatHigher_CI_lo", "AUC_TreatHigher_CI_hi")
f <- "outputs/paper/Supplementary_Table_S5_mouse_sensitivity.txt"
notes <- c(
  "# Supplementary Table S5. CCNB2 single-gene sensitivity analyses, mouse cohort GSE44783.",
  "# Expression: RMA log2 (Affymetrix Mouse Genome 430 2.0). Stratum-specific Welch t-test on log2 expression;",
  "# AUC_TreatHigher = ROC AUC under the convention that exposed samples have higher CCNB2 (direction '<');",
  "# AUC < 0.5 therefore indicates the opposite (repressed) direction, as seen for the non-inducing",
  "# non-hepatocarcinogens cyproterone acetate and propranolol. DeLong 95% CI. Per-compound rows compare each",
  "# compound (n_Treat) against all 120 vehicle controls; time/sex strata use stratum-matched controls.",
  "")
writeLines(notes, f)
write.table(tab, f, sep = "\t", quote = FALSE, row.names = FALSE, na = "", append = TRUE)
cat("written", f, "| rows:", nrow(tab), "\n")
print(tab, row.names = FALSE)
