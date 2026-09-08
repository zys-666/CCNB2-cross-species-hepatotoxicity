# 02_analysis.R : direction check + frozen transfer under corrected labels
suppressMessages({library(data.table); library(pROC); library(randomForest); library(xgboost); library(glmnet)})
setwd("D:/MLroute2")
out <- function(...) { cat(...); cat("\n") }
sink_file <- "outputs/results_stage1.txt"
sink(sink_file, split = TRUE)

# ---------- load rat ----------
r45 <- readRDS("outputs/rat_expr45.rds")
X45 <- r45$X; gsm <- r45$gsm; m <- r45$m
g3 <- rownames(X45)[toupper(rownames(X45)) %in% toupper(c("Ccnb2","Plk1","Ube2c"))]
X3 <- X45[g3, , drop = FALSE]
lab <- m$L1   # corrected label (all dose-0 = Control)
stopifnot(all(!is.na(lab)))
yrat <- ifelse(lab == "Treat", 1, 0)
stopifnot(length(yrat) == ncol(X45))

# sensitivity: L3 (3 vehicles only)
lab3 <- m$L3
keep3 <- !is.na(lab3)
yrat3 <- ifelse(lab3[keep3] == "Treat", 1, 0)
X3_3 <- X3[, keep3, drop = FALSE]; X45_3 <- X45[, keep3, drop = FALSE]

# ---------- load mouse ----------
mg <- fread("inputs/mouse_expr.csv", sep = ",", header = TRUE, check.names = FALSE)
mgenes <- mg[[1]]
M <- as.matrix(mg[, -1, with = FALSE]); rownames(M) <- mgenes
grp <- fread("inputs/mouse_group.csv")
setnames(grp, c("SampleID", "group"))
stopifnot(all(colnames(M) %in% grp$SampleID))
ymouse <- ifelse(grp$group[match(colnames(M), grp$SampleID)] == "Treat", 1, 0)
gm3 <- rownames(M)[toupper(rownames(M)) %in% toupper(g3)]
cat("mouse: n=", ncol(M), " control:", sum(ymouse==0), " treat:", sum(ymouse==1), "\n")
cat("mouse genes found for 3-gene set:", gm3, "\n")
shared45 <- intersect(toupper(rownames(X45)), toupper(rownames(M)))
cat("mouse: of 45 rat genes, found:", length(shared45), "\n")
M3 <- M[toupper(rownames(M)) %in% toupper(g3), , drop = FALSE]

# ============================================================
# 1. DIRECTION: mean(Treat) - mean(Control) per gene, log2 scale
# ============================================================
dir_tab <- function(X, y, subset_name, genes) {
  t1 <- rowMeans(X[, y == 1, drop = FALSE]); t0 <- rowMeans(X[, y == 0, drop = FALSE])
  d <- t1 - t0
  p <- sapply(seq_len(nrow(X)), function(i) {
    a <- X[i, y == 1]; b <- X[i, y == 0]
    tryCatch(t.test(a, b)$p.value, error = function(e) NA)
  })
  data.frame(genes = genes, mean_Treat = round(t1, 3), mean_Control = round(t0, 3),
             delta = round(d, 3), p = signif(p, 3), label = subset_name)
}
cat("\n========== DIRECTION (delta = meanTreat - meanControl, log2) ==========\n")
cat("\n-- rat L1 (all dose0 control; Treat=1837 vs Control=381) --\n")
print(dir_tab(X3, yrat, "rat_L1", g3))
cat("\n-- rat L3 sensitivity (3 vehicles only; Treat=", sum(yrat3==1), " vs Control=", sum(yrat3==0), ") --\n")
print(dir_tab(X3_3, yrat3, "rat_L3", g3))
cat("\n-- mouse (GSE44783; Treat=330 vs Control=120) --\n")
print(dir_tab(M3, ymouse, "mouse", gm3))

# ============================================================
# 2. FROZEN TRANSFER
# ============================================================
auc_ci <- function(pred, y) {
  r <- roc(y, pred, quiet = TRUE, ci = TRUE)
  c(auc = as.numeric(auc(r)), lo = as.numeric(ci.auc(r)[1]), hi = as.numeric(ci.auc(r)[3]))
}
tr_glm <- function(Xtr, ytr, Xte) { m0 <- glm.fit(t(Xtr), ytr, family = binomial()); plogis(as.numeric(t(Xte) %*% m0$coefficients)) }
# NB: glm.fit expects n x p matrix; X columns = samples -> t()

run_frozen <- function(Xtr, ytr, Xte, yte, tag) {
  res <- list()
  # logistic (unregularized)
  Xtr_t <- t(Xtr); Xte_t <- t(Xte)
  fit <- glm.fit(Xtr_t, ytr, family = binomial())
  pr <- plogis(as.numeric(Xte_t %*% fit$coefficients))
  res$glm <- auc_ci(pr, yte)
  # randomForest
  set.seed(1); rf <- randomForest(x = Xtr_t, y = factor(ytr), ntree = 500)
  res$rf <- auc_ci(predict(rf, Xte_t, type = "prob")[, "1"], yte)
  # xgboost (manuscript-like fixed params)
  set.seed(1)
  xb <- xgb.train(params = list(objective = "binary:logistic", eval_metric = "auc",
                                max_depth = 3, eta = 0.05, subsample = 0.8, colsample_bytree = 0.8),
                  data = xgb.DMatrix(Xtr_t, label = ytr), nrounds = 100, verbose = 0)
  res$xgb <- auc_ci(predict(xb, Xte_t), yte)
  for (nm in names(res)) cat(sprintf("  %-4s AUC=%.3f (95%%CI %.3f-%.3f)\n", nm, res[[nm]]["auc"], res[[nm]]["lo"], res[[nm]]["hi"]))
  invisible(res)
}

cat("\n========== FROZEN TRANSFER (trained frozen, no retraining on target) ==========\n")

cat("\n-- A. rat(TRAIN, L1, 3 genes) -> mouse(TEST) --\n")
ra <- run_frozen(X3, yrat, M3, ymouse, "rat3->mouse"); 
cat("-- A2. rat L3-trained (3 genes) -> mouse --\n")
run_frozen(X3_3, yrat3, M3, ymouse, "rat3L3->mouse")

cat("\n-- B. mouse(TRAIN, 3 genes) -> rat(TEST, L1) --\n")
rb <- run_frozen(M3, ymouse, X3, yrat, "mouse->rat3")

# 45-gene frozen transfer (only shared genes between rat and mouse)
cat("\n-- shared genes between rat45 and mouse:", length(shared45), "--\n")
X45s <- X45[toupper(rownames(X45)) %in% shared45, , drop = FALSE]
Ms <- M[toupper(rownames(M)) %in% shared45, , drop = FALSE]
ms <- toupper(rownames(Ms))
X45s <- X45s[match(shared45, toupper(rownames(X45s))), , drop = FALSE]
Ms <- Ms[match(shared45, ms), , drop = FALSE]
cat("-- C. rat(TRAIN, L1, shared45) -> mouse --\n")
rc <- run_frozen(X45s, yrat, Ms, ymouse, "rat45->mouse")
cat("-- D. mouse(TRAIN, shared45) -> rat(TEST L1) --\n")
rd <- run_frozen(Ms, ymouse, X45s, yrat, "mouse45->rat")

cat("\n========== done ==========\n")
sink()
cat("results written to", normalizePath(sink_file), "\n")
