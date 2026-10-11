############################################################
## nested_cv_v3.R —— 化合物分组的嵌套交叉验证(修正版,取代 nested_cv.R)
## 现名说明:本文件即原先的 rat_grouped_cv.R,v3 命名与原 nested_cv.R 对应。
## 相对原 nested_cv.R 的三处修正:
##   ① 外层折不再随机分配:按【化合物】分组(同化合物不跨折),对照组按实验批次分块;
##      内层调参折同样按化合物分组(真正的 grouped nested CV)。
##   ② 预处理(medianImpute / nearZeroVar)只在训练折内部拟合。
##   ③ 每样本只保留 1 个折外预测(重复取平均)后再算 AUC 与 bootstrap CI。
## 输出: outputs/v2/grouped/  (与训练主流程的输出分开存放)
## 回应审稿人对 nested_cv.R 的三点质疑:
##  (1) 外层折按【化合物】分组(同化合物不跨折),对照组按【实验批次】分块;
##      内层调参同样按化合物分组 = 真正的 grouped nested CV。
##  (2) 预处理(medianImpute / nearZeroVar)只在训练折内部拟合。
##  (3) 每样本只保留一个 OOF 预测(重复取平均)后再算 AUC 与 bootstrap CI。
## 输入: outputs/rat_expr45.rds ; inputs/rat_samples.tsv ; inputs/dm_labels.txt
## 输出: outputs/v2/grouped/
############################################################
suppressMessages({
  library(data.table); library(caret); library(pROC)
  library(randomForest); library(kernlab); library(xgboost)
  library(gbm); library(nnet); library(glmnet); library(rpart); library(doParallel)
})
setwd("D:/MLroute2")
out <- "outputs/v2/grouped"; dir.create(out, showWarnings = FALSE, recursive = TRUE)

r45 <- readRDS("outputs/rat_expr45.rds")
X   <- t(r45$X)
m   <- as.data.table(r45$m)
gsm <- as.character(m$gsm)
y   <- ifelse(m$L1 == "Treat", 1L, 0L)
yf  <- factor(ifelse(y == 1, "Treat", "Control"), levels = c("Control", "Treat"))
n   <- nrow(X)
cat("n =", n, "| Treat", sum(y == 1), "| Control", sum(y == 0), "| p =", ncol(X), "\n")

## ---------- 1. 化合物分组键 ----------
lab  <- fread("inputs/dm_labels.txt", sep = "\t", header = TRUE)
for (cc in c("ARRAY_ID","CHEMICAL","ANIMAL_ID","TOXICOLOGY_EXPERIMENT_DATE"))
  lab[[cc]] <- trimws(as.character(lab[[cc]]))
meta <- fread("inputs/rat_samples.tsv")
meta[, acc := trimws(as.character(acc))]
mi <- match(gsm, meta$acc); stopifnot(!anyNA(mi), !anyDuplicated(mi))
meta <- meta[mi]; stopifnot(all(gsm == meta$acc))   # 显式重排,列序与文件行序本就不同
bc <- paste0(as.character(meta$barcode), ".CEL")
li <- match(bc, lab$ARRAY_ID)
cat(sprintf("CEBS 匹配: %d/%d (%.1f%%)\n", sum(!is.na(li)), n, 100*mean(!is.na(li))))
stopifnot(!anyNA(li))

chem   <- lab$CHEMICAL[li]
animal <- lab$ANIMAL_ID[li]
batch  <- sub("-.*$", "", animal)
grp    <- ifelse(y == 1, paste0("CMP_", chem), paste0("VEH_", batch))
stopifnot(all(chem[y == 1] != ""))
cat("groups:", length(unique(grp)), "| 化合物组:", length(unique(grp[y==1])),
    "| 对照批次组:", length(unique(grp[y==0])), "\n")
fwrite(data.table(gsm, barcode = bc, group = grp, chemical = chem, animal, L1 = m$L1),
       file.path(out, "rat_compound_map.tsv"), sep = "\t")

## ---------- 2. 分组折 ----------
K <- 5; Rrep <- 2
grouped_folds <- function(g, K, seed = NULL) {
  ug <- unique(g); sz <- as.numeric(table(g)[ug])
  if (!is.null(seed)) { set.seed(seed); o <- order(-sz, runif(length(sz))) } else o <- order(-sz)
  fg <- integer(length(ug)); names(fg) <- ug
  fg[o] <- ((seq_along(o) - 1L) %% K) + 1L
  as.integer(fg[g])
}
folds <- lapply(seq_len(Rrep), function(r) grouped_folds(grp, K, seed = 300 + r))
for (r in seq_len(Rrep)) {
  cat("\nrepeat", r, " (fold x class):\n"); print(table(folds[[r]], y))
  bad <- sum(sapply(split(folds[[r]], grp), function(z) length(unique(z)) > 1))
  cat("  跨折组数(须为0):", bad, "\n")
}

## ---------- 3. 逐折拟合 ----------
specs <- list(
  RF    = list(method = "rf",        grid = data.frame(mtry = c(6, 15))),
  SVM   = list(method = "svmRadial", tuneLength = 2),
  GLM   = list(method = "glm"),
  GBM   = list(method = "gbm",       tuneLength = 2),
  KNN   = list(method = "knn",       grid = data.frame(k = c(5, 15))),
  NNET  = list(method = "nnet",      grid = expand.grid(size = 5, decay = 0.1)),
  LASSO = list(method = "glmnet",    tuneLength = 3),
  DT    = list(method = "rpart",     grid = data.frame(cp = c(0.01, 0.001))),
  XGB   = list(method = "xgb_native", nrounds = 100))
models <- names(specs)
P <- setNames(lapply(models, function(z) matrix(NA_real_, Rrep*K, n)), models)

inner_ctrl <- function(tr_glob) {
  g  <- grp[tr_glob]
  ff <- grouped_folds(g, 5, seed = 11)
  ii <- lapply(seq_len(5), function(k) which(ff != k))
  trainControl(method = "cv", number = 5, classProbs = TRUE,
               summaryFunction = twoClassSummary,
               index = ii, indexOut = lapply(ii, function(z) setdiff(seq_along(g), z)),
               allowParallel = FALSE)
}

cl <- makePSOCKcluster(6); registerDoParallel(cl); cat("workers: 6\n")
t0 <- Sys.time(); BT <- list()
for (r in seq_len(Rrep)) for (k in seq_len(K)) {
  teG <- which(folds[[r]] == k); trG <- which(folds[[r]] != k)
  trd <- as.data.frame(X[trG, , drop = FALSE]); pp <- preProcess(trd, method = "medianImpute")
  trd <- predict(pp, trd); ted <- predict(pp, as.data.frame(X[teG, , drop = FALSE]))
  keep <- setdiff(colnames(trd), nearZeroVar(trd)); if (!length(keep)) keep <- colnames(trd)
  dtr <- trd[, keep, drop = FALSE]; dte <- ted[, keep, drop = FALSE]
  dtr$Type <- yf[trG]
  ctrl <- inner_ctrl(trG)
  for (nm in models) {
    sp <- specs[[nm]]
    pr <- tryCatch({
      if (sp$method == "xgb_native") {
        set.seed(1)
        params <- list(objective="binary:logistic", eval_metric="auc", max_depth=3,
                       eta=0.05, subsample=0.8, colsample_bytree=0.8)
        b <- xgb.train(params = params,
                       data = xgb.DMatrix(as.matrix(dtr[, setdiff(colnames(dtr), "Type"), drop = FALSE]),
                                          label = y[trG]),
                       nrounds = sp$nrounds, verbose = 0)
        list(p = predict(b, as.matrix(dte)), best = paste0("nrounds=", sp$nrounds))
      } else {
        a <- list(Type ~ ., data = dtr, method = sp$method, metric = "ROC", trControl = ctrl)
        if (!is.null(sp$grid)) a$tuneGrid <- sp$grid else a$tuneLength <- sp$tuneLength
        if (sp$method == "gbm")  a$verbose <- FALSE
        if (sp$method == "nnet") a$trace <- FALSE
        if (sp$method == "glm")  a$family <- binomial()
        mo <- do.call(train, a)
        list(p = as.numeric(predict(mo, dte, type = "prob")[, "Treat"]),
             best = paste(unlist(mo$bestTune), collapse = "/"))
      }
    }, error = function(e) { cat("ERR", nm, conditionMessage(e), "\n"); NULL })
    if (!is.null(pr)) P[[nm]][(r-1)*K + k, teG] <- pr$p
    BT[[length(BT)+1]] <- data.frame(rep_idx = r, fold = k, model = nm,
                                     best = if (is.null(pr)) NA else pr$best)
  }
  cat(sprintf("  r%d k%d  nte=%d  (%.1f min)\n", r, k, length(teG),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
stopCluster(cl)
fwrite(do.call(rbind, BT), file.path(out, "rat_grouped_bestTune.txt"), sep = "\t")

## ---------- 4. 按样本平均 -> 每样本一个预测 ----------
oof <- sapply(models, function(nm) colMeans(P[[nm]], na.rm = TRUE))
rownames(oof) <- gsm
fwrite(as.data.table(cbind(gsm = gsm, as.data.frame(round(oof, 5)))),
       file.path(out, "rat_grouped_OOF_persample.txt"), sep = "\t")

per_fold <- do.call(rbind, lapply(models, function(nm) do.call(rbind, lapply(seq_len(Rrep), function(r)
  do.call(rbind, lapply(seq_len(K), function(k) {
    teG <- which(folds[[r]] == k); p <- P[[nm]][(r-1)*K + k, teG]
    if (all(is.na(p)) || length(unique(p[!is.na(p)])) < 2) return(NULL)
    data.frame(model = nm, rep_idx = r, fold = k,
               AUC = round(as.numeric(auc(roc(y[teG], p, quiet = TRUE))), 4))
  }))))))
fwrite(per_fold, file.path(out, "rat_grouped_perf_perFold.txt"), sep = "\t")

## 【v3-5】9 分类器的合并 OOF AUC 也用【化合物聚簇】bootstrap,
##   与单基因 CCNB2 的口径一致 —— 理由相同:同化合物样本不独立。
boot <- function(p, yy, gg, B = 2000) {
  ok <- !is.na(p); p <- p[ok]; yy <- yy[ok]; gg <- gg[ok]
  a <- as.numeric(auc(roc(yy, p, quiet = TRUE))); set.seed(99)
  ug <- unique(gg)
  bs <- replicate(B, {
    ix <- unlist(lapply(sample(ug, replace = TRUE), function(g) which(gg == g)))
    if (length(unique(yy[ix])) < 2) return(NA)
    as.numeric(auc(roc(yy[ix], p[ix], quiet = TRUE)))
  })
  c(a, quantile(bs, .025, na.rm = TRUE), quantile(bs, .975, na.rm = TRUE))
}
stat <- do.call(rbind, lapply(models, function(nm) {
  b <- boot(oof[, nm], y, grp); pf <- per_fold[per_fold$model == nm, ]
  data.frame(model = nm,
             groupedCV_AUC_mean = round(mean(pf$AUC, na.rm = TRUE), 3),
             groupedCV_AUC_sd   = round(sd(pf$AUC, na.rm = TRUE), 3),
             pooledOOF_AUC      = round(b[1], 3),
             boot_lo = round(b[2], 3), boot_hi = round(b[3], 3),
             n_fold_used = sum(!is.na(pf$AUC)))
}))
fwrite(stat, file.path(out, "rat_grouped_perf.txt"), sep = "\t")
cat("\n===== compound-grouped CV (rat, 45 genes) =====\n"); print(stat, row.names = FALSE)

## ---------- 5. 单基因 CCNB2(描述性;分组不改变其定义) ----------
gi <- which(toupper(colnames(X)) == "CCNB2")
if (length(gi) == 1) {
  v <- X[, gi]
  rr <- roc(y, v, quiet = TRUE, direction = "<", levels = c(0, 1)); ci <- ci.auc(rr)
  cat(sprintf("\nCCNB2 univariate AUC = %.3f (DeLong 95%%CI %.3f-%.3f)\n", auc(rr), ci[1], ci[3]))
  ug <- unique(grp); set.seed(7)
  idxs <- lapply(seq_len(2000), function(z) unlist(lapply(sample(ug, replace = TRUE),
                                                          function(g) which(grp == g))))
  bs  <- sapply(idxs, function(ix) { if (length(unique(y[ix])) < 2) NA else
    as.numeric(auc(roc(y[ix], v[ix], quiet = TRUE, direction = "<"))) })
  bs2 <- sapply(idxs, function(ix) { if (sum(y[ix]==1) < 2 || sum(y[ix]==0) < 2) NA else
    mean(v[ix][y[ix]==1]) - mean(v[ix][y[ix]==0]) })
  cat(sprintf("CCNB2 AUC, compound-clustered bootstrap 95%%CI = %.3f-%.3f\n",
              quantile(bs, .025, na.rm=TRUE), quantile(bs, .975, na.rm=TRUE)))
  cat(sprintf("CCNB2 log2FC = %+.3f (cluster-boot 95%%CI %.3f to %.3f)\n",
              mean(v[y==1]) - mean(v[y==0]), quantile(bs2, .025, na.rm=TRUE),
              quantile(bs2, .975, na.rm=TRUE)))
}
sink(file.path(out, "sessionInfo.txt")); print(sessionInfo()); sink()
cat("\nALL DONE ->", out, "\n")
