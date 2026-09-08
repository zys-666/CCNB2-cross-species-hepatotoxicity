# nested_cv.R : proper nested CV for the 9 screening classifiers (45 genes, rat L1 labels)
# outer: stratified 5-fold x 2 repeats; inner: 5-fold tuning inside each outer training part
suppressMessages({
  library(caret); library(pROC); library(randomForest); library(kernlab); library(xgboost)
  library(gbm); library(nnet); library(glmnet); library(rpart); library(data.table)
  library(doParallel)
})
setwd("D:/MLroute2")
out <- "outputs/v2"
r45 <- readRDS("outputs/rat_expr45.rds")
X <- t(r45$X)                       # n x 45
y <- ifelse(r45$m$L1 == "Treat", 1, 0)
yf <- factor(ifelse(y == 1, "Treat", "Control"), levels = c("Control", "Treat"))
dat <- as.data.frame(X); dat$Type <- yf
n <- nrow(dat)
cat("n =", n, "| Treat", sum(y), "Control", sum(y == 0), "| p =", ncol(X), "\n")

inner_ctrl <- trainControl(method = "cv", number = 5, classProbs = TRUE,
                           summaryFunction = twoClassSummary, allowParallel = FALSE)
models_spec <- list(
  RF    = list(method = "rf",        grid = data.frame(mtry = c(2, 7, 15))),
  SVM   = list(method = "svmRadial", tuneLength = 3),
  GLM   = list(method = "glm"),
  GBM   = list(method = "gbm",       tuneLength = 3),
  KNN   = list(method = "knn",       grid = data.frame(k = c(5, 9, 15))),
  NNET  = list(method = "nnet",      grid = expand.grid(size = c(3, 5), decay = c(0.01, 0.1))),
  LASSO = list(method = "glmnet",    tuneLength = 3),
  DT    = list(method = "rpart",     tuneLength = 3),
  XGB   = list(method = "xgb_native"))

Rrep <- 2; K <- 5
folds <- lapply(seq_len(Rrep), function(r) {
  set.seed(300 + r)
  sample(rep(seq_len(K), length.out = n))
})
foldOf <- unsplit(folds, rep(seq_len(Rrep), each = n))  # not used directly

fit_predict <- function(spec, tr_idx, te_idx) {
  tr <- dat[tr_idx, ]; te <- dat[te_idx, setdiff(colnames(dat), "Type")]
  if (spec$method == "xgb_native") {
    # inner: choose nrounds by 5-fold xgb.cv on the outer-train part
    set.seed(1)
    cvr <- xgb.cv(params = list(objective = "binary:logistic", eval_metric = "auc",
                                max_depth = 3, eta = 0.05, subsample = 0.8, colsample_bytree = 0.8),
                  data = xgb.DMatrix(as.matrix(tr[, -ncol(tr)]), label = y[tr_idx]),
                  nrounds = 200, nfold = 5, prediction = FALSE, verbose = 0)
    best <- which.max(cvr$evaluation_log$test_auc_mean)
    b <- xgb.train(params = list(objective = "binary:logistic", eval_metric = "auc",
                                 max_depth = 3, eta = 0.05, subsample = 0.8, colsample_bytree = 0.8),
                   data = xgb.DMatrix(as.matrix(tr[, -ncol(tr)]), label = y[tr_idx]),
                   nrounds = best, verbose = 0)
    return(predict(b, as.matrix(te)))
  }
  args <- list(Type ~ ., data = tr, method = spec$method, metric = "ROC", trControl = inner_ctrl)
  if (!is.null(spec$grid)) args$tuneGrid <- spec$grid else args$tuneLength <- spec$tuneLength
  if (spec$method == "gbm") args$verbose <- FALSE
  if (spec$method == "nnet") args$trace <- FALSE
  if (spec$method == "glm") args$family <- binomial()
  m <- do.call(train, args)
  as.numeric(predict(m, te, type = "prob")[, "Treat"])
}

cl <- makePSOCKcluster(6); registerDoParallel(cl)
cat("parallel workers: 6\n")
res <- foreach(r = seq_len(Rrep), .combine = rbind,
               .packages = c("caret", "pROC", "randomForest", "kernlab", "xgboost",
                             "gbm", "nnet", "glmnet", "rpart")) %dopar% {
  f <- folds[[r]]
  do.call(rbind, lapply(seq_len(K), function(k) {
    te <- f == k; tr <- !te
    sapply(names(models_spec), function(nm) {
      p <- tryCatch(fit_predict(models_spec[[nm]], which(tr), which(te)),
                    error = function(e) { cat("ERR", nm, conditionMessage(e), "\n"); rep(NA, sum(te)) })
      if (all(is.na(p))) return(NA)
      as.numeric(auc(roc(y[te], p, quiet = TRUE)))
    })
  }))
}
stopCluster(cl)
colnames(res) <- names(models_spec)
cat("\nNested CV AUC per outer fold (rows = fold x repeat):\n")
print(round(res, 3))
summ <- data.frame(Model = colnames(res),
                   NestedAUC_mean = round(colMeans(res, na.rm = TRUE), 3),
                   NestedAUC_sd   = round(apply(res, 2, sd, na.rm = TRUE), 3),
                   NestedAUC_min  = round(apply(res, 2, min, na.rm = TRUE), 3),
                   NestedAUC_max  = round(apply(res, 2, max, na.rm = TRUE), 3))
write.table(summ, file.path(out, "nestedCV_perf.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n===== nested CV summary =====\n"); print(summ)
saveRDS(list(res = res, summ = summ), file.path(out, "nestedCV_results.rds"))
cat("saved nestedCV_perf.txt & nestedCV_results.rds\n")
