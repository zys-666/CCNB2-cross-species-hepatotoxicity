# 03_mouse_cv.R : honest mouse-side evaluation (replaces in-sample 0.998)
suppressMessages({library(data.table); library(pROC); library(randomForest); library(xgboost)})
setwd("D:/MLroute2")
sink("outputs/results_mouse_cv.txt", split = TRUE)

M3 <- readRDS("outputs/rat_expr3.rds")  # placeholder not used; reload mouse below
grp <- fread("inputs/mouse_group.csv"); setnames(grp, c("SampleID","group"))
st <- fread("inputs/mouse_samples.tsv")  # geo parse: acc, treatment, time, gender
mg <- fread("inputs/mouse_expr.csv", sep=",", header=TRUE, check.names=FALSE)
mgenes <- mg[[1]]
M <- as.matrix(mg[, -1, with=FALSE]); rownames(M) <- mgenes
g3 <- c("Ccnb2","Plk1","Ube2c")
M3 <- M[toupper(rownames(M)) %in% toupper(g3), , drop=FALSE]
y <- ifelse(grp$group[match(colnames(M), grp$SampleID)] == "Treat", 1, 0)
trt <- st$treatment[match(colnames(M), st$acc)]
cat("n=", length(y), " Treat:", sum(y), " Control:", sum(y==0), "\n")
X <- t(M3)  # n x 3

set.seed(1)
auc1 <- function(p, yy) as.numeric(auc(roc(yy, p, quiet=TRUE)))

# ---- 1. in-sample refit (reproduce the manuscript 0.998) ----
fit_glm <- glm(y ~ ., data=as.data.frame(X), family=binomial)
cat("\nIN-SAMPLE glm AUC:", round(auc1(predict(fit_glm, type="response"), y), 4), "\n")
set.seed(1); fit_rf <- randomForest(x=X, y=factor(y), ntree=500)
cat("IN-SAMPLE rf  AUC:", round(auc1(predict(fit_rf, type="prob")[,"1"], y), 4), "\n")
set.seed(1)
fit_xgb <- xgb.train(params=list(objective="binary:logistic", eval_metric="auc",
                                 max_depth=6, eta=0.3),  # validation-script defaults
                     data=xgb.DMatrix(X, label=y), nrounds=50, verbose=0)
cat("IN-SAMPLE xgb AUC:", round(auc1(predict(fit_xgb, X), y), 4), "\n")

# ---- 2. repeated 5-fold CV (refit inside folds) ----
rep_cv <- function(k=5, R=5, seed=42) {
  Pg <- Pr <- Px <- rep(NA, length(y))
  for (r in seq_len(R)) {
    set.seed(seed + r)
    f <- sample(rep(seq_len(k), length.out=length(y)))
    for (i in seq_len(k)) {
      te <- f == i
      d <- as.data.frame(X)
      m1 <- glm(y[!te] ~ ., data=d[!te,,drop=FALSE], family=binomial)
      Pg[te] <- predict(m1, d[te,,drop=FALSE], type="response")
      m2 <- randomForest(x=X[!te,,drop=FALSE], y=factor(y[!te]), ntree=300)
      Pr[te] <- predict(m2, X[te,,drop=FALSE], type="prob")[,"1"]
      set.seed(seed+r+i)
      m3 <- xgb.train(params=list(objective="binary:logistic", eval_metric="auc",
                                  max_depth=3, eta=0.05, subsample=0.8, colsample_bytree=0.8),
                      data=xgb.DMatrix(X[!te,,drop=FALSE], label=y[!te]), nrounds=100, verbose=0)
      Px[te] <- predict(m3, X[te,,drop=FALSE])
    }
  }
  list(glm=Pg, rf=Pr, xgb=Px)
}
cv <- rep_cv()
cat("\nREPEATED 5-FOLD CV (refit, 5 repeats) pooled AUC:\n")
for (nm in names(cv)) {
  p <- cv[[nm]]; r <- roc(y, p, quiet=TRUE, ci=TRUE)
  th <- coords(r, "best", ret=c("threshold","sensitivity","specificity"))
  cat(sprintf("  %-4s AUC=%.3f (95%%CI %.3f-%.3f) | Youden sens=%.3f spec=%.3f\n",
              nm, auc(r), ci.auc(r)[1], ci.auc(r)[3], th$sensitivity[1], th$specificity[1]))
}

# ---- 3. leave-one-compound-out (13 compounds; vehicles in train) ----
compounds <- sort(unique(trt[y == 1]))
cat("\nLEAVE-ONE-COMPOUND-OUT (train excl. compound c, test on c), logistic:\n")
res <- sapply(compounds, function(c) {
  te <- !is.na(trt) & trt == c
  if (sum(te) < 3) return(NA)
  d <- as.data.frame(X)
  m <- glm(y[!te] ~ ., data=d[!te,,drop=FALSE], family=binomial)
  tryCatch(auc1(predict(m, d[te,,drop=FALSE], type="response"), y[te]), error=function(e) NA)
})
for (c in compounds) cat(sprintf("  %-28s n=%-3d AUC=%.3f\n", c, sum(!is.na(trt) & trt==c), res[c]))
cat("  mean(LOCO):", round(mean(res, na.rm=TRUE), 3), "\n")

cat("\nDONE\n"); sink()
cat("written to outputs/results_mouse_cv.txt\n")
