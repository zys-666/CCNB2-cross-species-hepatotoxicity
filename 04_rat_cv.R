# 04_rat_cv.R : honest repeated CV on rat under corrected labels (L1)
suppressMessages({library(pROC); library(randomForest); library(xgboost); library(glmnet)})
setwd("D:/MLroute2")
sink("outputs/results_rat_cv.txt", split = TRUE)
r45 <- readRDS("outputs/rat_expr45.rds")
X45 <- r45$X; m <- r45$m
lab <- m$L1; y <- ifelse(lab == "Treat", 1, 0)
g3 <- rownames(X45)[toupper(rownames(X45)) %in% toupper(c("Ccnb2","Plk1","Ube2c"))]
cat("rat L1: Treat", sum(y==1), "Control", sum(y==0), "| 45 genes; 3 genes:", g3, "\n")
X3 <- X45[g3, , drop=FALSE]
cat("gene means (log2) all samples: Ccnb2", round(mean(X3[1,]),3),
    "Plk1", round(mean(X3[2,]),3), "Ube2c", round(mean(X3[3,]),3), "\n")
auc1 <- function(p, yy) { if (length(unique(p)) < 2) return(NA); r <- roc(yy, p, quiet=TRUE, ci=TRUE)
  c(auc=as.numeric(auc(r)), lo=as.numeric(ci.auc(r)[1]), hi=as.numeric(ci.auc(r)[3])) }
run_cv <- function(Xm, tag) {
  cat("\n---", tag, "---\n")
  n <- ncol(Xm); Xt <- t(Xm)
  Pg <- Pr <- Px <- rep(NA, n)
  for (r in 1:3) {
    set.seed(100 + r); f <- sample(rep(1:5, length.out = n))
    for (i in 1:5) {
      te <- f == i
      m1 <- glm.fit(Xt[!te, , drop=FALSE], y[!te], family = binomial())
      Pg[te] <- plogis(as.numeric(Xt[te, , drop=FALSE] %*% m1$coefficients))
      m2 <- randomForest(x = Xt[!te, , drop=FALSE], y = factor(y[!te]), ntree = 300)
      Pr[te] <- predict(m2, Xt[te, , drop=FALSE], type="prob")[, "1"]
      set.seed(500 + r * 10 + i)
      m3 <- xgb.train(params = list(objective="binary:logistic", eval_metric="auc", max_depth=3,
                                    eta=0.05, subsample=0.8, colsample_bytree=0.8),
                      data = xgb.DMatrix(Xt[!te, , drop=FALSE], label = y[!te]), nrounds = 100, verbose = 0)
      Px[te] <- predict(m3, Xt[te, , drop=FALSE])
    }
  }
  for (nm in c("glm","rf","xgb")) {
    P <- get(paste0("P", substr(nm,1,1)))
    cc <- auc1(P, y)
    cat(sprintf("  %-4s CV-AUC=%.3f (95%%CI %.3f-%.3f)\n", nm, cc["auc"], cc["lo"], cc["hi"]))
  }
}
run_cv(X3, "rat 3 genes (L1) 5x3 CV")
run_cv(X45, "rat 45 genes (L1) 5x3 CV")
cat("\nDONE\n"); sink(); cat("written to outputs/results_rat_cv.txt\n")
