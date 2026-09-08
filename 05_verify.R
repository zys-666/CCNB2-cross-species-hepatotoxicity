# 05_verify.R : alignment check + positive-control genes + re-run CV with forced alignment
suppressMessages({library(data.table); library(pROC); library(randomForest); library(xgboost)})
setwd("D:/MLroute2")
sink("outputs/results_verify.txt", split = TRUE)

r45 <- readRDS("outputs/rat_expr45.rds")
X45 <- r45$X; gsm <- r45$gsm; m <- r45$m
# ---- 1. alignment check ----
cat("identical(expr column GSM order, m$gsm order):", identical(gsm, m$gsm), "\n")
y1 <- ifelse(m$L1 == "Treat", 1, 0)
y2 <- ifelse(m$L1[match(gsm, m$gsm)] == "Treat", 1, 0)
cat("y1 vs y2 (match-forced) identical:", identical(y1, y2), "\n")

# ---- 2. positive control: classic toxicant-response genes from full matrix ----
expr <- fread("inputs/normalize.txt", sep="\t", header=TRUE, check.names=FALSE)
gn <- toupper(expr[[1]])
ctrl_genes <- c("Cyp1a1","Mt1","Gstm1","Hmox1","Ahr")
found <- intersect(ctrl_genes, gn)
cat("positive-control genes found:", found, "\n")
auc1 <- function(p, yy) as.numeric(auc(roc(yy, p, quiet=TRUE)))
Xc <- as.matrix(expr[match(found, gn), -1, with=FALSE])
for (i in seq_len(nrow(Xc))) {
  d <- Xc[i, ] - 0
  cat(sprintf("  %-6s mean_dosed-mean_ctrl = %+.3f | single-gene AUC = %.3f\n",
      found[i], mean(d[y1==1]) - mean(d[y1==0]), auc1(d, y1)))
}

# ---- 3. re-run 45-gene CV with forced-alignment y2 (safety) ----
X45t <- t(X45); n <- ncol(X45)
Pg <- Pr <- Px <- rep(NA, n)
for (r in 1:2) {
  set.seed(200 + r); f <- sample(rep(1:5, length.out = n))
  for (i in 1:5) {
    te <- f == i
    m1 <- glm.fit(X45t[!te, , drop=FALSE], y2[!te], family = binomial())
    Pg[te] <- plogis(as.numeric(X45t[te, , drop=FALSE] %*% m1$coefficients))
    m2 <- randomForest(x = X45t[!te, , drop=FALSE], y = factor(y2[!te]), ntree = 300)
    Pr[te] <- predict(m2, X45t[te, , drop=FALSE], type="prob")[, "1"]
    set.seed(600 + r*10 + i)
    m3 <- xgb.train(params = list(objective="binary:logistic", eval_metric="auc", max_depth=3,
                                  eta=0.05, subsample=0.8, colsample_bytree=0.8),
                    data = xgb.DMatrix(X45t[!te, , drop=FALSE], label = y2[!te]), nrounds = 100, verbose = 0)
    Px[te] <- predict(m3, X45t[te, , drop=FALSE])
  }
}
for (nm in c("glm","rf","xgb")) {
  P <- get(paste0("P", substr(nm,1,1))); cc <- auc1(P, y2)
  cat(sprintf("forced-align 45-gene 5x2 CV: %-4s AUC=%.3f\n", nm, cc))
}
cat("\nDONE\n"); sink(); cat("written\n")
