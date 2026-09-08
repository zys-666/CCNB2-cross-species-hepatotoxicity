# 08_nine.R : 9-classifier consensus screening on rat L1 corrected labels
suppressMessages({library(data.table); library(pROC); library(randomForest); library(xgboost)
  library(gbm); library(kernlab); library(glmnet); library(nnet); library(rpart)})
setwd("D:/MLroute2")
sink("outputs/results_nine.txt", split = TRUE)
r45 <- readRDS("outputs/rat_expr45.rds")
X <- r45$X; m <- r45$m
y <- ifelse(m$L1 == "Treat", 1, 0)
genes <- rownames(X); Xt <- t(X)              # n x p
cat("rat L1:", sum(y==1), "vs", sum(y==0), "| p =", ncol(Xt), "genes\n")
set.seed(1)
auc1 <- function(p, yy) as.numeric(auc(roc(yy, p, quiet=TRUE)))

# ---------- fit 9 models once on full corrected data ----------
cat("fitting...\n")
m_rf  <- randomForest(x=Xt, y=factor(y), ntree=500, importance=TRUE)
set.seed(1); m_xgb <- xgb.train(params=list(objective="binary:logistic", max_depth=3, eta=0.05,
    subsample=0.8, colsample_bytree=0.8), data=xgb.DMatrix(Xt, label=y), nrounds=100, verbose=0)
set.seed(1); m_gbm <- gbm.fit(Xt, y, distribution="bernoulli", n.trees=300, interaction.depth=3,
    shrinkage=0.05, verbose=FALSE)
m_glm <- glm.fit(Xt, y, family=binomial())
set.seed(1); m_svm <- ksvm(Xt, factor(y), kernel="rbfdot", prob.model=TRUE)
m_las <- cv.glmnet(Xt, y, family="binomial", alpha=1)
m_knn <- knn_train <- NULL   # knn has no global importance -> permutation only
set.seed(1); m_nnet <- nnet(Xt, y, size=3, maxit=300, decay=0.01, trace=FALSE)
m_dt  <- rpart(factor(y) ~ ., data=as.data.frame(Xt), method="class")

predf <- list(
  rf  = function(d) predict(m_rf, d, type="prob")[, "1"],
  xgb = function(d) predict(m_xgb, d),
  gbm = function(d) predict(m_gbm, d, n.trees=300, type="response"),
  glm = function(d) plogis(as.numeric(d %*% m_glm$coefficients)),
  svm = function(d) predict(m_svm, d, type="probabilities")[, "1"],
  las = function(d) as.numeric(predict(m_las, d, s="lambda.min", type="response")),
  nnet= function(d) as.numeric(predict(m_nnet, d)),
  dt  = function(d) predict(m_dt, as.data.frame(d), type="prob")[, "1"])
# knn: kNN on full data is meaningless as a global model -> report via univariate only; skip in consensus

# ---------- permutation importance (uniform metric, train sample) ----------
set.seed(7); S <- sample(nrow(Xt), min(250, nrow(Xt)))
imp <- sapply(names(predf), function(nm) {
  p0 <- auc1(predf[[nm]](Xt[S, , drop=FALSE]), y[S])
  sapply(seq_len(ncol(Xt)), function(j) {
    Xp <- Xt[S, , drop=FALSE]; Xp[, j] <- sample(Xp[, j])
    p0 - auc1(predf[[nm]](Xp), y[S])
  })
})
rownames(imp) <- genes

# native importances
nat <- data.frame(gene=genes,
  rf_native  = importance(m_rf, type=1)[, 1],
  xgb_native = setNames(xgb.importance(model=m_xgb)$Gain, xgb.importance(model=m_xgb)$Feature)[genes],
  gbm_native = summary(m_gbm, plotit=FALSE)$rel.inf[match(genes, summary(m_gbm, plotit=FALSE)$var)],
  las_coef   = abs(as.numeric(coef(m_las, s="lambda.min")))[-1][match(genes, colnames(Xt))],
  dt_native  = m_dt$variable.importance[genes])
nat[is.na(nat)] <- 0
# ranks per metric
rk <- data.frame(gene=genes)
for (nm in colnames(imp)) rk[[paste0("perm_", nm)]] <- rank(-imp[, nm])
for (nm in c("rf_native","xgb_native","gbm_native","las_coef","dt_native"))
  rk[[nm]] <- rank(-nat[[nm]])
rk$uniAUC <- rank(-sapply(seq_len(nrow(X)), function(i) auc1(X[i, ], y)))
# mouse univariate rank for context
mg <- fread("inputs/mouse_expr.csv", sep=",", header=TRUE, check.names=FALSE)
M <- as.matrix(mg[, -1, with=FALSE]); rownames(M) <- mg[[1]]
grp <- fread("inputs/mouse_group.csv"); setnames(grp, c("SampleID","group"))
ym <- ifelse(grp$group[match(colnames(M), grp$SampleID)] == "Treat", 1, 0)
Ms <- M[toupper(rownames(M)) %in% toupper(genes), , drop=FALSE]; Ms <- Ms[match(toupper(genes), toupper(rownames(Ms))), , drop=FALSE]
rk$mouse_uniAUC <- rank(-sapply(seq_len(nrow(Ms)), function(i) auc1(Ms[i, ], ym)))

cc <- rk[toupper(rk$gene) == "CCNB2", ]
cat("\nCCNB2 rank (1 = best) across metrics/models:\n")
print(cc)
rk$meanRank <- rowMeans(rk[, -1])
cat("\nCCNB2 mean rank across all metrics:", rk$meanRank[toupper(rk$gene) == "CCNB2"], "/45\n")
cat("\ntop5 by mean rank:\n"); print(head(rk[order(rk$meanRank), c("gene","meanRank")], 5))
cat("\ntop5 by mean rank (excluding univariate):\n")
rk2 <- rk; rk2$meanRank2 <- rowMeans(rk2[, !names(rk2) %in% c("gene","uniAUC","mouse_uniAUC")])
print(head(rk2[order(rk2$meanRank2), c("gene","meanRank2")], 5))
fwrite(rk, "outputs/rank_all_metrics.tsv", sep="\t")
cat("\nDONE\n"); sink()
