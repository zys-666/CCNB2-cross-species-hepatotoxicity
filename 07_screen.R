# 07_screen.R : CCNB2 rank among 45 interGenes (rat L1 & mouse) + concordance context
suppressMessages({library(data.table); library(pROC); library(randomForest)})
setwd("D:/MLroute2")
sink("outputs/results_screen.txt", split = TRUE)
r45 <- readRDS("outputs/rat_expr45.rds")
X <- r45$X; m <- r45$m
y <- ifelse(m$L1 == "Treat", 1, 0)
mg <- fread("inputs/mouse_expr.csv", sep=",", header=TRUE, check.names=FALSE)
M <- as.matrix(mg[, -1, with=FALSE]); rownames(M) <- mg[[1]]
grp <- fread("inputs/mouse_group.csv"); setnames(grp, c("SampleID","group"))
ym <- ifelse(grp$group[match(colnames(M), grp$SampleID)] == "Treat", 1, 0)
genes <- rownames(X); shared <- intersect(toupper(genes), toupper(rownames(M)))
Xm <- M[toupper(rownames(M)) %in% shared, , drop=FALSE]
Xm <- Xm[match(toupper(genes)[toupper(genes) %in% shared], toupper(rownames(Xm))), , drop=FALSE]
auc1 <- function(v, yy) as.numeric(auc(roc(yy, v, quiet=TRUE)))
wt <- function(a, b) suppressWarnings(t.test(a, b)$p.value)

dA <- sapply(seq_len(nrow(X)), function(i) mean(X[i, y==1]) - mean(X[i, y==0]))
pA <- sapply(seq_len(nrow(X)), function(i) wt(X[i, y==1], X[i, y==0]))
aucA <- sapply(seq_len(nrow(X)), function(i) auc1(X[i, ], y))
dM <- sapply(seq_len(nrow(Xm)), function(i) mean(Xm[i, ym==1]) - mean(Xm[i, ym==0]))
pM <- sapply(seq_len(nrow(Xm)), function(i) wt(Xm[i, ym==1], Xm[i, ym==0]))
aucM <- sapply(seq_len(nrow(Xm)), function(i) auc1(Xm[i, ], ym))

tab <- data.frame(gene = genes, dRat = round(dA, 3), pRat = signif(pA, 2), aucRat = round(aucA, 3),
                  dMouse = round(dM, 3), pMouse = signif(pM, 2), aucMouse = round(aucM, 3))
tab$rank_aucRat <- rank(-tab$aucRat); tab$rank_aucMouse <- rank(-tab$aucMouse)
tab$concord_up <- tab$dRat > 0 & tab$dMouse > 0 & tab$pRat < 0.05 & tab$pMouse < 0.05
tab$concord_dir <- sign(tab$dRat) == sign(tab$dMouse)
setorder(setDT(tab), rank_aucRat)
print(tab, nrows = 50)
cat("\nCCNB2 (Ccnb2) row:\n")
print(tab[toupper(gene) == "CCNB2"])
cat("\nAmong 45 genes:\n")
cat("  concordant UP both species (p<0.05):", sum(tab$concord_up), "/45\n")
cat("  concordant direction (any):", sum(tab$concord_dir), "/45\n")
cat("  genes UP in rat (p<0.05):", sum(tab$dRat > 0 & tab$pRat < 0.05),
    "| UP in mouse (p<0.05):", sum(tab$dMouse > 0 & tab$pMouse < 0.05), "\n")
cat("  rat AUC >= CCNB2's:", sum(tab$aucRat >= tab$aucRat[toupper(gene)=="CCNB2"]), "/45\n")
cat("  mouse AUC >= CCNB2's:", sum(tab$aucMouse >= tab$aucMouse[toupper(gene)=="CCNB2"]), "/45\n")
# RF importance rank (rat, L1)
set.seed(1); rf <- randomForest(x = t(X), y = factor(y), ntree = 500, importance = TRUE)
imp <- importance(rf, type = 1)
tab$rfImpRank <- rank(-imp[, 1])
cat("  CCNB2 RF(MeanDecreaseAccuracy) rank:", tab$rfImpRank[toupper(gene) == "CCNB2"], "/45\n")
cat("\nDONE\n"); sink()
