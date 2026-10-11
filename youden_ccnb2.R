# youden_ccnb2.R : single-gene CCNB2 Youden-optimal sens/spec (rat L1 & mouse)
suppressMessages({library(data.table); library(pROC)})
setwd("D:/MLroute2")
r45 <- readRDS("outputs/rat_expr45.rds")
m <- r45$m; y <- ifelse(m$L1 == "Treat", 1, 0)
g <- rownames(r45$X)[toupper(rownames(r45$X)) == "CCNB2"]
xr <- r45$X[g, ]

mg <- fread("inputs/mouse_expr.csv", sep = ",", header = TRUE, check.names = FALSE)
M <- as.matrix(mg[, -1, with = FALSE]); rownames(M) <- mg[[1]]
grp <- fread("inputs/mouse_group.csv"); setnames(grp, c("SampleID", "group"))
ym <- ifelse(grp$group[match(colnames(M), grp$SampleID)] == "Treat", 1, 0)
xm <- M[toupper(rownames(M)) == "CCNB2", ]

youden <- function(v, yy, tag) {
  r <- roc(yy, v, quiet = TRUE, direction = "<", ci = TRUE)
  th <- coords(r, "best", ret = c("threshold", "sensitivity", "specificity"))
  cat(sprintf("%-24s AUC=%.3f (%.3f-%.3f) | Youden thr=%.3f | Sens=%.3f Spec=%.3f | nT=%d nC=%d\n",
      tag, auc(r), ci.auc(r)[1], ci.auc(r)[3], th$threshold[1],
      th$sensitivity[1], th$specificity[1], sum(yy == 1), sum(yy == 0)))
  data.frame(Cohort = tag, AUC = round(as.numeric(auc(r)), 3),
             CI_lo = round(ci.auc(r)[1], 3), CI_hi = round(ci.auc(r)[3], 3),
             Youden_thr = round(th$threshold[1], 3),
             Sens = round(th$sensitivity[1], 3), Spec = round(th$specificity[1], 3))
}
out <- rbind(
  youden(xr, y, "Rat GSE57815 (CCNB2)"),
  youden(xm, ym, "Mouse GSE44783 (CCNB2)"))
write.table(out, "outputs/paper/CCNB2_youden.txt", sep = "\t", quote = FALSE, row.names = FALSE)
cat("written outputs/paper/CCNB2_youden.txt\n")
