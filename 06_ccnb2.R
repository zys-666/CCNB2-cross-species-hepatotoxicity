# 06_ccnb2.R : single-gene CCNB2 cross-cohort consistency (rat L1 / mouse)
suppressMessages({library(data.table); library(pROC)})
setwd("D:/MLroute2")
sink("outputs/results_ccnb2.txt", split = TRUE)

r45 <- readRDS("outputs/rat_expr45.rds")
X45 <- r45$X; m <- r45$m
g <- rownames(X45)[toupper(rownames(X45)) == "CCNB2"]
x <- X45[g, ]
y <- ifelse(m$L1 == "Treat", 1, 0)

stat1 <- function(v, yv, tag) {
  a <- v[yv == 1]; b <- v[yv == 0]
  tt <- t.test(a, b)
  r <- tryCatch(roc(yv, v, quiet = TRUE, ci = TRUE), error = function(e) NULL)
  cat(sprintf("%-22s nT=%4d nC=%4d | Treat=%.3f Ctrl=%.3f | d=%.3f (95%%CI %.3f-%.3f) p=%.2e",
      tag, sum(yv == 1), sum(yv == 0), mean(a), mean(b), tt$estimate[1] - tt$estimate[2],
      tt$conf.int[1], tt$conf.int[2], tt$p.value))
  if (!is.null(r)) cat(sprintf(" | AUC=%.3f (%.3f-%.3f)", auc(r), ci.auc(r)[1], ci.auc(r)[3]))
  cat("\n")
}
cat("========== CCNB2 single-gene, log2 expression ==========\n")
cat("-- rat GSE57815 --\n")
stat1(x, y, "rat L1 (all dose0=Ctrl)")
m3 <- m[!is.na(m$L3), ]; stat1(x[!is.na(m$L3)], y[!is.na(m$L3)], "rat L3 (3 vehicles)")
# by timepoint: dosed vs same-timepoint vehicle controls
cat("-- rat L1 by timepoint (dosed vs same-time vehicle controls) --\n")
for (tm in sort(unique(m$time))) {
  ii <- m$time == tm & !is.na(m$time)
  if (sum(ii & y == 1) > 5 && sum(ii & y == 0) > 5) stat1(x[ii], y[ii], paste("rat time", tm))
}
# by vehicle: dosed vs vehicle controls within same vehicle
cat("-- rat L1 by vehicle (dosed vs same-vehicle controls) --\n")
for (vh in c("Corn Oil","CMC","Water","Saline")) {
  ii <- m$vehicle == vh
  if (sum(ii & y == 1) > 5 && sum(ii & y == 0) > 5) stat1(x[ii], y[ii], paste("rat veh", vh))
}

cat("\n-- mouse GSE44783 --\n")
mg <- fread("inputs/mouse_expr.csv", sep=",", header=TRUE, check.names=FALSE)
mm <- as.matrix(mg[, -1, with=FALSE]); rownames(mm) <- mg[[1]]
grp <- fread("inputs/mouse_group.csv"); setnames(grp, c("SampleID","group"))
ym <- ifelse(grp$group[match(colnames(mm), grp$SampleID)] == "Treat", 1, 0)
xm <- mm[toupper(rownames(mm)) == "CCNB2", ]
stat1(xm, ym, "mouse all")
st <- fread("inputs/mouse_samples.tsv")
trt <- st$treatment[match(colnames(mm), st$acc)]
tm_ <- st$time[match(colnames(mm), st$acc)]
for (tmv in sort(unique(tm_))) { ii <- tm_ == tmv; stat1(xm[ii], ym[ii], paste("mouse time", tmv)) }
cat("-- mouse per compound (each vs all 120 controls) --\n")
for (c in sort(unique(trt[ym == 1]))) {
  ii <- !is.na(trt) & trt == c & ym == 1
  jj <- ym == 0
  if (sum(ii) >= 3) {
    a <- xm[ii]; b <- xm[jj]; tt <- t.test(a, b)
    r <- tryCatch(auc(roc(c(rep(1,length(a)),rep(0,length(b))), c(a,b), quiet=TRUE)), error=function(e) NA)
    cat(sprintf("  %-28s n=%2d d=%+.3f p=%.2e AUC=%.3f\n", c, sum(ii),
        mean(a)-mean(b), tt$p.value, r))
  }
}
cat("\nDONE\n"); sink(); cat("written\n")
