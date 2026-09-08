# 03b_loco.R : corrected leave-one-compound-out for mouse (control half held out per fold)
suppressMessages({library(data.table); library(pROC)})
setwd("D:/MLroute2")
sink("outputs/results_loco.txt", split = TRUE)
grp <- fread("inputs/mouse_group.csv"); setnames(grp, c("SampleID","group"))
st <- fread("inputs/mouse_samples.tsv")
mg <- fread("inputs/mouse_expr.csv", sep=",", header=TRUE, check.names=FALSE)
mgenes <- mg[[1]]
M <- as.matrix(mg[, -1, with=FALSE]); rownames(M) <- mgenes
g3 <- c("Ccnb2","Plk1","Ube2c")
M3 <- M[toupper(rownames(M)) %in% toupper(g3), , drop=FALSE]
y <- ifelse(grp$group[match(colnames(M), grp$SampleID)] == "Treat", 1, 0)
trt <- st$treatment[match(colnames(M), st$acc)]
X <- t(M3)
auc1 <- function(p, yy) { if (length(unique(yy)) < 2) return(NA); as.numeric(auc(roc(yy, p, quiet=TRUE))) }
compounds <- sort(unique(trt[y == 1]))
ctrl_ids <- which(y == 0)
cat("LOCO design: test = compound-c treated + held-out half of controls; train = other treated + other controls; logistic, 3 genes\n")
mat <- sapply(compounds, function(c) {
  te_treat <- which(!is.na(trt) & trt == c & y == 1)
  a <- replicate(5, {
    set.seed(1000 + match(c, compounds) * 10 + sample(1:100, 1)); # random control split
    set.seed(sample(1:1e5, 1))
    h <- sample(ctrl_ids, floor(length(ctrl_ids) / 2))
    tr <- setdiff(seq_along(y), c(te_treat, h))
    m <- glm(y[tr] ~ ., data = as.data.frame(X[tr, , drop=FALSE]), family = binomial)
    auc1(predict(m, as.data.frame(X[c(te_treat, h), , drop=FALSE]), type="response"), y[c(te_treat, h)])
  })
  mean(a, na.rm = TRUE)
})
for (c in compounds) cat(sprintf("  %-28s n_trt=%-3d LOCO AUC=%.3f\n", c, sum(trt==c & y==1), mat[c]))
cat("  mean LOCO (logistic, 3 genes):", round(mean(mat, na.rm=TRUE), 3), "\n")
sink(); cat("written\n")
