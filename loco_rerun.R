# loco_rerun.R : leave-one-compound-out, mouse GSE44783, control-half holdout design
# 3 gene sets: single CCNB2 | 3-gene panel | 45-gene panel ; logistic, 5 control-split repeats per compound
suppressMessages({library(data.table); library(pROC)})
setwd("D:/MLroute2")
grp <- fread("inputs/mouse_group.csv"); setnames(grp, c("SampleID", "group"))
st  <- fread("inputs/mouse_samples.tsv")
mg  <- fread("inputs/mouse_expr.csv", sep = ",", header = TRUE, check.names = FALSE)
M   <- as.matrix(mg[, -1, with = FALSE]); rownames(M) <- mg[[1]]
y   <- ifelse(grp$group[match(colnames(M), grp$SampleID)] == "Treat", 1, 0)
trt <- st$treatment[match(colnames(M), st$acc)]
compounds <- sort(unique(trt[y == 1]))
ctrl_ids <- which(y == 0)

g45 <- NULL   # 小鼠阵列无 rat-45 全映射；LOCO 用单基因 + 3基因两套（与可用基因集一致）
gene_sets <- list(
  CCNB2 = c("Ccnb2"),
  g3    = c("Ccnb2", "Plk1", "Ube2c"))
cat("LOCO gene sets: single CCNB2 + 3-gene panel (mouse array symbols)\n")

auc1 <- function(p, yy) { if (length(unique(yy)) < 2) return(NA); as.numeric(auc(roc(yy, p, quiet = TRUE))) }

out_tab <- lapply(names(gene_sets), function(gs) {
  X <- t(M[toupper(rownames(M)) %in% toupper(gene_sets[[gs]]), , drop = FALSE])
  per <- sapply(compounds, function(c) {
    te_treat <- which(!is.na(trt) & trt == c & y == 1)
    a <- replicate(5, {
      set.seed(5000 + match(c, compounds) * 7 + sample(1:97, 1))
      h <- sample(ctrl_ids, floor(length(ctrl_ids) / 2))
      tr <- setdiff(seq_along(y), c(te_treat, h))
      m <- glm(y[tr] ~ ., data = as.data.frame(X[tr, , drop = FALSE]), family = binomial)
      auc1(predict(m, as.data.frame(X[c(te_treat, h), , drop = FALSE]), type = "response"),
           y[c(te_treat, h)])
    })
    mean(a, na.rm = TRUE)
  })
  data.frame(compound = compounds, n_trt = as.integer(table(trt[y == 1])[compounds]),
             AUC = round(as.numeric(per), 3), stringsAsFactors = FALSE)
})
names(out_tab) <- names(gene_sets)

cat("\n=========== LOCO results ===========\n")
for (gs in names(out_tab)) {
  cat(sprintf("\n-- gene set: %s --\n", gs))
  print(out_tab[[gs]], row.names = FALSE)
  cat(sprintf("mean LOCO AUC (%s): %.3f\n", gs, round(mean(out_tab[[gs]]$AUC, na.rm = TRUE), 3)))
}
comb <- Reduce(function(a, b) merge(a, b, by = c("compound", "n_trt"), sort = FALSE),
               out_tab)
names(comb)[3:4] <- paste0("LOCO_AUC_", names(out_tab))
dir.create("outputs/paper", showWarnings = FALSE)
write.table(comb, "outputs/paper/GSE44783_LOCO.txt", sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nwritten outputs/paper/GSE44783_LOCO.txt\n")
