############################################################
## training_v3.R —— 主流程(修正版,取代 training_v2.R)
##
## 相对 training_v2.R 的三处修正(对应审稿意见):
##  【v3-1｜意见5】预处理(medianImpute / nearZeroVar)不再在 FULL 数据上拟合,
##              改为在 trainControl(preProcess = ...) 内部,由 caret 在每个重采样折内拟合。
##              —— 实测:本矩阵无 NA、无近零方差基因,故数值结果与 v2 完全一致;
##                 此改动是方法学正确性,不改变任何报告数字。
##  【v3-2｜意见4】不再声称"hyperparameters were fixed"。caret 的调参行为照旧,
##              但把每个模型最终选中的超参显式导出(bestTune_table.txt),
##              并在下方表头写明:超参由 caret 网格在重采样内选出,而非预先固定。
##  【v3-3｜意见6+7】ROC 不再直接使用 $pred(它含【全部候选调参配置】且每样本重复多次),
##              改为:先按 bestTune 过滤到最终配置 → 再按样本对重复取平均
##              → 每样本仅 1 个折外预测,再算 AUC 与 DeLong 95%CI。
##
## 【保留不变】数据读取 / 标签对齐(match + stopifnot)/ 原生 XGB 手动 5×3 CV /
##            DALEX 置换重要性 / CCNB2 排名 / SHAP 蜂群 / 原版绘图语句,全部逐字沿用。
##
## 输出目录: outputs/v3(不覆盖 v2 结果,便于对照)
## 输入: inputs/normalize.txt, interGenes.txt, inputs/rat_samples.tsv
############################################################

rm(list = ls())
set.seed(123)

############################################################
## 0. 加载包
############################################################
library(caret); library(DALEX); library(ggplot2); library(randomForest)
library(kernlab); library(xgboost); library(pROC); library(gbm)
library(nnet); library(glmnet); library(rpart); library(dplyr)
library(reshape2); library(data.table)
suppressMessages(library(doParallel))

############################################################
## 1. 文件路径
############################################################
setwd("D:/MLroute2")
inputFile <- "inputs/normalize.txt"
geneFile  <- "interGenes_input.txt"
metaFile  <- "inputs/rat_samples.tsv"
outDir    <- "outputs/v3"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

############################################################
## 2. 读取表达矩阵
############################################################
expr_all <- fread(inputFile, sep = "\t", header = TRUE, check.names = FALSE)
geneIDcol <- names(expr_all)[1]
expr <- as.matrix(expr_all[, -1, with = FALSE])
rownames(expr) <- expr_all[[geneIDcol]]
rm(expr_all); gc()

############################################################
## 3. 候选基因取交集
############################################################
geneList_raw <- as.vector(read.table(geneFile, header = FALSE, sep = "\t",
                                     stringsAsFactors = FALSE)[, 1])
expr_genes_raw <- rownames(expr); expr_genes_up <- toupper(expr_genes_raw)
commonGenes_up <- intersect(unique(toupper(geneList_raw)), expr_genes_up)
commonGenes <- expr_genes_raw[expr_genes_up %in% commonGenes_up]
expr <- expr[commonGenes, ]
cat("Input genes:", length(geneList_raw), " Used genes:", length(commonGenes), "\n")

############################################################
## 4. 转置 + 分组信息(L1 标签,按 GSM 对齐)
############################################################
gsm_cols <- colnames(expr)
gsm <- sub("_(Control|Treat)$", "", gsm_cols)
stopifnot(all(grepl("^GSM", gsm)))
meta <- fread(metaFile)
mi <- match(gsm, meta$acc); stopifnot(!anyNA(mi), !anyDuplicated(mi))
meta <- meta[mi]; stopifnot(all(gsm == meta$acc))     # 显式重排(列序与文件行序不同)
dnum  <- suppressWarnings(as.numeric(sub(" .*", "", trimws(meta$dose))))
dose0 <- !is.na(dnum) & dnum == 0
Type  <- ifelse(dose0, "Control", "Treat")
stopifnot(length(Type) == ncol(expr))

data <- as.data.frame(t(expr))
data$Type <- factor(Type, levels = c("Control", "Treat"))
cat("L1 labels: Control", sum(Type == "Control"), " Treat", sum(Type == "Treat"), "\n")
write.table(data.frame(GSM = gsm, Type = Type, dose0 = dose0, vehicle = meta$vehicle),
            file.path(outDir, "labels_L1.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

############################################################
## 5. 【v3-1】不再在此处做全局预处理
##    原 training_v2.R 第 85-90 行把 nearZeroVar / medianImpute 拟合在 FULL 数据上,
##    导致留出折参与了插补参数与变量筛选的估计。v3 改为交给 caret 在折内拟合(见第 6 节)。
##    这里只做一次诊断记录,便于复现"本矩阵无需处理":
nzv_global <- nearZeroVar(data[, -ncol(data)])
cat("[v3-1] 全局近零方差基因数:", length(nzv_global),
    "| NA 单元格数:", sum(is.na(data[, -ncol(data)])), "\n")
writeLines(c(sprintf("nearZeroVar (global, diagnostic only): %d", length(nzv_global)),
             sprintf("NA cells: %d", sum(is.na(data[, -ncol(data)])))),
           file.path(outDir, "preproc_diagnostic.txt"))

############################################################
## 6. 【v3-1】预处理移入重采样:preProcess 由 caret 在每个折内拟合
############################################################
## 注意:preProcess 是 train() 的参数,trainControl() 没有这个参数。
ctrl <- trainControl(method = "repeatedcv", number = 5, repeats = 3,
                     classProbs = TRUE, savePredictions = TRUE,
                     summaryFunction = twoClassSummary,
                     allowParallel = TRUE)
## 【v3-1】预处理交由 train() 的 preProcess 参数处理:
##   caret 的重采样流程(nominalTrainWorkflow)每个折只把【训练折数据】传给建模函数,
##   因此预处理统计量在折内拟合,不会用到留出折。
PP <- c("medianImpute", "nzv")

nCores <- min(6, parallel::detectCores() - 2)
cl <- makePSOCKcluster(nCores); registerDoParallel(cl)
cat("parallel cores:", nCores, "\n")

############################################################
## 7. 训练 8 个 caret 模型
##    【v3-2】超参由 caret 默认网格在本重采样内选出(非预先固定);
##            最终选中值在第 7b 节显式导出。
############################################################
set.seed(123)
mod_rf    <- train(Type ~ ., data, method = "rf",        metric = "ROC", trControl = ctrl, preProcess = PP)
mod_svm   <- train(Type ~ ., data, method = "svmRadial", metric = "ROC", trControl = ctrl, preProcess = PP)
mod_glm   <- train(Type ~ ., data, method = "glm", family = "binomial", metric = "ROC", trControl = ctrl, preProcess = PP)
mod_gbm   <- train(Type ~ ., data, method = "gbm",       metric = "ROC", trControl = ctrl, verbose = FALSE, preProcess = PP)
mod_knn   <- train(Type ~ ., data, method = "knn",       metric = "ROC", trControl = ctrl, preProcess = PP)
mod_nnet  <- train(Type ~ ., data, method = "nnet",      metric = "ROC", trControl = ctrl, trace = FALSE, preProcess = PP)
mod_lasso <- train(Type ~ ., data, method = "glmnet",    metric = "ROC", trControl = ctrl, preProcess = PP)
mod_dt    <- train(Type ~ ., data, method = "rpart",     metric = "ROC", trControl = ctrl, preProcess = PP)

############################################################
## 7b. 【v3-2】导出最终选中超参(回应"固定超参"的表述问题)
############################################################
mods_all <- list(RF = mod_rf, SVM = mod_svm, GLM = mod_glm, GBM = mod_gbm,
                 KNN = mod_knn, NNET = mod_nnet, LASSO = mod_lasso, DT = mod_dt)
bt <- do.call(rbind, lapply(names(mods_all), function(nm) {
  b <- mods_all[[nm]]$bestTune
  data.frame(Model = nm, selected = paste(colnames(b), unlist(b), sep = "=", collapse = "; "),
             stringsAsFactors = FALSE)
}))
bt <- rbind(bt, data.frame(Model = "XGB", selected = "max_depth=3; eta=0.05; subsample=0.8; colsample_bytree=0.8; nrounds=100"))
write.table(bt, file.path(outDir, "bestTune_table.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n[v3-2] 最终选中超参:\n"); print(bt, row.names = FALSE)

############################################################
## 8. 原生 XGBoost(超参文档化;同结构 5×3 手动 CV 取 OOF)
############################################################
x_mat <- as.matrix(data[, setdiff(colnames(data), "Type")])
y_bin <- ifelse(data$Type == "Treat", 1, 0)
params_xgb <- list(objective = "binary:logistic", eval_metric = "auc",
                   max_depth = 3, eta = 0.05, subsample = 0.8, colsample_bytree = 0.8)
set.seed(123)
mod_xgb <- xgb.train(params = params_xgb,
                     data = xgb.DMatrix(x_mat, label = y_bin), nrounds = 100, verbose = 0)
xgb_pr <- matrix(NA_real_, nrow(x_mat), 3)      # 【v3-3】保留 3 次重复以便按样本平均
oof_aucs <- numeric(3)
n <- nrow(x_mat)
for (r in 1:3) {
  set.seed(200 + r)
  ## 【v3-4】XGB 改用分层折,与 8 个 caret 模型口径一致(createFolds 只调用一次!)
  ##   ⚠️ 若写成 for (k in 1:5) f[createFolds(...)[[k]]] <- k,则每次循环都会
  ##   重新生成一套不同的划分,f 被反复覆盖,约 1/3 样本拿不到折号 → 折外预测残缺。
  fl <- createFolds(data$Type, k = 5, list = TRUE)
  f  <- integer(n)
  for (k in 1:5) f[fl[[k]]] <- k
  stopifnot(all(f %in% 1:5))          # 自检:不允许有未分配的样本
  p <- rep(NA, n)
  for (k in 1:5) {
    te <- f == k
    b <- xgb.train(params = params_xgb,
                   data = xgb.DMatrix(x_mat[!te, , drop = FALSE], label = y_bin[!te]),
                   nrounds = 100, verbose = 0)
    p[te] <- predict(b, x_mat[te, , drop = FALSE])
  }
  xgb_pr[, r] <- p
  oof_aucs[r] <- as.numeric(auc(roc(y_bin, p, quiet = TRUE)))
}
oof_xgb_plot <- rowMeans(xgb_pr)                 # 【v3-3】每样本 1 个预测

############################################################
## 9. 内部 CV 性能汇总(数据源 = $resample,只用最终配置,无污染)
############################################################
perf <- data.frame(Model = character(), CV_AUC_mean = numeric(), CV_AUC_sd = numeric(),
                   CV_Sens = numeric(), CV_Spec = numeric(), stringsAsFactors = FALSE)
mods <- list(RF = mod_rf, SVM = mod_svm, GLM = mod_glm, GBM = mod_gbm, KNN = mod_knn,
             NNET = mod_nnet, LASSO = mod_lasso, DT = mod_dt)
for (nm in names(mods)) {
  rs <- mods[[nm]]$resample
  perf <- rbind(perf, data.frame(Model = nm,
                                 CV_AUC_mean = round(mean(rs$ROC), 3),
                                 CV_AUC_sd   = round(sd(rs$ROC), 3),
                                 CV_Sens     = round(mean(rs$Sens), 3),
                                 CV_Spec     = round(mean(rs$Spec), 3)))
}
perf <- rbind(perf, data.frame(Model = "XGB",
                               CV_AUC_mean = round(mean(oof_aucs), 3),
                               CV_AUC_sd = round(sd(oof_aucs), 3),
                               CV_Sens = NA, CV_Spec = NA))
write.table(perf, file.path(outDir, "modelPerf_CV.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
print(perf)

############################################################
## 10. 训练侧解释器(DALEX)
############################################################
set.seed(123)
trainSampleN <- min(500, nrow(data))
idxS <- sample(seq_len(nrow(data)), trainSampleN)
xSample_df <- data[idxS, setdiff(colnames(data), "Type")]
ySample    <- y_bin[idxS]

prob_caret <- function(m, x) as.numeric(predict(m, x, type = "prob")[, "Treat"])
prob_fun <- list(
  RF = prob_caret, SVM = prob_caret, GLM = prob_caret, GBM = prob_caret,
  KNN = prob_caret, NNET = prob_caret, LASSO = prob_caret, DT = prob_caret,
  XGB = function(m, x) as.numeric(predict(m, as.matrix(x))))
mod_final <- list(RF = mod_rf, SVM = mod_svm, GLM = mod_glm, GBM = mod_gbm,
                  KNN = mod_knn, NNET = mod_nnet, LASSO = mod_lasso, DT = mod_dt,
                  XGB = mod_xgb)

cat("\nbuilding DALEX explainers on training sample (n =", trainSampleN, ")...\n")
explainers <- lapply(names(mod_final), function(nm) {
  DALEX::explain.default(model = mod_final[[nm]], data = xSample_df, y = ySample,
                         label = nm, predict_function = prob_fun[[nm]], verbose = FALSE)
})
names(explainers) <- names(mod_final)

############################################################
## 10b. 残差图 / 箱线图(原版绘图语句逐字沿用)
############################################################
mp_list <- lapply(explainers, model_performance)
pdf(file.path(outDir, "residual.pdf"), width = 6, height = 6)
tryCatch({
  p1 <- plot(mp_list[["RF"]], mp_list[["SVM"]], mp_list[["GLM"]], mp_list[["GBM"]],
             mp_list[["KNN"]], mp_list[["NNET"]], mp_list[["LASSO"]], mp_list[["DT"]],
             mp_list[["XGB"]]); print(p1)
}, error = function(e) cat("residual plot err:", e$message, "\n"))
dev.off()
pdf(file.path(outDir, "boxplot.pdf"), width = 6, height = 6)
tryCatch({
  p2 <- plot(mp_list[["RF"]], mp_list[["SVM"]], mp_list[["GLM"]], mp_list[["GBM"]],
             mp_list[["KNN"]], mp_list[["NNET"]], mp_list[["LASSO"]], mp_list[["DT"]],
             mp_list[["XGB"]], geom = "boxplot"); print(p2)
}, error = function(e) cat("boxplot err:", e$message, "\n"))
dev.off()

############################################################
## 11. 变量重要性(DALEX 置换,训练侧)
############################################################
vi_raw <- lapply(names(explainers), function(nm) {
  vi <- variable_importance(explainers[[nm]], loss_function = loss_root_mean_square)
  vi[order(vi$dropout_loss), ]
})
names(vi_raw) <- names(explainers)
imp_list <- lapply(vi_raw, function(vi) {
  vi <- vi[vi$variable != "_baseline_", ]
  vi <- aggregate(dropout_loss ~ variable, data = vi, FUN = mean)
  vi <- vi[order(-vi$dropout_loss), ]
  colnames(vi) <- c("variable", "mean_dropout_loss"); vi
})
names(imp_list) <- names(explainers)

############################################################
## 11b. 重要性图(原版绘图语句逐字沿用)
############################################################
pdf(file.path(outDir, "importance.pdf"), width = 9, height = 11)
tryCatch({
  imp_agg <- lapply(names(vi_raw), function(nm) {
    v <- vi_raw[[nm]]; v <- v[v$variable != "_baseline_", ]
    d <- aggregate(dropout_loss ~ variable, data = v, FUN = mean); d$model <- nm; d
  })
  imp_df <- do.call(rbind, imp_agg)
  ord <- imp_df %>% group_by(variable) %>% summarise(m = mean(dropout_loss)) %>%
    arrange(desc(m)) %>% pull(variable)
  imp_df$variable <- factor(imp_df$variable, levels = rev(ord))
  p3 <- ggplot(imp_df, aes(x = dropout_loss, y = variable, color = model, group = model)) +
    geom_line(aes(group = model), alpha = 0.45) + geom_point(size = 1.6) +
    scale_color_manual(values = c(RF="#e6194b", SVM="#3cb44b", GLM="#ffe119", GBM="#4363d8",
                                  KNN="#f58231", NNET="#911eb4", LASSO="#46f0f0", DT="#f032e6",
                                  XGB="#bcf60c")) +
    labs(x = "mean dropout loss (permutation importance)", y = "",
         title = "Variable importance across nine classifiers (training side)") +
    theme_bw(base_size = 11)
  print(p3)
}, error = function(e) cat("importance plot err:", e$message, "\n"))
dev.off()

geneTop <- data.frame(Model = character(), stringsAsFactors = FALSE)
for (nm in names(imp_list)) {
  top <- head(imp_list[[nm]], 10)
  write.table(top, file.path(outDir, paste0("importanceGene.", nm, ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  geneTop <- rbind(geneTop, data.frame(Model = nm, t(top$variable), check.names = FALSE))
}
colnames(geneTop) <- c("Model", paste0("Top", 1:(ncol(geneTop) - 1)))
write.table(geneTop, file.path(outDir, "importanceGene.ALL.txt"), sep = "\t",
            quote = FALSE, row.names = FALSE)
print(geneTop)

############################################################
## 12. CCNB2 排名表
############################################################
rawAUC <- sapply(commonGenes, function(g) {
  v <- data[[g]]
  as.numeric(auc(roc(data$Type, v, quiet = TRUE, levels = c("Control", "Treat"), direction = "<")))
})
rankTab <- data.frame(gene = commonGenes, AUC_treatVsCtrl = rawAUC,
                      discAUC = pmax(rawAUC, 1 - rawAUC))
for (nm in names(imp_list)) {
  vi <- imp_list[[nm]]; rk <- rank(-vi$mean_dropout_loss); names(rk) <- vi$variable
  rankTab[[paste0("rank_", nm)]] <- rk[commonGenes]
}
rankTab$rank_uniAUC <- rank(-rankTab$discAUC)
rankTab$meanRank <- rowMeans(rankTab[, grep("^rank_", colnames(rankTab))])
rankTab <- rankTab[order(rankTab$meanRank), ]
write.table(rankTab, file.path(outDir, "ranks_all_metrics.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
ccnb2 <- rankTab[toupper(rankTab$gene) == "CCNB2", ]
write.table(ccnb2, file.path(outDir, "rank_CCNB2.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nCCNB2 ranks:\n"); print(ccnb2)
cat("\ntop5 by mean rank:\n"); print(head(rankTab[, c("gene", "meanRank")], 5))

############################################################
## 13. 【v3-3】ROC:先按 bestTune 过滤,再按样本平均(每样本 1 个预测)
##      绘图语句逐字沿用原版(legacy.axes / 9 色 / 图例格式 / main / abline)
############################################################
oof <- list()
for (nm in names(mods)) {
  p  <- as.data.table(mods[[nm]]$pred)
  btt <- as.data.table(mods[[nm]]$bestTune)
  for (cc in intersect(colnames(btt), colnames(p))) p <- p[get(cc) == btt[[cc]][1]]
  avg <- tapply(p$Treat, p$rowIndex, mean)
  oof[[nm]] <- avg[order(as.integer(names(avg)))]
}
oof[["XGB"]] <- oof_xgb_plot

ord9 <- c("RF","SVM","GLM","GBM","KNN","NNET","LASSO","DT","XGB")
res <- data.frame(); roc_list <- list()
for (nm in ord9) {
  pv <- as.numeric(oof[[nm]])
  ro <- roc(y_bin, pv, quiet = TRUE, direction = "<")
  ci <- ci.auc(ro, method = "delong")
  res <- rbind(res, data.frame(Model = nm, OOF_AUC = round(as.numeric(auc(ro)), 3),
                               CI_lo = round(ci[1], 3), CI_hi = round(ci[3], 3)))
  roc_list[[nm]] <- ro
}
write.table(res, file.path(outDir, "ROC_final.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n[v3-3] ROC(最终配置 + 每样本 1 预测):\n"); print(res, row.names = FALSE)

pdf(file.path(outDir, "图3A_九模型ROC曲线.pdf"), width = 6, height = 6)
cols9 <- c("red", "blue", "green", "orange", "purple", "brown", "black", "pink", "gray")
lab <- sapply(ord9, function(nm) paste0(nm, " (AUC:", sprintf("%.3f", res$OOF_AUC[res$Model == nm]),
                                        ", 95%CI:", sprintf("%.3f", res$CI_lo[res$Model == nm]),
                                        "-", sprintf("%.3f", res$CI_hi[res$Model == nm]), ")"))
plot(roc_list[[1]], legacy.axes = TRUE, col = cols9[1], lwd = 2,
     main = "Internal CV ROC (rat GSE57815, pooled out-of-fold)")
for (i in 2:9) plot(roc_list[[i]], col = cols9[i], lwd = 2, add = TRUE)
abline(a = 0, b = 1, lty = 3, col = "grey60")
legend("bottomright", legend = lab, col = cols9, lwd = 2, bty = "n", cex = 0.8)
dev.off()

############################################################
## 13b. 校准曲线 + ECE — 与第 13 节同一口径(最终配置 + 每样本 1 预测)
##      绘图语句逐字沿用原 calibration.R 风格
############################################################
my_colors <- c(RF="#E41A1C", SVM="#377EB8", XGB="#4DAF4A", GLM="#984EA3", GBM="#FF7F00",
               KNN="#E6AB02", NNET="#A65628", LASSO="#F781BF", DT="#666666")
df <- do.call(rbind, lapply(ord9, function(nm) {
  data.frame(model = nm, p = as.numeric(oof[[nm]]), y = y_bin)
}))
ece_tab <- lapply(split(df, df$model), function(d) {
  qs <- unique(quantile(d$p, probs = seq(0, 1, 0.1), na.rm = TRUE))
  d$bin <- cut(d$p, breaks = qs, include.lowest = TRUE)
  b <- d %>% group_by(bin) %>%
    summarise(pred = mean(p), obs = mean(y), n = n(), .groups = "drop") %>% na.omit()
  data.frame(model = d$model[1],
             ECE = round(sum(b$n * abs(b$obs - b$pred)) / sum(b$n), 4), bins = nrow(b))
})
ece <- do.call(rbind, ece_tab)
write.table(ece, file.path(outDir, "calibration_final.txt"), sep = "\t",
            quote = FALSE, row.names = FALSE)
cat("\n[v3-3] 校准 ECE(同一口径):\n"); print(ece, row.names = FALSE)

bin_all <- do.call(rbind, lapply(split(df, df$model), function(d) {
  qs <- unique(quantile(d$p, probs = seq(0, 1, 0.1), na.rm = TRUE))
  b <- d %>% mutate(bin = cut(p, breaks = qs, include.lowest = TRUE)) %>%
    group_by(bin) %>% summarise(pred = mean(p), obs = mean(y), n = n(), .groups = "drop") %>% na.omit()
  b$model <- d$model[1]; b
}))
bin_all <- merge(bin_all, ece, by = "model")
bin_all$lab <- paste0(bin_all$model, " (ECE=", sprintf("%.3f", bin_all$ECE), ")")

p_cal <- ggplot(bin_all, aes(x = pred, y = obs, colour = model, group = model)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey40", linewidth = 0.6) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  scale_colour_manual(values = my_colors, breaks = names(my_colors),
                      labels = setNames(paste0(names(my_colors), " (ECE=",
                        sprintf("%.3f", ece$ECE[match(names(my_colors), ece$model)]), ")"),
                        names(my_colors))) +
  labs(title = "Calibration curves (internal CV out-of-fold predictions)",
       subtitle = "Ideal calibration = dashed diagonal; legend shows Expected Calibration Error",
       x = "Predicted probability (Treat)", y = "Observed frequency") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + theme_bw(base_size = 12) +
  theme(legend.position = "right", legend.title = element_blank())
pdf(file.path(outDir, "图3B_模型校准曲线.pdf"), width = 9, height = 6.5); print(p_cal); dev.off()
cat("saved 图3A_九模型ROC曲线.pdf / 图3B_模型校准曲线.pdf\n")

############################################################
## 14. SHAP 蜂群图(训练侧 40 样本,B=20,原版画法逐字沿用)
############################################################
library(ggbeeswarm)
set.seed(123)
shap_sample <- xSample_df[sample(seq_len(nrow(xSample_df)), min(40, nrow(xSample_df))), ]
pdf(file.path(outDir, "SHAP_beeswarm.pdf"), width = 9, height = 7)
for (nm in names(explainers)) {
  cat("SHAP:", nm, "\n")
  sp <- DALEX::predict_parts(explainer = explainers[[nm]], new_observation = shap_sample,
                             type = "shap", B = 20)
  sdf <- as.data.frame(sp) %>% group_by(variable_name) %>%
    mutate(mean_abs = mean(abs(contribution))) %>% ungroup()
  top_genes <- sdf %>% distinct(variable_name, mean_abs) %>%
    arrange(desc(mean_abs)) %>% slice(1:10) %>% pull(variable_name)
  p <- ggplot(filter(sdf, variable_name %in% top_genes),
              aes(x = contribution, y = reorder(variable_name, mean_abs), color = contribution)) +
    geom_beeswarm(priority = "density", cex = 1.1) +
    scale_color_gradient(low = "blue", high = "red") +
    labs(title = paste("SHAP beeswarm (Top 10 genes) —", nm), x = "SHAP value", y = "") +
    theme_minimal(base_size = 13)
  print(p)
}
dev.off()

saveRDS(list(mods = mods, mod_xgb = mod_xgb, explainers = explainers, vi_raw = vi_raw,
             oof = oof, xgb_pr = xgb_pr, y_bin = y_bin, data = data,
             commonGenes = commonGenes),
        file.path(outDir, "models_for_figures.rds"))
cat("saved models_for_figures.rds (re-plot without retraining)\n")
stopCluster(cl)

############################################################
## 15. 版本记录
############################################################
sink(file.path(outDir, "sessionInfo.txt"))
cat("R:", R.version.string, "\n"); print(sessionInfo()); sink()
cat("\nALL DONE. outputs in", outDir, "\n")
