############################################################
## 机器学习训练集代码 v2（单基因 CCNB2 · 共识筛选版）
## 由原版《机器学习训练集代码.txt》升级而来。改动清单：
## 【改1】样本标签：不再解析 normalize.txt 列名的 _Control/_Treat 后缀
##        （GEO 剂量复核：数值 dose==0 → 279 对照 / dose>0 → 1939 处理，与原文件标签一致）；
##        改为读取 GEO 样本元数据 rat_samples.tsv，
##        规则 L1：数值 dose == 0 → Control（279），否则 Treat（1939，含低剂量臂 0.01–0.65 mg/kg）；
##        并按 GSM 与表达矩阵列显式对齐（merge 会改变行序，必须 match 对齐）。
## 【改2】删除"70/30 随机划分 + 30% 测试集 AUC"的验证式表述。
##        模型在修正标签的 FULL 数据上训练；性能=分层 5 折 x 3 次重复内部 CV，
##        仅作描述性内部判别报告（训练物种内部，非跨物种验证）。
## 【改3】DALEX 重要性 & SHAP 的计算对象由"30% 测试集"改为"训练侧样本"
##        （无信息泄漏：不再用测试样本参与任何基因排序）。
## 【改4】新增输出：CCNB2 在全部模型/指标中的排名表 rank_CCNB2.txt。
## 【改5】补全随机种子管理 & sessionInfo() 版本记录；并行调参（doParallel）。
## 【保留】9 个模型结构（8 caret + 原生 XGB）、SHAP 蜂群图全套、重要性文件命名。
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
## 1. 文件路径（改为本项目工作目录；输入文件随补充材料提供）
############################################################
setwd("D:/MLroute2")
inputFile  <- "inputs/normalize.txt"                 # GSE57815 表达矩阵（探针注释后 log2 值）
geneFile   <- "interGenes_input.txt"                 # = interGenes.txt（53 候选基因）
metaFile   <- "inputs/rat_samples.tsv"               # 【改1】GEO 元数据(acc,dose,vehicle,time,route)
outDir     <- "outputs/v2"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

############################################################
## 2. 读取表达矩阵（fread 加速；等价于原 read.table）
############################################################
expr_all <- fread(inputFile, sep = "\t", header = TRUE, check.names = FALSE)
geneIDcol <- names(expr_all)[1]
expr <- as.matrix(expr_all[, -1, with = FALSE])
rownames(expr) <- expr_all[[geneIDcol]]
rm(expr_all); gc()

############################################################
## 3. 读取基因列表并取交集（沿用原版大小写无关匹配）
############################################################
geneList_raw <- as.vector(read.table(geneFile, header = FALSE, sep = "\t",
                                     stringsAsFactors = FALSE)[, 1])
commonGenes <- expr_genes_up <- NULL
expr_genes_raw <- rownames(expr); expr_genes_up <- toupper(expr_genes_raw)
commonGenes_up <- intersect(unique(toupper(geneList_raw)), expr_genes_up)
commonGenes <- expr_genes_raw[expr_genes_up %in% commonGenes_up]
expr <- expr[commonGenes, ]
cat("Input genes:", length(geneList_raw), " Used genes:", length(commonGenes), "\n")

############################################################
## 4. 转置 + 分组信息（【改1】L1 标签，按 GSM 对齐元数据）
############################################################
gsm_cols <- colnames(expr)                            # 形如 GSM1392188_Control
gsm <- sub("_(Control|Treat)$", "", gsm_cols)
stopifnot(all(grepl("^GSM", gsm)))
meta <- fread(metaFile)
mi <- match(gsm, meta$acc); stopifnot(!anyNA(mi))
meta <- meta[mi]; stopifnot(identical(gsm, meta$acc))
dnum <- suppressWarnings(as.numeric(sub(" .*", "", trimws(meta$dose))))
dose0 <- !is.na(dnum) & dnum == 0           # 数值剂量==0 → 对照（279）；不能用 grepl("^0")
Type <- ifelse(dose0, "Control", "Treat")   # 修正后 L1 = 279 Control / 1939 Treat（与原文件标签一致）
stopifnot(length(Type) == ncol(expr))

data <- as.data.frame(t(expr))
data$Type <- factor(Type, levels = c("Control", "Treat"))
cat("L1 labels: Control", sum(Type == "Control"), " Treat", sum(Type == "Treat"), "\n")
labelAudit <- data.frame(GSM = gsm, Type = Type, dose0 = dose0, vehicle = meta$vehicle)
write.table(labelAudit, file.path(outDir, "labels_L1.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

############################################################
## 5. 零方差 & 缺失值处理（在 FULL 训练侧做，fit 与 transform 同一数据）
############################################################
nzv <- nearZeroVar(data[, -ncol(data)])
if (length(nzv) > 0) data <- data[, -nzv]
preProc <- preProcess(data[, -ncol(data)], method = "medianImpute")
data <- predict(preProc, data)
nzv <- nearZeroVar(data[, -ncol(data)])
if (length(nzv) > 0) data <- data[, -nzv]

############################################################
## 6. caret 交叉验证设置（【改2】5 折 x 3 重复内部 CV；保留 savePredictions 供 OOF ROC）
############################################################
ctrl <- trainControl(method = "repeatedcv", number = 5, repeats = 3,
                     classProbs = TRUE, savePredictions = TRUE,
                     summaryFunction = twoClassSummary,
                     allowParallel = TRUE)

nCores <- min(6, parallel::detectCores() - 2)
cl <- makePSOCKcluster(nCores); registerDoParallel(cl)
cat("parallel cores:", nCores, "\n")

############################################################
## 7. 训练 8 个 caret 模型（结构同原版；metric=ROC）
############################################################
set.seed(123)
mod_rf    <- train(Type ~ ., data, method = "rf",    metric = "ROC", trControl = ctrl)
mod_svm   <- train(Type ~ ., data, method = "svmRadial", metric = "ROC", trControl = ctrl)
mod_glm   <- train(Type ~ ., data, method = "glm", family = "binomial", metric = "ROC", trControl = ctrl)
mod_gbm   <- train(Type ~ ., data, method = "gbm", metric = "ROC", trControl = ctrl, verbose = FALSE)
mod_knn   <- train(Type ~ ., data, method = "knn", metric = "ROC", trControl = ctrl)
mod_nnet  <- train(Type ~ ., data, method = "nnet", metric = "ROC", trControl = ctrl, trace = FALSE)
mod_lasso <- train(Type ~ ., data, method = "glmnet", metric = "ROC", trControl = ctrl)
mod_dt    <- train(Type ~ ., data, method = "rpart", metric = "ROC", trControl = ctrl)

############################################################
## 8. 原生 XGBoost（超参沿用原版：max_depth=3, eta=0.05, subsample=0.8,
##    colsample_bytree=0.8, nrounds=100；另做同结构 5x3 手动 CV 取 OOF 预测）
############################################################
x_mat <- as.matrix(data[, setdiff(colnames(data), "Type")])
y_bin <- ifelse(data$Type == "Treat", 1, 0)
params_xgb <- list(objective = "binary:logistic", eval_metric = "auc",
                   max_depth = 3, eta = 0.05, subsample = 0.8, colsample_bytree = 0.8)
set.seed(123)
mod_xgb <- xgb.train(params = params_xgb,
                     data = xgb.DMatrix(x_mat, label = y_bin), nrounds = 100, verbose = 0)
# XGB 内部 CV：手动分层 5 折 x 3 重复（保证 OOF 预测与行序对齐）
oof_aucs <- numeric(3); oof_xgb_plot <- NULL
n <- nrow(x_mat)
for (r in 1:3) {
  set.seed(200 + r)
  f <- sample(rep(1:5, length.out = n))
  p <- rep(NA, n)
  for (k in 1:5) {
    te <- f == k
    b <- xgb.train(params = params_xgb,
                   data = xgb.DMatrix(x_mat[!te, , drop = FALSE], label = y_bin[!te]),
                   nrounds = 100, verbose = 0)
    p[te] <- predict(b, x_mat[te, , drop = FALSE])
  }
  oof_aucs[r] <- as.numeric(auc(roc(y_bin, p, quiet = TRUE)))
  if (r == 1) oof_xgb_plot <- p
}

############################################################
## 9. 内部 CV 性能汇总（【改2】描述性；caret 取 15 次重采样的 ROC/Sens/Spec）
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
## 10. 训练侧样本：解释器 & 重要性 & SHAP 的统一数据（【改3】不再是测试集）
############################################################
set.seed(123)
trainSampleN <- min(500, nrow(data))
idxS <- sample(seq_len(nrow(data)), trainSampleN)
xSample_df <- data[idxS, setdiff(colnames(data), "Type")]
ySample    <- y_bin[idxS]

prob_caret <- function(m, x) as.numeric(predict(m, x, type = "prob")[, "Treat"])
prob_fun <- list(
  RF    = prob_caret, SVM = prob_caret, GLM = prob_caret, GBM = prob_caret,
  KNN   = prob_caret, NNET = prob_caret, LASSO = prob_caret, DT = prob_caret,
  XGB   = function(m, x) as.numeric(predict(m, as.matrix(x))))
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
## 10b. 模型性能残差图（原版绘图风格，逐字还原；数据=训练侧解释器）
############################################################
mp_list <- lapply(explainers, model_performance)
pdf(file.path(outDir, "residual.pdf"), width = 6, height = 6)
tryCatch({
  p1 <- plot(mp_list[["RF"]], mp_list[["SVM"]], mp_list[["GLM"]], mp_list[["GBM"]],
             mp_list[["KNN"]], mp_list[["NNET"]], mp_list[["LASSO"]], mp_list[["DT"]],
             mp_list[["XGB"]])
  print(p1)
}, error = function(e) cat("residual plot err:", e$message, "\n"))
dev.off()
pdf(file.path(outDir, "boxplot.pdf"), width = 6, height = 6)
tryCatch({
  p2 <- plot(mp_list[["RF"]], mp_list[["SVM"]], mp_list[["GLM"]], mp_list[["GBM"]],
             mp_list[["KNN"]], mp_list[["NNET"]], mp_list[["LASSO"]], mp_list[["DT"]],
             mp_list[["XGB"]], geom = "boxplot")
  print(p2)
}, error = function(e) cat("boxplot err:", e$message, "\n"))
dev.off()

############################################################
## 11. 变量重要性（DALEX 置换，训练侧）→ importanceGene.*.txt
############################################################
vi_raw <- lapply(names(explainers), function(nm) {
  vi <- variable_importance(explainers[[nm]], loss_function = loss_root_mean_square)
  vi[order(vi$dropout_loss), ]                       # baseline 损失最小，自然在首行；保留原始行类供绘图
})
names(vi_raw) <- names(explainers)
imp_list <- lapply(vi_raw, function(vi) {
  vi <- vi[vi$variable != "_baseline_", ]
  # DALEX 默认每变量多次置换：按变量聚合为平均 dropout loss
  vi <- aggregate(dropout_loss ~ variable, data = vi, FUN = mean)
  vi <- vi[order(-vi$dropout_loss), ]
  colnames(vi) <- c("variable", "mean_dropout_loss")
  vi
})
names(imp_list) <- names(explainers)

############################################################
## 11b. 重要性图（原版绘图风格：baseline + 每模型 top10 基因，9 模型叠加）
############################################################
pdf(file.path(outDir, "importance.pdf"), width = 9, height = 11)
tryCatch({
  imp_agg <- lapply(names(vi_raw), function(nm) {
    v <- vi_raw[[nm]]
    v <- v[v$variable != "_baseline_", ]
    d <- aggregate(dropout_loss ~ variable, data = v, FUN = mean)
    d$model <- nm
    d
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
## 12. 【改4】CCNB2 排名表（每模型置换重要性 + 单基因判别 AUC）
############################################################
rawAUC <- sapply(commonGenes, function(g) {
  v <- data[[g]]
  as.numeric(auc(roc(data$Type, v, quiet = TRUE, levels = c("Control", "Treat"), direction = "<")))
})
rankTab <- data.frame(gene = commonGenes, AUC_treatVsCtrl = rawAUC,
                      discAUC = pmax(rawAUC, 1 - rawAUC))
for (nm in names(imp_list)) {
  vi <- imp_list[[nm]]
  rk <- rank(-vi$mean_dropout_loss)
  names(rk) <- vi$variable
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
## 13. ROC 图（原版绘图风格逐字还原；数据=内部 CV 合并 OOF 预测，图例含 AUC+95%CI）
############################################################
auc_ci_label <- function(name, rocobj) {
  ci <- ci.auc(rocobj)
  paste0(name, " (AUC:", sprintf("%.3f", auc(rocobj)), ", 95%CI:",
         sprintf("%.3f", ci[1]), "-", sprintf("%.3f", ci[3]), ")")
}
roc_list <- lapply(names(mods), function(nm) {
  p <- mods[[nm]]$pred
  roc(p$obs, p$Treat, quiet = TRUE, levels = c("Control", "Treat"), direction = "<")
})
roc_xgb <- roc(y_bin, oof_xgb_plot, quiet = TRUE, direction = "<")
roc_list <- c(roc_list, list(roc_xgb))
names(roc_list) <- c(names(mods), "XGB")

pdf(file.path(outDir, "ROC.pdf"), width = 6, height = 6)
cols9 <- c("red", "blue", "green", "orange", "purple", "brown", "black", "pink", "gray")
plot(roc_list[[1]], legacy.axes = TRUE, col = cols9[1], lwd = 2)
for (i in 2:9) plot(roc_list[[i]], col = cols9[i], lwd = 2, add = TRUE)
legend("bottomright",
       legend = sapply(names(roc_list), function(nm) auc_ci_label(nm, roc_list[[nm]])),
       col = cols9, lwd = 2, bty = "n", cex = 0.8)
dev.off()

############################################################
## 14. SHAP 蜂群图全套（【改3】在训练侧 40 样本上，B=20，沿用原版画法）
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
             oof_xgb_plot = oof_xgb_plot, y_bin = y_bin, data = data, commonGenes = commonGenes),
        file.path(outDir, "models_for_figures.rds"))
cat("saved models_for_figures.rds (re-plot without retraining)\n")
stopCluster(cl)

############################################################
## 15. 【改5】版本记录
############################################################
sink(file.path(outDir, "sessionInfo.txt"))
cat("R:", R.version.string, "\n")
print(sessionInfo())
sink()
cat("\nALL DONE. outputs in", outDir, "\n")
