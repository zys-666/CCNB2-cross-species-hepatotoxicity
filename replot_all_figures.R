# replot_all_figures.R : regenerate residual/boxplot (original DALEX multi style) + importance (ggplot)
suppressMessages({library(DALEX); library(ggplot2); library(dplyr)})
setwd("D:/MLroute2")
o <- "outputs/v2"
rds <- readRDS(file.path(o, "models_for_figures.rds"))
mp <- lapply(rds$explainers, model_performance)

# residual.pdf —— 原版显式调用
pdf(file.path(o, "residual.pdf"), width = 6, height = 6)
p1 <- plot(mp[["RF"]], mp[["SVM"]], mp[["GLM"]], mp[["GBM"]], mp[["KNN"]],
           mp[["NNET"]], mp[["LASSO"]], mp[["DT"]], mp[["XGB"]])
print(p1)
dev.off()
# boxplot.pdf —— 原版显式调用 + geom=boxplot
pdf(file.path(o, "boxplot.pdf"), width = 6, height = 6)
p2 <- plot(mp[["RF"]], mp[["SVM"]], mp[["GLM"]], mp[["GBM"]], mp[["KNN"]],
           mp[["NNET"]], mp[["LASSO"]], mp[["DT"]], mp[["XGB"]], geom = "boxplot")
print(p2)
dev.off()

# importance.pdf —— 聚合每模型 dropout loss，45 基因按跨模型均值排序，9 模型折线叠加
imp_agg <- lapply(names(rds$vi_raw), function(nm) {
  vi <- rds$vi_raw[[nm]]
  vi <- vi[vi$variable != "_baseline_", ]
  d <- aggregate(dropout_loss ~ variable, data = vi, FUN = mean)
  d$model <- nm
  d
})
imp_df <- do.call(rbind, imp_agg)
imp_df$gene <- imp_df$variable
ord <- imp_df %>% group_by(variable) %>% summarise(m = mean(dropout_loss)) %>%
  arrange(desc(m)) %>% pull(variable)
imp_df$variable <- factor(imp_df$variable, levels = rev(ord))   # 最重要的在最上
cols9 <- c("RF"="#e6194b","SVM"="#3cb44b","GLM"="#ffe119","GBM"="#4363d8","KNN"="#f58231",
           "NNET"="#911eb4","LASSO"="#46f0f0","DT"="#f032e6","XGB"="#bcf60c")
pdf(file.path(o, "importance.pdf"), width = 9, height = 11)
p3 <- ggplot(imp_df, aes(x = dropout_loss, y = variable, color = model, group = model)) +
  geom_line(aes(group = model), alpha = 0.45) +
  geom_point(size = 1.6) +
  scale_color_manual(values = cols9) +
  labs(x = "mean dropout loss (permutation importance)", y = "",
       title = "Variable importance across nine classifiers (training side)") +
  theme_bw(base_size = 11)
print(p3)
dev.off()

sizes <- file.info(file.path(o, c("residual.pdf", "boxplot.pdf", "importance.pdf")))$size
print(setNames(round(sizes), c("residual.pdf", "boxplot.pdf", "importance.pdf")))
cat(ifelse(all(sizes > 10000), "OK\n", "WARN\n"))
