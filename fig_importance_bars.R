# fig_importance_bars.R : single-canvas horizontal bars, 9 model colors, per-model own top10 genes
suppressMessages({library(ggplot2); library(dplyr)})
setwd("D:/MLroute2")
rds <- readRDS("outputs/v2/models_for_figures.rds")

my_colors <- c(RF="#E41A1C", SVM="#377EB8", XGB="#4DAF4A", GLM="#984EA3", GBM="#FF7F00",
               KNN="#E6AB02", NNET="#A65628", LASSO="#F781BF", DT="#666666")
model_order <- c("RF","SVM","XGB","GLM","GBM","KNN","NNET","LASSO","DT")

agg <- lapply(names(rds$vi_raw), function(nm) {
  v <- rds$vi_raw[[nm]]
  v <- v[v$variable != "_baseline_", ]
  d <- aggregate(dropout_loss ~ variable, data = v, FUN = mean)
  d$model <- nm
  d <- d[order(d$dropout_loss), ]
  d <- tail(d, 10)                       # 每模型自己的 top10（loss 升序尾10）
  d
})
df <- do.call(rbind, agg)
df$ycat <- paste0(df$model, "::", df$variable)
# 视觉顺序 top→bottom：每个模型块内 loss 最大的基因在最上；模型块按 model_order
top_down <- unlist(lapply(model_order, function(m) {
  d <- df[df$model == m, ]
  d <- d[order(-d$dropout_loss), ]       # 块内：大 loss 在顶
  d$ycat
}))
df$ycat <- factor(df$ycat, levels = rev(top_down))   # ggplot y 第1个level在底部 → 反转
df$lbl <- sub(".*::", "", as.character(df$ycat))

p <- ggplot(df, aes(x = dropout_loss, y = ycat, fill = model)) +
  geom_col(width = 0.72, alpha = 0.92) +
  scale_fill_manual(values = my_colors, breaks = model_order) +
  scale_y_discrete(labels = setNames(df$lbl, df$ycat)) +
  labs(title = "Feature Importance",
       subtitle = "created for the DT, GBM, GLM, KNN, LASSO, NNET, RF, SVM, XGB model",
       x = "Root mean square error (RMSE) loss after permutations", y = NULL) +
  theme_bw(base_size = 10) +
  theme(axis.text.y = element_text(size = 6.5),
        legend.position = "right", legend.title = element_blank())

nrow_all <- nrow(df)
pdf("outputs/v2/importance.pdf", width = 9, height = max(12, nrow_all * 0.135))
print(p)
dev.off()
cat("rows:", nrow_all, "| size:", file.info("outputs/v2/importance.pdf")$size, "\n")
