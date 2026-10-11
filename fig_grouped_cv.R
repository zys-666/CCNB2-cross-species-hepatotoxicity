############################################################
## fig_grouped_cv.R — 新增补充图:compound-grouped CV 的证据图
##   图S-1  三口径对比(5×3 CV / 嵌套 CV / 化合物分组 CV),点-误差棒
##   图S-2  化合物分组 CV 的逐折 AUC(散点 + 均值线),展示折间稳定性
## 绘图风格逐字沿用原 fig_style2.R(theme_classic 变体、点-误差棒、配色体系)
############################################################
suppressMessages({library(data.table); library(ggplot2)})
setwd("D:/MLroute2")
out <- "outputs/paper/figs_grouped"; dir.create(out, showWarnings = FALSE, recursive = TRUE)

th <- theme_classic(base_size = 10) +
  theme(axis.line = element_line(linewidth = 0.4, colour = "grey20"),
        axis.ticks = element_line(linewidth = 0.4, colour = "grey20"),
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5, colour = "grey30"),
        legend.key.size = unit(0.35, "cm"))

cv  <- fread("outputs/v2/modelPerf_CV.txt")
ncv <- fread("outputs/v2/nestedCV_perf.txt")
grp <- fread("outputs/v2/grouped/rat_grouped_perf.txt")
pf  <- fread("outputs/v2/grouped/rat_grouped_perf_perFold.txt")

## ---------------- 图S-1:三口径对比 ----------------
m <- merge(cv[, .(Model, CV_AUC_mean, CV_AUC_sd)],
           ncv[, .(Model, NestedAUC_mean, NestedAUC_sd)], by = "Model")
m <- merge(m, grp[, .(Model = model, groupedCV_AUC_mean, groupedCV_AUC_sd)], by = "Model")
m <- m[order(-groupedCV_AUC_mean)]
lvl <- rev(m$Model)

d <- rbind(
  data.table(Model = m$Model, metric = "5\u00d73 CV",              auc = m$CV_AUC_mean,        sd = m$CV_AUC_sd),
  data.table(Model = m$Model, metric = "Nested CV",                auc = m$NestedAUC_mean,     sd = m$NestedAUC_sd),
  data.table(Model = m$Model, metric = "Compound-grouped CV",      auc = m$groupedCV_AUC_mean, sd = m$groupedCV_AUC_sd))
d[, Model  := factor(Model, levels = lvl)]
d[, metric := factor(metric, levels = c("5\u00d73 CV", "Nested CV", "Compound-grouped CV"))]

p1 <- ggplot(d, aes(x = auc, y = Model, colour = metric, shape = metric)) +
  geom_errorbarh(aes(xmin = auc - sd, xmax = auc + sd), height = 0.16, linewidth = 0.45,
                 position = position_dodge(width = 0.62)) +
  geom_point(size = 2.1, position = position_dodge(width = 0.62)) +
  scale_colour_manual(values = c("5\u00d73 CV" = "#3B6EA5",
                                 "Nested CV" = "#C0504D",
                                 "Compound-grouped CV" = "#2E7D5B")) +
  scale_shape_manual(values = c("5\u00d73 CV" = 16, "Nested CV" = 17, "Compound-grouped CV" = 15)) +
  scale_x_continuous(limits = c(0.90, 1.005), breaks = seq(0.90, 1.00, 0.02)) +
  labs(x = "AUC (mean \u00b1 SD across resampling folds)", y = NULL, colour = NULL, shape = NULL,
       title = "Compound-grouped validation of the nine rat screening classifiers",
       subtitle = "Rat GSE57815; 258 groups = 200 compound + 58 control-batch") +
  th + theme(legend.position = c(0.02, 0.06), legend.justification = c(0, 0),
             legend.background = element_rect(fill = alpha("white", 0.75), colour = NA))
pdf(file.path(out, "FigS1_CV_grouped_comparison.pdf"), width = 5.8, height = 3.8); print(p1); dev.off()

## ---------------- 图S-2:逐折 AUC 散点 + 均值 ----------------
pf[, Model := factor(model, levels = lvl)]
mu <- pf[, .(m = mean(AUC)), by = .(Model)]
p2 <- ggplot(pf, aes(x = Model, y = AUC)) +
  geom_point(aes(colour = Model), size = 1.9, alpha = 0.75, stroke = 0,
             position = position_jitter(width = 0.14, height = 0, seed = 7)) +
  geom_segment(data = mu, aes(x = as.numeric(Model) - 0.30, xend = as.numeric(Model) + 0.30,
                              y = m, yend = m), inherit.aes = FALSE,
               colour = "grey15", linewidth = 0.55) +
  scale_colour_manual(values = c(RF="#E41A1C", SVM="#377EB8", XGB="#4DAF4A", GLM="#984EA3",
                                 GBM="#FF7F00", KNN="#E6AB02", NNET="#A65628",
                                 LASSO="#F781BF", DT="#666666")) +
  scale_y_continuous(limits = c(0.88, 1.008), breaks = seq(0.90, 1.00, 0.02)) +
  labs(x = NULL, y = "AUC",
       title = "Fold-level stability of compound-grouped cross-validation",
       subtitle = "10 points per model (5 folds \u00d7 2 repeats); bar = mean") +
  th + theme(legend.position = "none",
             axis.text.x = element_text(angle = 0))

pdf(file.path(out, "FigS2_groupedCV_perfold.pdf"), width = 5.4, height = 3.4); print(p2); dev.off()

sz <- file.info(file.path(out, c("FigS1_CV_grouped_comparison.pdf", "FigS2_groupedCV_perfold.pdf")))$size
print(setNames(round(sz), c("FigS1", "FigS2")))
cat(ifelse(all(sz > 4000), "OK\n", "WARN: 文件过小,可能空白\n"))
cat("fold AUC range:", sprintf("%.3f - %.3f", min(pf$AUC), max(pf$AUC)), "\n")
