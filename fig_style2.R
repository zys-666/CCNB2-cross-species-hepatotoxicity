# fig_style2.R : 三张图改为期刊风格（矢量 PDF）→ S1 CV性能(点-误差棒) / S2 森林图 / S3 鱼剂量响应
suppressMessages({library(data.table); library(ggplot2)})
setwd("D:/MLroute2"); out <- "outputs/paper"
th <- theme_classic(base_size = 10) +
  theme(axis.line = element_line(linewidth = 0.4, colour = "grey20"),
        axis.ticks = element_line(linewidth = 0.4, colour = "grey20"),
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5, colour = "grey30"),
        legend.key.size = unit(0.35, "cm"))

## ---------------- S1：九模型 5x3 CV vs 嵌套 CV（点-误差棒） ----------------
cv  <- fread("outputs/v2/modelPerf_CV.txt")
ncv <- fread("outputs/v2/nestedCV_perf.txt")
m <- merge(cv[, .(Model, CV_AUC_mean, CV_AUC_sd)], ncv[, .(Model, NestedAUC_mean, NestedAUC_sd)], by = "Model")
m[, ordv := NestedAUC_mean]; m <- m[order(-NestedAUC_mean)]
lvl <- rev(m$Model)
d <- rbind(data.table(Model = m$Model, metric = "5\u00d73 CV",    auc = m$CV_AUC_mean,     sd = m$CV_AUC_sd),
           data.table(Model = m$Model, metric = "Nested CV", auc = m$NestedAUC_mean, sd = m$NestedAUC_sd))
d[, Model := factor(Model, levels = lvl)]
d[, metric := factor(metric, levels = c("5\u00d73 CV", "Nested CV"))]
p1 <- ggplot(d, aes(x = auc, y = Model, colour = metric, shape = metric)) +
  geom_errorbarh(aes(xmin = auc - sd, xmax = auc + sd), height = 0.18, linewidth = 0.45,
                 position = position_dodge(width = 0.45)) +
  geom_point(size = 2.1, position = position_dodge(width = 0.45)) +
  scale_colour_manual(values = c("5\u00d73 CV" = "#3B6EA5", "Nested CV" = "#C0504D")) +
  scale_shape_manual(values = c("5\u00d73 CV" = 16, "Nested CV" = 17)) +
  scale_x_continuous(limits = c(0.925, 1.002), breaks = seq(0.93, 1.00, 0.02)) +
  labs(x = "AUC (mean \u00b1 SD)", y = NULL, colour = NULL, shape = NULL,
       title = "Internal discrimination of nine classifiers",
       subtitle = "Rat GSE57815; 5\u00d73 CV vs fully nested CV") +
  th + theme(legend.position = c(0.02, 0.06), legend.justification = c(0, 0),
             legend.background = element_rect(fill = alpha("white", 0.75), colour = NA))
pdf(file.path(out, "FigS_CV_performance.pdf"), width = 5.4, height = 3.6); print(p1); dev.off()

## ---------------- S2：三物种森林图（期刊风格） ----------------
ts <- fread(file.path(out, "CCNB2_three_species.txt"))
ts[, cohort := c("Crucian carp (NTH  H vs C)", "Rat (GSE57815)", "Mouse (GSE44783)")]
ts[, nlab := c("n = 3 vs 3", "n = 1,939 vs 279", "n = 330 vs 120")]
ts[, cohort := factor(cohort, levels = rev(ts$cohort))]
p2 <- ggplot(ts, aes(x = delta_log2, y = cohort, colour = cohort)) +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey55", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = CI_lo, xmax = CI_hi), height = 0.10, linewidth = 0.6, colour = "grey25") +
  geom_point(size = 2.6) +
  geom_text(aes(label = sprintf("%+.2f  (%.2f to %.2f)", delta_log2, CI_lo, CI_hi)),
            hjust = -0.14, size = 2.9, colour = "grey15", show.legend = FALSE) +
  geom_text(aes(label = nlab, x = -0.45), hjust = 1.05, size = 2.7, colour = "grey45", show.legend = FALSE) +
  scale_colour_manual(values = c("Crucian carp (NTH  H vs C)" = "#B03A2E",
                                 "Rat (GSE57815)" = "#1F618D", "Mouse (GSE44783)" = "#196F3D")) +
  scale_x_continuous(limits = c(-1.5, 6.6), breaks = c(0, 2, 4, 6)) +
  labs(x = expression(log[2]~"fold change of CCNB2 (treated vs control)"), y = NULL,
       title = "CCNB2 up-regulation across three cohorts\n(one discovery cohort and two independent rodent cohorts)") +
  th + theme(legend.position = "none", axis.line.y = element_blank(),
             axis.ticks.y = element_blank())
pdf(file.path(out, "FigS_forest.pdf"), width = 6.4, height = 2.7); print(p2); dev.off()

## ---------------- S3：鱼 CCNB2 剂量响应（期刊风格） ----------------
g <- data.table(group = factor(rep(c("Control", "Low NTH", "High NTH"), each = 3),
                               levels = c("Control", "Low NTH", "High NTH")),
                fpkm = c(2.385, 2.072, 2.210, 0.811, 0.929, 0.215, 56.196, 65.060, 69.061))
sm <- g[, .(mean = mean(fpkm), sd = sd(fpkm)), by = group]
p3 <- ggplot(g, aes(x = group, y = fpkm)) +
  geom_col(data = sm, aes(x = group, y = mean), width = 0.55, fill = "grey88",
           colour = "grey30", linewidth = 0.35, inherit.aes = FALSE) +
  geom_errorbar(data = sm, aes(x = group, ymin = mean - sd, ymax = mean + sd), width = 0.13,
                linewidth = 0.4, colour = "grey25", inherit.aes = FALSE) +
  geom_jitter(width = 0.08, size = 1.7, shape = 21, fill = "white", colour = "grey15", stroke = 0.6) +
  scale_y_log10(limits = c(0.12, 320), breaks = c(0.2, 1, 5, 20, 60),
                labels = c("0.2", "1", "5", "20", "60")) +
  annotation_logticks(sides = "l", linewidth = 0.3, colour = "grey40") +
  annotate("segment", x = 1, xend = 2, y = 3.3, yend = 3.3, linewidth = 0.35) +
  annotate("text", x = 1.5, y = 4.2, label = "*", size = 4.5) +
  annotate("segment", x = 1, xend = 3, y = 145, yend = 145, linewidth = 0.35) +
  annotate("text", x = 2, y = 185, label = "*", size = 4.5) +
  annotate("text", x = 2, y = 0.16, size = 2.7, colour = "grey30",
           label = "Low: log2FC = -1.77, FDR = 8.9e-07") +
  annotate("text", x = 2, y = 0.135, size = 2.7, colour = "grey30",
           label = "High: log2FC = +4.84, FDR = 7.4e-67") +
  labs(x = NULL, y = "CCNB2 expression (FPKM, log scale)",
       title = "Crucian carp liver: CCNB2 dose response") +
  th + theme(axis.line.x = element_line(linewidth = 0.4, colour = "grey20"))
pdf(file.path(out, "FigS_fish_dose.pdf"), width = 4.0, height = 4.2); print(p3); dev.off()

cat("written: FigS_CV_performance.pdf / FigS_forest.pdf / FigS_fish_dose.pdf\n")
cat("sizes:", sapply(c("FigS_CV_performance.pdf","FigS_forest.pdf","FigS_fish_dose.pdf"),
                    function(f) file.info(file.path(out, f))$size), "\n")
