# fig3H.R —— Fig. 3H:大鼠队列 CCNB2 诱导在【时间点】与【溶媒】两个分层内的一致性
# 风格沿用 fig_style2.R(th 主题);数据源为 pipeline 自身的输出,不重新拟合任何模型。
suppressMessages({library(data.table); library(ggplot2); library(pROC)})
setwd("D:/MLroute2"); out <- "outputs/paper"

th <- theme_classic(base_size = 10) +
  theme(axis.line = element_line(linewidth = 0.4, colour = "grey20"),
        axis.ticks = element_line(linewidth = 0.4, colour = "grey20"),
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5, colour = "grey30"),
        legend.key.size = unit(0.35, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 9.5, face = "bold", hjust = 0))

## ---- 1. 取 CCNB2 表达 + 样本元数据 ----
x   <- readRDS("outputs/rat_expr45.rds")
gi  <- which(toupper(rownames(x$X)) == "CCNB2"); stopifnot(length(gi) == 1)
v   <- as.numeric(x$X[gi, ])
gsm <- as.character(x$gsm)
sam <- read.delim("inputs/rat_samples.tsv", stringsAsFactors = FALSE)
mi  <- match(gsm, sam$acc); stopifnot(all(!is.na(mi)))
d <- data.table(gsm = gsm, expr = v,
                dose = sam$dose[mi], time = sam$time[mi], vehicle = sam$vehicle[mi])
d[, grp := ifelse(dose == "0 mg/kg", "Control", "Treated")]
# 组别顺序:分组折里对 Control 用的是投药批次键,这里同一分层内 Control 即该分层的空白对照

## ---- 2. 两个分层:时间点 / 溶媒 ----
mk <- function(key, levels, panel) {
  rbindlist(lapply(levels, function(k) {
    dd <- d[get(key) == k & grp %in% c("Control", "Treated")]
    if (nrow(dd) < 4 || length(unique(dd$grp)) < 2) return(NULL)
    r  <- roc(as.integer(dd$grp == "Treated"), dd$expr, quiet = TRUE, direction = "<")
    ci <- as.numeric(ci.auc(r))
    dd[, `:=`(stratum = k, panel = panel)]
    attr(dd, "stat") <- data.table(stratum = k, panel = panel,
      nT = sum(dd$grp == "Treated"), nC = sum(dd$grp == "Control"),
      log2FC = mean(dd$expr[dd$grp == "Treated"]) - mean(dd$expr[dd$grp == "Control"]),
      AUC = as.numeric(auc(r)), lo = ci[1], hi = ci[3])
    dd
  }), fill = TRUE)
}
## 每个分层单取统计量(避免 attr 在 rbindlist 中丢失)
stat_of <- function(key, levels, panel) rbindlist(lapply(levels, function(k) {
  dd <- d[get(key) == k & grp %in% c("Control", "Treated")]
  if (nrow(dd) < 4 || length(unique(dd$grp)) < 2) return(NULL)
  r <- roc(as.integer(dd$grp == "Treated"), dd$expr, quiet = TRUE, direction = "<")
  ci <- as.numeric(ci.auc(r))
  data.table(stratum = k, panel = panel, nT = sum(dd$grp == "Treated"),
             nC = sum(dd$grp == "Control"),
             log2FC = mean(dd$expr[dd$grp == "Treated"]) - mean(dd$expr[dd$grp == "Control"]),
             AUC = as.numeric(auc(r)), lo = ci[1], hi = ci[3])
}))
tlev <- c("0.25 d", "1 d", "3 d", "5 d", "7 d")
vlev <- c("Corn Oil", "CMC", "Water", "Saline")

dst <- rbind(stat_of("time", tlev, "By exposure duration"),
             stat_of("vehicle", vlev, "By control vehicle"))
fwrite(dst, file.path(out, "Fig3H_stats.txt"), sep = "\t")
print(dst)

## ---- 3. 绘图 ----
dat <- rbind(d[time %in% tlev, .(stratum = time, panel = "By exposure duration",
                                 grp, expr)],
             d[vehicle %in% vlev, .(stratum = vehicle, panel = "By control vehicle",
                                    grp, expr)])
dat[, stratum := factor(stratum, levels = c(tlev, vlev))]
dat[, grp := factor(grp, levels = c("Control", "Treated"))]
dat[, panel := factor(panel, levels = c("By exposure duration", "By control vehicle"))]

ann <- copy(dst)
ann[, stratum := factor(stratum, levels = c(tlev, vlev))]
ann[, panel := factor(panel, levels = c("By exposure duration", "By control vehicle"))]
ann[, lab := sprintf("log2FC %+0.2f\nAUC %.2f", log2FC, AUC)]

p <- ggplot(dat, aes(x = stratum, y = expr)) +
  geom_boxplot(aes(fill = grp), width = 0.62, outlier.shape = NA, linewidth = 0.32,
               colour = "grey25") +
  geom_jitter(aes(colour = grp), width = 0.16, size = 0.28, alpha = 0.30) +
  geom_text(data = ann, aes(x = stratum, y = Inf, label = lab), vjust = 1.15,
            size = 2.35, colour = "grey25", lineheight = 0.95, inherit.aes = FALSE) +
  facet_grid(panel ~ ., scales = "free", space = "free") +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.16))) +
  scale_fill_manual(values = c("Control" = "grey88", "Treated" = "#C0504D")) +
  scale_colour_manual(values = c("Control" = "grey45", "Treated" = "#8E2A24")) +
  labs(x = NULL, y = expression(log[2]~"CCNB2 expression"),
       title = "Consistency of CCNB2 induction in the rat screening cohort",
       subtitle = "Treated vs matched controls within each stratum; rat GSE57815",
       fill = NULL, colour = NULL) +
  th + theme(legend.position = c(0.99, 0.99), legend.justification = c(1, 1),
             legend.background = element_rect(fill = alpha("white", 0.8), colour = NA),
             panel.spacing = unit(0.35, "lines"))
pdf(file.path(out, "Fig3H_rat_CCNB2_strata.pdf"), width = 6.4, height = 4.6)
print(p); dev.off()
cat("\nwritten: Fig3H_rat_CCNB2_strata.pdf  size:",
    file.info(file.path(out, "Fig3H_rat_CCNB2_strata.pdf"))$size, "\n")
