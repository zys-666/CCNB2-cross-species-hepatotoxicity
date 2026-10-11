# 09_fish_summary.R : three-species CCNB2 evidence summary + forest plot (dynamic parse from results_ccnb2.txt)
# 鱼侧数据来源: D:/user data/Desktop/轮胎/5.机器学习/基因表达量/基因表达量.xlsx (H_vs_C 块, 已人工核对)
suppressMessages({library(data.table)})
setwd("D:/MLroute2")
dir.create("outputs/paper", showWarnings = FALSE, recursive = TRUE)

# fish CCNB2 (FPKM, n=3/3) from xlsx: C: 2.385 2.072 2.21 ; H: 56.196 65.06 69.061
C <- c(2.385, 2.072, 2.210); H <- c(56.196, 65.060, 69.061)
tt <- t.test(log2(H), log2(C))
fish_d <- mean(log2(H)) - mean(log2(C)); fish_ci <- tt$conf.int

# 动态读取 results_ccnb2.txt 中的大鼠 L1 与小鼠 all 两行（06_ccnb2.R 每次重跑后自动更新）
rc <- readLines("outputs/results_ccnb2.txt")
parse_row <- function(line) {
  m <- regexec("nT= *([0-9]+) nC= *([0-9]+).*d= *(-?[0-9.]+) \\(95%CI (-?[0-9.]+)-(-?[0-9.]+)\\) p=([0-9.eE+-]+)", line)
  v <- regmatches(line, m)[[1]]
  if (length(v) != 7) stop("parse failed: ", substr(line, 1, 80))
  list(nT = as.integer(v[2]), nC = as.integer(v[3]), d = as.numeric(v[4]),
       lo = as.numeric(v[5]), hi = as.numeric(v[6]), p = as.numeric(v[7]))
}
rat <- parse_row(rc[grep("^rat L1 ", rc)])
mus <- parse_row(rc[grep("^mouse all ", rc)])

tab <- data.frame(
  Species_Cohort = c("Crucian carp (NTH H vs C, RNA-seq)",
                     "Rat (GSE57815, 200 compounds vs vehicle)",
                     "Mouse (GSE44783, 13 compounds vs vehicle)"),
  n_Treat = c(3L, rat$nT, mus$nT), n_Control = c(3L, rat$nC, mus$nC),
  delta_log2 = c(fish_d, rat$d, mus$d),
  CI_lo = c(fish_ci[1], rat$lo, mus$lo), CI_hi = c(fish_ci[2], rat$hi, mus$hi),
  p = c(tt$p.value, rat$p, mus$p), stringsAsFactors = FALSE)
tab$delta_log2 <- round(tab$delta_log2, 3); tab$CI_lo <- round(tab$CI_lo, 3)
tab$CI_hi <- round(tab$CI_hi, 3); tab$p <- signif(tab$p, 2)
write.table(tab, "outputs/paper/CCNB2_three_species.txt", sep = "\t", quote = FALSE, row.names = FALSE)
print(tab)

png("outputs/paper/CCNB2_forest.png", width = 2400, height = 1400, res = 300)
par(mar = c(4.5, 13, 2.5, 1))
y <- 3:1
plot(NA, xlim = c(0, 6.2), ylim = c(0.5, 3.5), yaxt = "n", ylab = "", xlab = "",
     main = "CCNB2 upregulation across species (log2 scale)")
abline(v = 0, lty = 2, col = "grey50")
for (i in 1:3) {
  xi <- tab$delta_log2[i]; lo <- tab$CI_lo[i]; hi <- tab$CI_hi[i]
  points(xi, y[i], pch = 18, cex = 1.8, col = c("#c0392b", "#2471a3", "#1e8449")[i])
  segments(lo, y[i], hi, y[i], lwd = 2.5, col = c("#c0392b", "#2471a3", "#1e8449")[i])
}
axis(2, at = y, labels = tab$Species_Cohort, las = 1, cex.axis = 0.62)
axis(1); mtext("log2(mean treated / mean control)", side = 1, line = 2.6)
text(tab$delta_log2 + 0.18, y, labels = paste0(sprintf("%+.2f", tab$delta_log2),
     "  (95%CI ", sprintf("%.2f", tab$CI_lo), "-", sprintf("%.2f", tab$CI_hi), ")"),
     cex = 0.58, pos = 4)
dev.off()
cat("wrote outputs/paper/CCNB2_three_species.txt and CCNB2_forest.png\n")
