suppressMessages({library(data.table)})
setwd("D:/MLroute2")
r <- fread("inputs/rat_samples.tsv")
r[, dose := trimws(dose)]
r[, dose0_bug := grepl("^0", dose)]
r[, num := suppressWarnings(as.numeric(sub(" .*", "", dose)))]
r[, dose0_ok := !is.na(num) & num == 0]
cat("=== 剂量文本里以 0 开头但数值 >0 的样本（被 bug 误标为 Control 的候选）===\n")
bad <- r[dose0_bug == TRUE & dose0_ok == FALSE]
print(bad[, .N, by = .(dose)][order(-N)])
cat("\n误标候选总数:", nrow(bad), "\n")
cat("\n=== 对照侧(dose0_bug)总数:", sum(r$dose0_bug),
    "| 修正后真正 dose==0:", sum(r$dose0_ok),
    "| 差异:", sum(r$dose0_bug) - sum(r$dose0_ok), "\n")
cat("\n=== 全部 dose 文本取值一览 ===\n")
print(r[, .N, by = dose][order(dose)])
cat("\n=== 修正规则下 L1 计数 ===\n")
cat("Control:", sum(r$dose0_ok), " Treat:", nrow(r) - sum(r$dose0_ok), "\n")
# 误标样本的溶媒/时间分布
cat("\n=== 误标候选 by vehicle ===\n")
print(bad[, .N, by = .(vehicle, time)][order(vehicle)])
