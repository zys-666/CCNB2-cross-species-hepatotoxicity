# driver_rerun.R : 修正标签(数值dose==0 → 279/1939)后的大鼠侧全管线重跑
suppressMessages({library(data.table)})
steps <- c(
  "01_prep.R",            # 重build rat_expr45.rds (L1/L3 修正)
  "02_analysis.R",        # 方向表(3基因, 诊断用)
  "06_ccnb2.R",           # CCNB2 单基因全分层(大鼠新标签 + 小鼠)
  "09_fish_summary.R",    # 三物种表 + 森林图(动态读取06结果)
  "training_v2.R",        # 正式脚本: 9模型CV/重要性/排名/图 + models_for_figures.rds (~40min)
  "fig_importance_bars.R",# importance.pdf 柱状版(读新rds)
  "fig_shap_polar.R",     # SHAP 极坐标全套 9 模型(读新rds)
  "nested_cv.R",          # 嵌套CV (读新 rat_expr45.rds)
  "calibration.R",        # 校准图(读新 models_for_figures.rds)
  "youden_ccnb2.R",       # Youden 大鼠+小鼠(读新 rat_expr45.rds)
  "appendix_tables.R",    # 补充表 S1/S2(大鼠标签计数修正)
  "format_supp_tables.R") # 正式表 S3/S4/S6(从新 v2 输出生成)
t0 <- Sys.time()
for (s in steps) {
  cat("\n################################################\n")
  cat("STEP:", s, "| start:", format(Sys.time(), "%H:%M:%S"), "\n")
  cat("################################################\n")
  ok <- tryCatch({ source(file.path("D:/MLroute2", s), local = FALSE, echo = FALSE); TRUE },
                 error = function(e) { cat("\n!!! STEP FAILED:", s, "|", conditionMessage(e), "\n"); FALSE })
  if (!ok) { cat("ABORT at", s, "\n"); quit(status = 1) }
  cat("STEP DONE:", s, "|", format(Sys.time(), "%H:%M:%S"), "\n")
}
cat("\nALL STEPS COMPLETED in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
