## verify_resample.R — 验证 caret $resample 是否只反映【最终调参配置】
## 若 $resample$ROC 与"只用 bestTune 的 $pred 逐折算出的 AUC"一致 → $resample 未被污染
suppressMessages({library(pROC); library(data.table)})
r <- readRDS("D:/MLroute2/outputs/v2/models_for_figures.rds")
for (nm in c("RF", "LASSO", "GBM", "DT")) {
  p  <- as.data.table(r$mods[[nm]]$pred)
  bt <- as.data.table(r$mods[[nm]]$bestTune)
  n_all <- nrow(p)
  for (cc in intersect(colnames(bt), colnames(p))) p <- p[get(cc) == bt[[cc]][1]]
  ## 只用最终配置的预测,逐 Resample 算 AUC
  byres <- p[, .(auc_final = as.numeric(auc(roc(obs, Treat, quiet = TRUE,
                                                levels = c("Control","Treat"), direction = "<")))), by = Resample]
  rs <- as.data.table(r$mods[[nm]]$resample)[, .(Resample, ROC_resample = ROC)]
  m <- merge(byres, rs, by = "Resample")
  cat(sprintf("%-6s rows %6d -> %5d (final) | fold AUC: max|diff| = %.4f ; corr = %.4f\n",
              nm, n_all, nrow(p), max(abs(m$auc_final - m$ROC_resample)), cor(m$auc_final, m$ROC_resample)))
  cat("       resample mean =", round(mean(rs$ROC_resample), 4),
      "| final-config per-fold mean =", round(mean(byres$auc_final), 4), "\n")
}
