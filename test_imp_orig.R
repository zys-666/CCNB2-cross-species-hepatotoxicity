suppressMessages({library(DALEX)})
setwd("D:/MLroute2")
rds <- readRDS("outputs/v2/models_for_figures.rds")
# replicate THEIR original importance plot code verbatim (their row slicing, their model order)
data <- rds$data                       # 45 genes + Type -> ncol = 46
cat("ncol(data):", ncol(data), "\n")
importance_rf    <- rds$vi_raw[["RF"]]; importance_svm <- rds$vi_raw[["SVM"]]
importance_xgb   <- rds$vi_raw[["XGB"]]; importance_gbm <- rds$vi_raw[["GBM"]]
importance_knn   <- rds$vi_raw[["KNN"]]; importance_nnet <- rds$vi_raw[["NNET"]]
importance_lasso <- rds$vi_raw[["LASSO"]]; importance_dt <- rds$vi_raw[["DT"]]
importance_glm   <- rds$vi_raw[["GLM"]]
cat("rows per importance obj:", nrow(importance_rf), "| class:", class(importance_rf)[1], "\n")
pdf("outputs/v2/_imp_orig.pdf", width = 7, height = 12)
tryCatch({
  plot(
    importance_rf[c(1,(ncol(data)-9):(ncol(data)+1)),],
    importance_svm[c(1,(ncol(data)-9):(ncol(data)+1)),],
    importance_xgb[c(1,(ncol(data)-9):(ncol(data)+1)),],
    importance_gbm[c(1,(ncol(data)-9):(ncol(data)+1)),],
    importance_knn[c(1,(ncol(data)-9):(ncol(data)+1)),],
    importance_nnet[c(1,(ncol(data)-9):(ncol(data)+1)),],
    importance_lasso[c(1,(ncol(data)-9):(ncol(data)+1)),],
    importance_dt[c(1,(ncol(data)-9):(ncol(data)+1)),],
    importance_glm[c(1,(ncol(data)-9):(ncol(data)+1)),]
  )
}, error = function(e) cat("ERR:", conditionMessage(e), "\n"))
dev.off()
cat("size:", file.info("outputs/v2/_imp_orig.pdf")$size, "\n")
