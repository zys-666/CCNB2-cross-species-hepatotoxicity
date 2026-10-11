## test_caret_foldwise.R — 判定 caret 的 preProcess 是否在每个折内拟合
## 设计:toy 数据 + preProcess="range";取第一个留出折,比较两种缩放下同一 glm 的预测
##   (a) range 由【全量 200 样本】估计 → 若有泄漏,caret 的折外预测应与 (a) 一致
##   (b) range 由【训练折 159 样本】估计 → 若 caret 正确,应与 (b) 一致
suppressMessages(library(caret))
set.seed(1)
n <- 200
x1 <- rnorm(n, 5, 1); x2 <- rnorm(n, 100, 20)
y <- factor(ifelse(x1 + x2 / 25 + rnorm(n) > 9, "b", "a"), levels = c("a", "b"))
d <- data.frame(x1 = x1, x2 = x2, y = y)

ctrl <- trainControl(method = "cv", number = 5, classProbs = TRUE,
                     savePredictions = TRUE, summaryFunction = twoClassSummary)
set.seed(2)
m <- train(y ~ ., data = d, method = "glm", family = "binomial",
           metric = "ROC", trControl = ctrl, preProcess = "range")

pr <- m$pred
pr <- pr[pr$Resample == pr$Resample[1], ]       # 取第一个折(每样本只出现一次)
te <- pr$rowIndex
tr <- setdiff(seq_len(n), te)
cat("留出折:", pr$Resample[1], "| 测试", length(te), "训练", length(tr), "\n")

pred_case <- function(fit_on) {
  rng <- apply(as.matrix(d[fit_on, c("x1", "x2")]), 2, range)
  sc  <- function(mm) sweep(sweep(mm, 2, rng[1, ], "-"), 2, rng[2, ] - rng[1, ], "/")
  trd <- as.data.frame(sc(as.matrix(d[tr, c("x1", "x2")]))); trd$y <- d$y[tr]
  ted <- as.data.frame(sc(as.matrix(d[te, c("x1", "x2")])))
  g <- glm(y ~ ., data = trd, family = binomial)
  predict(g, newdata = ted, type = "response")
}
pa <- pred_case(seq_len(n))        # (a) 全量估计 range
pb <- pred_case(tr)                # (b) 训练折估计 range

da <- max(abs(pa - pr$b)); db <- max(abs(pb - pr$b))
cat(sprintf("\n(a) 全量 range   与 caret 折外预测 max|diff| = %.8f\n", da))
cat(sprintf("(b) 训练折 range 与 caret 折外预测 max|diff| = %.8f\n", db))
cat("\n结论:",
    if (db < 1e-8) "caret 在【折内】拟合预处理 ✅(无泄漏)"
    else if (da < 1e-8) "caret 用【全量】拟合预处理 ❌(有泄漏)"
    else "两者都不完全一致,需进一步核查", "\n")
cat("\ncaret 为最终模型保存的 range(全量拟合,用于最终模型):\n"); print(m$preProcess$ranges)
