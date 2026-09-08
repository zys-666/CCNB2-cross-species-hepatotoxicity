# calibration.R : reliability curves from internal-CV OOF predictions (stored models)
suppressMessages({library(ggplot2); library(dplyr)})
setwd("D:/MLroute2")
out <- "outputs/v2"
rds <- readRDS(file.path(out, "models_for_figures.rds"))
my_colors <- c(RF="#E41A1C", SVM="#377EB8", XGB="#4DAF4A", GLM="#984EA3", GBM="#FF7F00",
               KNN="#E6AB02", NNET="#A65628", LASSO="#F781BF", DT="#666666")

preds <- lapply(names(rds$mods), function(nm) {
  p <- rds$mods[[nm]]$pred
  data.frame(model = nm, p = p$Treat, y = as.integer(p$obs == "Treat"))
})
preds[["XGB"]] <- data.frame(model = "XGB", p = rds$oof_xgb_plot, y = rds$y_bin)
df <- do.call(rbind, preds)

ece_tab <- lapply(split(df, df$model), function(d) {
  qs <- unique(quantile(d$p, probs = seq(0, 1, 0.1), na.rm = TRUE))
  d$bin <- cut(d$p, breaks = qs, include.lowest = TRUE)
  b <- d %>% group_by(bin) %>%
    summarise(pred = mean(p), obs = mean(y), n = n(), .groups = "drop") %>% na.omit()
  ece <- sum(b$n * abs(b$obs - b$pred)) / sum(b$n)
  data.frame(model = d$model[1], ECE = round(ece, 4), bins = nrow(b))
})
ece <- do.call(rbind, ece_tab)
write.table(ece, file.path(out, "calibration_summary.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("ECE per model:\n"); print(ece)

bin_all <- do.call(rbind, lapply(split(df, df$model), function(d) {
  qs <- unique(quantile(d$p, probs = seq(0, 1, 0.1), na.rm = TRUE))
  b <- d %>% mutate(bin = cut(p, breaks = qs, include.lowest = TRUE)) %>%
    group_by(bin) %>% summarise(pred = mean(p), obs = mean(y), n = n(), .groups = "drop") %>% na.omit()
  b$model <- d$model[1]; b
}))
bin_all <- merge(bin_all, ece, by = "model")
bin_all$lab <- paste0(bin_all$model, " (ECE=", sprintf("%.3f", bin_all$ECE), ")")

p <- ggplot(bin_all, aes(x = pred, y = obs, colour = model, group = model)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey40", linewidth = 0.6) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  scale_colour_manual(values = my_colors, breaks = names(my_colors)) +
  labs(title = "Calibration curves (internal CV out-of-fold predictions)",
       subtitle = "Ideal calibration = dashed diagonal; legend shows Expected Calibration Error",
       x = "Predicted probability (Treat)", y = "Observed frequency") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right", legend.title = element_blank())
pdf(file.path(out, "calibration_CV.pdf"), width = 9, height = 6.5)
print(p); dev.off()
cat("saved calibration_CV.pdf | size:", file.info(file.path(out, "calibration_CV.pdf"))$size, "\n")
