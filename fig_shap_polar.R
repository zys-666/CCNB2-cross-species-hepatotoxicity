# fig_shap_polar.R —— 按 旋图.txt(500-1083行) 的 fastshap 极坐标蜂巢图全套，
# 数据源改为训练侧：background = 修正标签全量训练矩阵, shap_sample = 训练侧 40 样本(seed 123)
# 用法: Rscript fig_shap_polar.R [RF|SVM|...]   (不传则跑全部9个)
suppressMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(tibble); library(stringr)
  library(scales); library(purrr); library(fastshap)
  library(randomForest); library(kernlab); library(gbm); library(nnet)
  library(glmnet); library(rpart); library(xgboost)   # predict 方法注册
})
setwd("D:/MLroute2")
out <- "outputs/v2"
rds <- readRDS(file.path(out, "models_for_figures.rds"))
data <- rds$data
genes <- setdiff(colnames(data), "Type")
x_train <- as.data.frame(data[, genes, drop = FALSE])   # 背景 = 全量修正标签训练数据

model_list <- list(
  RF    = rds$mods$RF$finalModel,    SVM = rds$mods$SVM$finalModel,
  GLM   = rds$mods$GLM$finalModel,   GBM = rds$mods$GBM$finalModel,
  KNN   = rds$mods$KNN$finalModel,   NNET = rds$mods$NNET$finalModel,
  LASSO = rds$mods$LASSO$finalModel, DT = rds$mods$DT$finalModel,
  XGB   = rds$mod_xgb)
pred_fun_list <- list(
  RF    = function(m, x) as.numeric(predict(m, x, type = "prob")[, 2]),
  SVM   = function(m, x) as.numeric(predict(m, x, type = "probabilities")[, 2]),
  GLM   = function(m, x) as.numeric(predict(m, x, type = "response")),
  GBM   = function(m, x) as.numeric(predict(m, x, type = "response")),
  KNN   = function(m, x) as.numeric(predict(m, x, type = "prob")[, 2]),
  NNET  = function(m, x) {
    p <- predict(m, x, type = "raw")
    if (is.matrix(p)) {
      cn <- colnames(p)
      if (!is.null(cn) && "Treat" %in% cn) as.numeric(p[, "Treat"]) else as.numeric(p[, ncol(p)])
    } else as.numeric(p)
  },
  LASSO = function(m, x) as.numeric(predict(m, as.matrix(x), type = "response")[, 2]),
  DT    = function(m, x) as.numeric(predict(m, x, type = "prob")[, 2]),
  XGB   = function(m, x) as.numeric(predict(m, as.matrix(x))))

rescale_clip <- function(x, q_low = 0.05, q_high = 0.95) {
  x <- as.numeric(x); qs <- quantile(x, probs = c(q_low, q_high), na.rm = TRUE, names = FALSE)
  lo <- qs[1]; hi <- qs[2]
  if (abs(hi - lo) < 1e-12) { qs <- quantile(x, probs = c(0.01, 0.99), na.rm = TRUE, names = FALSE)
    lo <- qs[1]; hi <- qs[2]
    if (abs(hi - lo) < 1e-12) { lo <- min(x, na.rm = TRUE); hi <- max(x, na.rm = TRUE) } }
  if (abs(hi - lo) < 1e-12) return(rep(0.5, length(x)))
  out <- (x - lo) / (hi - lo); scales::squish(out, range = c(0, 1))
}
simple_beeswarm <- function(y, nbins = NULL, width = 0.12) {
  y <- as.numeric(y); n <- length(y)
  if (n <= 1) return(rep(0, n))
  if (is.null(nbins)) nbins <- max(floor(n / 6), 1)
  nbins <- max(1L, as.integer(nbins))
  if (diff(range(y, na.rm = TRUE)) < 1e-12) {
    step <- width / max(ceiling(n / 2), 1); pos <- numeric(n)
    if (n > 1) { vals <- rep(seq_len(ceiling((n - 1) / 2)), each = 2)[seq_len(n - 1)]
      sgn <- rep(c(1, -1), length.out = n - 1); pos[-1] <- vals * sgn }
    return(pos * step)
  }
  bins <- cut(y, breaks = nbins, include.lowest = TRUE, labels = FALSE)
  xoff <- numeric(n); nmax <- max(tabulate(bins, nbins))
  step <- width / max(ceiling(nmax / 2), 1)
  pattern_fun <- function(k) { out <- numeric(k)
    if (k > 1) { vals <- rep(seq_len(ceiling((k - 1) / 2)), each = 2)[seq_len(k - 1)]
      sgn <- rep(c(1, -1), length.out = k - 1); out[-1] <- vals * sgn }; out }
  for (b in sort(unique(bins))) { idx <- which(bins == b); idx <- idx[order(y[idx])]
    xoff[idx] <- pattern_fun(length(idx)) * step }
  xoff
}
circle_df <- function(r, n = 500, group = "circle") {
  theta <- seq(0, 2 * pi, length.out = n)
  tibble(group = group, x = r * sin(theta), y = r * cos(theta))
}
wedge_df <- function(theta, width, r1, group, r0 = 0, n = 120) {
  ang <- seq(theta - width / 2, theta + width / 2, length.out = n)
  tibble(group = group,
         x = c(r0 * sin(theta - width / 2), r1 * sin(ang), r0 * sin(theta + width / 2)),
         y = c(r0 * cos(theta - width / 2), r1 * cos(ang), r0 * cos(theta + width / 2)))
}

top_n_plot <- 15; base_family <- "serif"; rose_power <- 1.0
shap_nsim <- 50; shap_sample_n <- 40

args <- commandArgs(trailingOnly = TRUE)
run_models <- if (length(args) > 0) args else names(model_list)

set.seed(123)
shap_sample_idx <- sample(seq_len(nrow(x_train)), min(shap_sample_n, nrow(x_train)))
shap_sample <- x_train[shap_sample_idx, , drop = FALSE]

for (nm in run_models) {
  cat("\n========== SHAP polar:", nm, "==========\n")
  current_model <- model_list[[nm]]; current_pred_fun <- pred_fun_list[[nm]]
  pred_wrapper_fastshap <- function(object, newdata) {
    newdata <- as.data.frame(newdata)
    as.numeric(current_pred_fun(object, newdata))
  }
  set.seed(123)
  shap_mat <- fastshap::explain(object = current_model, X = x_train,
                                newdata = shap_sample,
                                pred_wrapper = pred_wrapper_fastshap,
                                nsim = shap_nsim, adjust = TRUE)
  shap_mat <- as.matrix(shap_mat); colnames(shap_mat) <- colnames(shap_sample)
  cat("SHAP dim:", dim(shap_mat)[1], "x", dim(shap_mat)[2], "\n")

  shap_long_all <- as.data.frame(shap_mat, check.names = FALSE) %>% as_tibble() %>%
    mutate(row_id = row_number()) %>%
    pivot_longer(cols = -row_id, names_to = "feature", values_to = "shap") %>%
    left_join(as.data.frame(shap_sample, check.names = FALSE) %>% as_tibble() %>%
                mutate(row_id = row_number()) %>%
                pivot_longer(cols = -row_id, names_to = "feature", values_to = "feature_value"),
              by = c("row_id", "feature"))
  feature_importance_all <- shap_long_all %>%
    group_by(feature) %>%
    summarise(mean_abs_shap = mean(abs(shap), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_abs_shap))
  top_n_plot_actual <- min(top_n_plot, nrow(feature_importance_all))
  feature_levels <- feature_importance_all$feature[seq_len(top_n_plot_actual)]
  shap_long <- shap_long_all %>% filter(feature %in% feature_levels) %>%
    group_by(feature) %>% mutate(feature_value_norm = rescale_clip(feature_value)) %>% ungroup() %>%
    mutate(feature = factor(feature, levels = feature_levels))
  feature_importance <- feature_importance_all %>% filter(feature %in% feature_levels) %>%
    mutate(feature = factor(feature, levels = feature_levels))
  n_feat <- length(feature_levels)

  shap_min <- min(shap_long$shap, na.rm = TRUE); shap_max <- max(shap_long$shap, na.rm = TRUE)
  max_abs_shap <- max(abs(shap_long$shap), na.rm = TRUE)
  max_mean <- max(feature_importance$mean_abs_shap, na.rm = TRUE)
  angle_step <- 2 * pi / n_feat
  feature_angles <- tibble(feature = factor(feature_levels, levels = feature_levels),
                           theta = seq(0, 2 * pi - angle_step, length.out = n_feat))
  bar_scale <- if (is.finite(max_abs_shap) && max_abs_shap > 0) max_abs_shap * 0.88 else 1
  inner_gap <- max(max_abs_shap * 0.22, 0.12)
  shift <- if (shap_min < 0) abs(shap_min) + bar_scale + inner_gap else bar_scale + inner_gap
  zero_radius <- shift; r_min_data <- shift + shap_min; r_max_data <- shift + shap_max
  base_span <- max(max_abs_shap, bar_scale)
  outer_radius <- max(r_max_data, zero_radius) + base_span * 0.12
  label_radius <- outer_radius + base_span * 0.35
  mean_to_radius <- function(v) {
    v <- as.numeric(v)
    if (is.finite(max_mean) && max_mean > 0) (v / max_mean)^rose_power * bar_scale else rep(0, length(v))
  }
  feature_importance_polar <- feature_importance %>% select(feature, mean_abs_shap) %>%
    mutate(feature = factor(as.character(feature), levels = feature_levels)) %>%
    left_join(feature_angles, by = "feature") %>% mutate(bar_r = mean_to_radius(mean_abs_shap))
  bar_width <- angle_step * 0.82
  bar_df <- purrr::map_dfr(seq_len(nrow(feature_importance_polar)), function(i) {
    wedge_df(theta = feature_importance_polar$theta[i], width = bar_width,
             r1 = feature_importance_polar$bar_r[i], group = paste0("bar_", i), r0 = 0)
  })
  shap_polar_df <- shap_long %>%
    mutate(feature = factor(as.character(feature), levels = feature_levels)) %>%
    left_join(feature_angles, by = "feature") %>%
    group_by(feature) %>%
    mutate(theta_offset = simple_beeswarm(shap, nbins = 100, width = angle_step * 0.26),
           theta_plot = theta + theta_offset, r = shap + shift) %>% ungroup() %>%
    mutate(x = r * sin(theta_plot), y = r * cos(theta_plot))
  if (shap_min < 0) { neg_ticks <- pretty(c(shap_min, 0), n = 4); neg_ticks <- neg_ticks[neg_ticks < 0 & neg_ticks >= shap_min]
  } else neg_ticks <- numeric(0)
  if (shap_max > 0) { pos_ticks <- pretty(c(0, shap_max), n = 4); pos_ticks <- pos_ticks[pos_ticks > 0 & pos_ticks <= shap_max]
  } else pos_ticks <- numeric(0)
  full_ticks <- sort(unique(c(neg_ticks, 0, pos_ticks)))
  tick_circle_df <- purrr::map_dfr(seq_along(full_ticks), function(i)
    circle_df(r = shift + full_ticks[i], group = paste0("tick_circle_", i)))
  zero_circle_df  <- circle_df(zero_radius, group = "zero_circle")
  outer_circle_df <- circle_df(outer_radius, group = "outer_circle")
  spoke_df <- feature_angles %>%
    transmute(x = 0, y = 0, xend = outer_radius * sin(theta), yend = outer_radius * cos(theta))
  outer_tick_len <- base_span * 0.05
  outer_tick_df <- feature_angles %>% transmute(
    x = outer_radius * sin(theta), y = outer_radius * cos(theta),
    xend = (outer_radius + outer_tick_len) * sin(theta), yend = (outer_radius + outer_tick_len) * cos(theta))
  label_gap <- outer_radius * 0.08
  shap_tick_len <- base_span * 0.035
  empty_seg_df <- tibble(x = numeric(), y = numeric(), xend = numeric(), yend = numeric())
  empty_text_df <- tibble(x = numeric(), y = numeric(), label = character())
  neg_axis_df <- if (shap_min < 0) tibble(x = -zero_radius, y = 0, xend = -r_min_data, yend = 0) else empty_seg_df
  if (shap_min < 0) { neg_ticks <- sort(unique(c(round(shap_min, 2), round(shap_min / 2, 2), 0)))
    neg_ticks <- neg_ticks[neg_ticks < 0 & neg_ticks >= shap_min] } else neg_ticks <- numeric(0)
  neg_tick_df <- tibble(value = neg_ticks, x = -(shift + neg_ticks), y = 0)
  neg_tick_seg_df <- neg_tick_df %>% transmute(x = x, y = 0, xend = x, yend = -label_gap * 0.35)
  neg_label_df <- neg_tick_df %>% transmute(x = x, y = -label_gap * 1.1, label = sprintf("%.2f", value))
  neg_axis_title_df <- if (shap_min < 0) tibble(
    x = mean(c(-zero_radius, -r_min_data)), y = -label_gap * 2.6,
    label = "Negative SHAP\n(Decrease Treat Prob)") else empty_text_df
  pos_axis_df <- if (shap_max > 0) tibble(x = zero_radius, y = 0, xend = r_max_data, yend = 0) else empty_seg_df
  if (shap_max > 0) { pos_ticks <- sort(unique(c(0, round(shap_max / 2, 2), round(shap_max, 2))))
    pos_ticks <- pos_ticks[pos_ticks > 0 & pos_ticks <= shap_max] } else pos_ticks <- numeric(0)
  pos_tick_df <- tibble(value = pos_ticks, x = shift + pos_ticks, y = 0)
  pos_tick_seg_df <- pos_tick_df %>% transmute(x = x, y = 0, xend = x, yend = label_gap * 0.35)
  pos_label_df <- pos_tick_df %>% transmute(x = x, y = label_gap * 1.1, label = sprintf("%.2f", value))
  pos_axis_title_df <- if (shap_max > 0) tibble(
    x = mean(c(zero_radius, r_max_data)), y = label_gap * 2.6,
    label = "Positive SHAP\n(Increase Treat Prob)") else empty_text_df
  mean_ticks <- pretty(c(0, max_mean), n = 3); mean_ticks <- mean_ticks[mean_ticks >= 0 & mean_ticks <= max_mean]
  mean_tick_df <- tibble(value = mean_ticks, r = mean_to_radius(mean_ticks))
  mean_tick_len <- base_span * 0.03
  mean_axis_df <- tibble(x = 0, y = 0, xend = 0, yend = bar_scale)
  mean_tick_seg_df <- mean_tick_df %>% transmute(x = -mean_tick_len, y = r, xend = mean_tick_len, yend = r)
  mean_label_df <- mean_tick_df %>% transmute(x = -mean_tick_len * 1.8, y = r, label = sprintf("%.3f", value), hjust = 1)
  label_df <- feature_angles %>% mutate(
    x = label_radius * sin(theta), y = label_radius * cos(theta),
    label = stringr::str_wrap(as.character(feature), width = 14),
    hjust = case_when(abs(sin(theta)) < 0.12 ~ 0.5, sin(theta) > 0 ~ 0, TRUE ~ 1))

  p_polar <- ggplot() +
    geom_path(data = tick_circle_df, aes(x = x, y = y, group = group), colour = "grey85", linetype = "dashed", linewidth = 0.38) +
    geom_segment(data = spoke_df, aes(x = x, y = y, xend = xend, yend = yend), colour = "grey90", linewidth = 0.30) +
    geom_segment(data = neg_axis_df, aes(x = x, y = y, xend = xend, yend = yend), colour = "black", linewidth = 0.62) +
    geom_segment(data = neg_tick_seg_df, aes(x = x, y = y, xend = xend, yend = yend), colour = "black", linewidth = 0.50) +
    geom_text(data = neg_label_df, aes(x = x, y = y, label = label), family = base_family, size = 3.7, fontface = "bold", vjust = 1) +
    geom_text(data = neg_axis_title_df, aes(x = x, y = y, label = label), family = base_family, size = 3.8, fontface = "bold") +
    geom_segment(data = pos_axis_df, aes(x = x, y = y, xend = xend, yend = yend), colour = "black", linewidth = 0.62) +
    geom_segment(data = pos_tick_seg_df, aes(x = x, y = y, xend = xend, yend = yend), colour = "black", linewidth = 0.50) +
    geom_text(data = pos_label_df, aes(x = x, y = y, label = label), family = base_family, size = 3.7, fontface = "bold", vjust = 1) +
    geom_text(data = pos_axis_title_df, aes(x = x, y = y, label = label), family = base_family, size = 3.8, fontface = "bold") +
    geom_polygon(data = bar_df, aes(x = x, y = y, group = group), fill = "#D7191C", alpha = 0.94, colour = NA) +
    geom_path(data = zero_circle_df, aes(x = x, y = y), colour = "black", linewidth = 0.65) +
    geom_point(data = shap_polar_df, aes(x = x, y = y, colour = feature_value_norm), size = 2.35, alpha = 0.74) +
    geom_path(data = outer_circle_df, aes(x = x, y = y), colour = "black", linewidth = 0.65) +
    geom_segment(data = outer_tick_df, aes(x = x, y = y, xend = xend, yend = yend), colour = "black", linewidth = 0.48) +
    geom_segment(data = mean_axis_df, aes(x = x, y = y, xend = xend, yend = yend), colour = "black", linewidth = 0.55) +
    geom_segment(data = mean_tick_seg_df, aes(x = x, y = y, xend = xend, yend = yend), colour = "black", linewidth = 0.48) +
    geom_text(data = mean_label_df, aes(x = x, y = y, label = label, hjust = hjust), family = base_family, size = 3.7, fontface = "bold") +
    annotate("text", x = 0, y = bar_scale + mean_tick_len * 4.6, label = "Mean |SHAP|", family = base_family, fontface = "bold", size = 4.6) +
    geom_text(data = label_df, aes(x = x, y = y, label = label, hjust = hjust), family = base_family, size = 4.5, fontface = "bold", lineheight = 0.92) +
    scale_colour_gradientn(colours = c("#3B4CC0", "#8DB0FE", "#F7F7F7", "#F4987A", "#B40426"),
                           values = scales::rescale(c(0, 0.25, 0.5, 0.75, 1)), limits = c(0, 1),
                           breaks = c(0, 1), labels = c("Low", "High"), name = "Feature value") +
    annotate("text", x = -label_radius * 0.74, y = -label_radius * 1.03, label = "Red rose = Mean |SHAP|",
             colour = "#D7191C", family = base_family, fontface = "bold", size = 4.4) +
    annotate("text", x = 0, y = zero_radius + base_span * 0.06, label = "SHAP = 0",
             family = base_family, fontface = "bold", size = 4.0) +
    coord_equal(xlim = c(-label_radius * 1.12, label_radius * 1.12),
                ylim = c(-label_radius * 1.12, label_radius * 1.12), clip = "off") +
    theme_void(base_family = base_family) +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 11, hjust = 0.5),
          plot.margin = ggplot2::margin(20, 45, 55, 45)) +
    labs(title = paste("SHAP Polar Beeswarm Plot —", nm),
         subtitle = "Center = global importance; outer ring = local SHAP effects on predicted Treat probability") +
    guides(colour = guide_colourbar(title.position = "top", title.hjust = 0.5,
                                    barwidth = grid::unit(10, "cm"), barheight = grid::unit(0.45, "cm")))
  pdf(file.path(out, paste0("SHAP_Polar_", nm, ".pdf")), width = 10, height = 10)
  print(p_polar); dev.off()
  cat("saved SHAP_Polar_", nm, ".pdf\n", sep = "")

  # ---- beeswarm (per model) ----
  feature_y_lookup <- tibble(feature = feature_levels, y_base = seq(n_feat, 1, by = -1))
  shap_bee_df <- shap_long %>% mutate(feature = as.character(feature)) %>%
    left_join(feature_y_lookup, by = "feature") %>%
    group_by(feature) %>%
    mutate(y_offset = simple_beeswarm(shap, nbins = 80, width = 0.72), y_plot = y_base + y_offset) %>% ungroup()
  p_beeswarm <- ggplot(shap_bee_df, aes(x = shap, y = y_plot, colour = feature_value_norm)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, colour = "grey45") +
    geom_point(size = 2.1, alpha = 0.78) +
    scale_y_continuous(breaks = seq_len(n_feat), labels = rev(feature_levels),
                       expand = expansion(add = c(0.6, 0.6))) +
    scale_colour_gradientn(colours = c("#3B4CC0", "#8DB0FE", "#F7F7F7", "#F4987A", "#B40426"),
                           values = scales::rescale(c(0, 0.25, 0.5, 0.75, 1)), limits = c(0, 1),
                           breaks = c(0, 1), labels = c("Low", "High"), name = "Feature value") +
    labs(title = paste("SHAP Beeswarm Plot —", nm),
         subtitle = "SHAP > 0 increases predicted Treat Prob; SHAP < 0 decreases predicted Treat Prob",
         x = "SHAP value for Treat probability prediction", y = NULL) +
    theme_classic(base_family = base_family) +
    theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 11, hjust = 0.5),
          axis.text.y = element_text(size = 12, face = "bold"), axis.text.x = element_text(size = 11),
          axis.title.x = element_text(size = 13, face = "bold"), legend.position = "right",
          legend.title = element_text(size = 11, face = "bold"), legend.text = element_text(size = 10),
          panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.35),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank()) +
    guides(colour = guide_colourbar(title.position = "top", title.hjust = 0.5,
                                    barheight = grid::unit(4.8, "cm")))
  pdf(file.path(out, paste0("SHAP_Beeswarm_", nm, ".pdf")), width = 10, height = max(6, n_feat * 0.45))
  print(p_beeswarm); dev.off()
  cat("saved SHAP_Beeswarm_", nm, ".pdf\n", sep = "")

  # ---- bar (per model) ----
  bar_plot_df <- feature_importance %>%
    mutate(feature = factor(as.character(feature), levels = rev(feature_levels)))
  p_bar <- ggplot(bar_plot_df, aes(x = feature, y = mean_abs_shap)) +
    geom_col(fill = "#D7191C", width = 0.72) + coord_flip() +
    labs(title = paste("SHAP Feature Importance —", nm),
         subtitle = "Global importance measured by mean absolute SHAP value", x = NULL,
         y = "Mean |SHAP value|") +
    theme_classic(base_family = base_family) +
    theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 11, hjust = 0.5),
          axis.text.y = element_text(size = 12, face = "bold"), axis.text.x = element_text(size = 11),
          axis.title.x = element_text(size = 13, face = "bold"),
          panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.35),
          panel.grid.minor = element_blank())
  pdf(file.path(out, paste0("SHAP_Bar_", nm, ".pdf")), width = 8, height = max(5, n_feat * 0.4))
  print(p_bar); dev.off()
  cat("saved SHAP_Bar_", nm, ".pdf\n", sep = "")
}
cat("\nALL DONE\n")
