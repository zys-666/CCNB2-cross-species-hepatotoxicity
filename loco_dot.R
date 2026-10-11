# loco_dot.R : 13-compound LOCO dot plot (single-gene CCNB2), colored by induction status
suppressMessages({library(ggplot2); library(data.table)})
setwd("D:/MLroute2")
l <- fread("outputs/paper/GSE44783_LOCO.txt")
setnames(l, c("compound", "n_trt", "AUC_CCNB2", "AUC_g3"))
# induction 状态（来自结果文件：10/13 显著上调；3 个非诱导 = Cefuroxime/Propranolol/CPA）
no_ind <- c("Cefuroxime", "Propranolol", "Cyproterone acetate")
l[, induced := ifelse(compound %in% no_ind, "No CCNB2 induction (n=3)", "CCNB2 induced (n=10)")]
l[, induced := factor(induced, levels = c("CCNB2 induced (n=10)", "No CCNB2 induction (n=3)"))]
l <- l[order(AUC_CCNB2)]
l[, compound := factor(compound, levels = compound)]

p <- ggplot(l, aes(x = AUC_CCNB2, y = compound, colour = induced)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey45", linewidth = 0.5) +
  geom_segment(aes(xend = 0.5, yend = compound), linewidth = 0.6, alpha = 0.5) +
  geom_point(size = 4) +
  geom_text(aes(label = sprintf("%.2f", AUC_CCNB2)), hjust = -0.5, size = 3.4, colour = "black") +
  scale_colour_manual(values = c("CCNB2 induced (n=10)" = "#D7191C", "No CCNB2 induction (n=3)" = "grey55")) +
  scale_x_continuous(limits = c(0.5, 1.03), breaks = seq(0.5, 1, 0.1)) +
  labs(title = "Leave-one-compound-out AUC for CCNB2 (mouse GSE44783)",
       subtitle = "Model trained on all other compounds; test = held-out compound + half of controls (5 repeats, mean)",
       x = "LOCO AUC (CCNB2 single gene, logistic)", y = NULL, colour = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
        plot.subtitle = element_text(hjust = 0.5, size = 9.5))
pdf("outputs/paper/GSE44783_LOCO_dot.pdf", width = 8, height = 6)
print(p); dev.off()
cat("saved outputs/paper/GSE44783_LOCO_dot.pdf | size:", file.info("outputs/paper/GSE44783_LOCO_dot.pdf")$size, "\n")
