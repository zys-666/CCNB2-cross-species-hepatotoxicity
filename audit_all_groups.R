# audit_all_groups.R : full label/group audit — rat & mouse, all grouping files
suppressMessages({library(data.table)})
setwd("D:/MLroute2")
cat("################ RAT GSE57815 ################\n")
meta <- fread("inputs/rat_samples.tsv")
meta[, dose := trimws(dose)]
meta[, num := suppressWarnings(as.numeric(sub(" .*", "", dose)))]
meta[, dose0 := !is.na(num) & num == 0]

# normalize.txt header -> file labels
hdr <- fread("inputs/normalize.txt", nrows = 0)
cn  <- names(hdr)[-1]                    # 第1列为基因名行名
flab <- sub(".*_", "", cn)
cn0 <- sub("_(Control|Treat)$", "", cn)  # 纯 GSM id
cat("normalize ncol:", length(cn), "| rat_samples nrow:", nrow(meta), "\n")
cat("GSM set identical:", setequal(cn0, meta$acc),
    "| same order:", identical(cn0, meta$acc), "\n")
if (!identical(cn0, meta$acc)) {
  cat("first order diff at:", which(cn0 != meta$acc)[1],
      "| in meta but not normalize:", setdiff(meta$acc, cn0)[1:5],
      "| in normalize but not meta:", setdiff(cn0, meta$acc)[1:5], "\n")
}
# 关键：按 normalize 列序对齐 meta 后再挂 file_label
meta2 <- meta[match(cn0, meta$acc)]
stopifnot(!anyNA(meta2$acc), identical(meta2$acc, cn0))
meta2[, file_label := flab]
meta <- meta2
rm(meta2); gc()
cat("\n--- file_label counts (normalize suffix) ---\n"); print(table(meta$file_label))
cat("\n--- TRUE dose0 counts ---\n"); print(table(meta$dose0))
cat("\n--- cross-tab file_label x dose0(TRUE rule) ---\n"); print(table(flab = meta$file_label, dose0_true = meta$dose0))
cat("\n--- file_label='Control' samples whose dose text != '0 mg/kg' ---\n")
print(meta[file_label == "Control" & dose != "0 mg/kg", .N, by = dose])
cat("\n--- file_label='Treat' samples with dose NA/missing/empty ---\n")
print(meta[file_label == "Treat" & (is.na(dose) | dose == "" | is.na(num)), .N, by = .(dose, vehicle)])
cat("\n--- file_label='Treat' but dose text exactly '0 mg/kg' ---\n")
print(meta[file_label == "Treat" & dose == "0 mg/kg", .N, by = .(vehicle, time)])
cat("\n--- any duplicated GSM in meta? ---\n", anyDuplicated(meta$acc), "\n")

# labels_L1.txt 与修正规则比对
l1 <- fread("outputs/v2/labels_L1.txt")
setnames(l1, "dose0", "dose0_l1")
cat("\n--- labels_L1.txt: rows:", nrow(l1), "| counts ---\n"); print(table(l1$Type))
mrg <- merge(meta[, .(acc, dose0, file_label)], l1[, .(GSM, Type, dose0_l1)],
             by.x = "acc", by.y = "GSM")
cat("--- labels_L1 误标为 Control 但真实 dose>0 的样本数（^0 bug 受害者）---\n")
cat(mrg[Type == "Control" & dose0 == FALSE, .N], "\n")
cat("labels_L1 Control:", sum(mrg$Type == "Control"), "| TRUE dose0 Control:", sum(mrg$dose0), "\n")

cat("\n################ MOUSE GSE44783 ################\n")
grp <- fread("inputs/mouse_group.csv"); setnames(grp, c("SampleID", "group"))
st  <- fread("inputs/mouse_samples.tsv")
cat("group file rows:", nrow(grp), "| counts:\n"); print(table(grp$group))
cat("mouse_samples rows:", nrow(st), "| unique acc:", uniqueN(st$acc), "\n")
cat("group vs mouse_samples acc: same set:", setequal(grp$SampleID, st$acc), "\n")
cat("mouse_samples acc order == group order:", identical(st$acc, grp$SampleID), "\n")
# expression colnames
mg <- fread("inputs/mouse_expr.csv", sep = ",", nrows = 0)
mcn <- names(mg)[-1]
cat("mouse_expr cols:", length(mcn), "| == group order:", identical(mcn, grp$SampleID),
    "| == mouse_samples order:", identical(mcn, st$acc), "\n")
# 对照组的元数据核查: control 不应有化合物/剂量
st2 <- copy(st); st2[, grp2 := grp$group[match(acc, grp$SampleID)]]
cat("\n--- group x (treatment is NA?) ---\n")
print(st2[, .N, by = .(grp2, has_trt = !is.na(treatment) & treatment != "")])
cat("\n--- 对照组 treatment/剂量/来源抽查 ---\n")
print(st2[grp2 == "Control", .N, by = .(treatment, source, time)])
cat("\n--- 对照组 vehicle 词汇(从 source/title 提取前 2 词) ---\n")
print(st2[grp2 == "Control", .N, by = sub(" mice.*|,.*", "", source)][order(-N)])
cat("\n--- 处理组 compound 计数 ---\n")
print(st2[grp2 == "Treat", .N, by = treatment][order(-N)])
cat("\n--- 组别 x 性别 x 时间 ---\n")
print(st2[, .N, by = .(grp2, gender, time)])
cat("\n--- 450 组别与样本文件无缺失 ---\n")
cat("grp unmatched in mouse_samples:", sum(!grp$SampleID %in% st$acc), "\n")
