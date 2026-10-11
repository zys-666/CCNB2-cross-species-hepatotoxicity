############################################################
## ortholog_map.R — 直系同源基因映射工作流(离线可复现)
##
## 目的:对候选基因面板做正式的跨物种直系同源映射,给出
##       ① 每个基因在 大鼠/小鼠/斑马鱼/人 中的对应符号
##       ② 映射类型(1:1:1:1 / 一对多 / 未映射 / 大鼠内歧义)
##       ③ 预指定基因 CCNB2 的逐物种核验
##
## 数据源(默认 A):
##   A. homologene 本地库(NCBI HomoloGene;随包版本锁定,完全离线)
##   B. biomaRt 在线查询(必须锁定 Ensembl archive 版本,否则不可复现;见文末注释块)
##
## 输入 : inputs/interGenes.txt   (53 个候选基因,上游注释符号)
##        outputs/rat_expr45.rds  (实际进入建模的 45 个基因)
## 输出 : outputs/ortholog/ortholog_map.tsv
##        outputs/ortholog/ortholog_summary.txt
##        outputs/ortholog/sessionInfo.txt
############################################################
suppressMessages({library(data.table); library(homologene)})
setwd("D:/MLroute2")
out <- "outputs/ortholog"; dir.create(out, showWarnings = FALSE, recursive = TRUE)

TAX <- c(rat = 10116, mouse = 10090, zebrafish = 7955, human = 9606)

## ---------- 1. 本地同源库 ----------
e <- new.env(); data("homologeneData", package = "homologene", envir = e)
hg <- as.data.table(get("homologeneData", envir = e))
cat("homologene 库: rows =", nrow(hg), "| 原始列名 =", paste(colnames(hg), collapse = ","), "\n")
## 实际列名是 HID / Gene.ID / Gene.Symbol / Taxonomy —— 不要凭记忆重命名
stopifnot(all(c("HID", "Gene.Symbol", "Taxonomy") %in% colnames(hg)))
hg[, Gene.Symbol := toupper(trimws(Gene.Symbol))]
hg[, Taxonomy := as.integer(Taxonomy)]
cat("库内物种数:", length(unique(hg$Taxonomy)),
    "| 目标物种均在库中:", all(TAX %in% hg$Taxonomy), "\n")
cat("homologene 包版本:", as.character(packageVersion("homologene")), "\n\n")

## ---------- 2. 待映射基因 ----------
genes53 <- unique(as.vector(read.table("inputs/interGenes.txt", header = FALSE,
                                       sep = "\t", stringsAsFactors = FALSE)[, 1]))
genes45 <- rownames(readRDS("outputs/rat_expr45.rds")$X)
cat("候选清单 interGenes.txt:", length(genes53),
    "| 实际建模面板:", length(genes45), "\n\n")

## ---------- 3. 逐基因映射 ----------
map_one <- function(sym) {
  s <- toupper(trimws(sym))
  hids <- unique(hg[Taxonomy == TAX[["rat"]] & Gene.Symbol == s, HID])
  got <- lapply(TAX, function(tx) {
    if (!length(hids)) return(character(0))
    sort(unique(hg[Taxonomy == tx & HID %in% hids, Gene.Symbol]))
  })
  names(got) <- names(TAX)
  n <- vapply(got, length, integer(1))
  cls <- if (n[["rat"]] == 0) "unmapped_in_homologene"
         else if (length(hids) > 1) "ambiguous_in_rat"
         else if (all(n[c("mouse", "zebrafish", "human")] == 1)) "1:1:1:1"
         else if (any(n[c("mouse", "zebrafish", "human")] > 1)) "one_to_many"
         else "partial"
  data.frame(symbol_input = sym,
             symbol_rat = paste(got$rat, collapse = ";"),
             n_rat_hid = length(hids),
             symbol_mouse = paste(got$mouse, collapse = ";"),
             symbol_zebrafish = paste(got$zebrafish, collapse = ";"),
             symbol_human = paste(got$human, collapse = ";"),
             n_rat = n[["rat"]], n_mouse = n[["mouse"]],
             n_zebrafish = n[["zebrafish"]], n_human = n[["human"]],
             mapping_class = cls, stringsAsFactors = FALSE)
}
tab <- rbindlist(lapply(genes53, map_one))
tab[, in_model_panel := symbol_input %in% genes45]
fwrite(tab, file.path(out, "ortholog_map.tsv"), sep = "\t")

## ---------- 4. 汇总 ----------
sink(file.path(out, "ortholog_summary.txt"))
cat("=== 直系同源映射汇总(homologene 本地库, 包版本",
    as.character(packageVersion("homologene")), ")===\n")
cat("输入基因:", nrow(tab), "| 进入建模面板:", sum(tab$in_model_panel), "\n\n")
cat("-- 映射类型分布(全部 53)--\n");  print(table(tab$mapping_class))
cat("\n-- 映射类型分布(建模面板 45)--\n"); print(table(tab[in_model_panel == TRUE, mapping_class]))
cat("\n-- 每物种唯一映射数 --\n")
for (tx in c("rat", "mouse", "zebrafish", "human"))
  cat(sprintf("  %-10s %d / %d\n", tx, sum(tab[[paste0("n_", tx)]] == 1), nrow(tab)))
cat("\n-- 预指定基因 CCNB2 核验 --\n")
print(as.data.frame(tab[toupper(symbol_input) == "CCNB2",
      .(symbol_input, symbol_rat, symbol_mouse, symbol_zebrafish, symbol_human, mapping_class)]),
      row.names = FALSE)
cat("\n-- 建模面板中非 1:1:1:1 的基因(需人工确认)--\n")
print(as.data.frame(tab[in_model_panel == TRUE & mapping_class != "1:1:1:1",
      .(symbol_input, mapping_class, symbol_mouse, symbol_zebrafish)]), row.names = FALSE)
sink()

## ---------- 5. 可选:biomaRt 在线路径(须锁定 archive 才可复现)----------
## library(biomaRt)
## ENSEMBL_ARCHIVE <- "https://nov2020.archive.ensembl.org"   # 必须写死版本
## rn <- useEnsembl("ensembl", dataset = "rnorvegicus_gene_ensembl", host = ENSEMBL_ARCHIVE)
## mm <- useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl",  host = ENSEMBL_ARCHIVE)
## bm <- getLDS(attributes = c("rgd_symbol"), filters = "rgd_symbol", values = genes45,
##              mart = rn, attributesL = c("mgi_symbol"), martL = mm)
## fwrite(bm, file.path(out, "ortholog_map_biomart.tsv"), sep = "\t")

sink(file.path(out, "sessionInfo.txt")); print(sessionInfo()); sink()
cat("ALL DONE ->", out, "\n")
