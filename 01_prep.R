# 01_prep.R : label audit + expression subsetting (Route-2 stage 1)
suppressMessages({library(data.table)})
setwd("D:/MLroute2")

# ---------- 1. header audit: what are the 279/1939 labels really? ----------
con <- file("inputs/normalize.txt", "rt")
header <- readLines(con, n = 1); close(con)
cols <- strsplit(header, "\t")[[1]][-1]
gsm <- sub("_(Control|Treat)$", "", cols)
flab <- sub("^.*_(Control|Treat)$", "\\1", cols)
cat("n columns:", length(cols), " unique GSM:", length(unique(gsm)), "\n")
stopifnot(all(grepl("^GSM", gsm)))

meta <- fread("inputs/rat_samples.tsv")
mi <- match(gsm, meta$acc)                      # align metadata to file column order directly
stopifnot(!anyNA(mi))
meta <- meta[mi]
stopifnot(identical(gsm, meta$acc))
m <- data.table(gsm = gsm, file_label = flab, dose = meta$dose, vehicle = meta$vehicle,
                time = meta$time, route = meta$route)
cat("unmatched GSM:", 0L, "\n")
m[, dnum := suppressWarnings(as.numeric(sub(" .*", "", trimws(dose))))]
m[, dose0 := !is.na(dnum) & dnum == 0]     # 数值剂量==0 才是对照（不能用 grepl("^0")：会吞 0.5 mg/kg 等低剂量）
m[, dnum := NULL]

cat("\n--- file label x dose0 ---\n"); print(table(m$file_label, m$dose0, useNA = "ifany"))
cat("\n--- file label x vehicle (rows=dose0 cols=label? ) ---\n")
print(table(vehicle = m$vehicle, file_label = m$file_label, dose0 = m$dose0))
cat("\nTreat samples that are dose==0 (contamination candidates):\n")
print(m[dose0 == TRUE & file_label == "Treat", .N, by = .(vehicle, time)])
cat("Control samples that are dose>0 (if any):\n")
print(m[dose0 == FALSE & file_label == "Control", .N, by = .(vehicle)])

# ---------- 2. candidate corrected labels ----------
m[, L1 := ifelse(dose0, "Control", "Treat")]                       # all dose-0 = control
m[, keep3 := vehicle %in% c("Corn Oil", "CMC", "Water")]           # 3-vehicle design (manuscript story)
m[, L3 := ifelse(!keep3, NA, ifelse(dose0, "Control", "Treat"))]   # subset to 3 vehicles
cat("\nL1 (all dose0): Control", sum(m$L1 == "Control"), " Treat", sum(m$L1 == "Treat"), "\n")
cat("L3 (3 vehicles only): Control", sum(m$L3 == "Control", na.rm = TRUE),
    " Treat", sum(m$L3 == "Treat", na.rm = TRUE), " dropped", sum(is.na(m$L3)), "\n")

fwrite(m[, .(gsm, file_label, dose0, vehicle, time, route, L1, keep3, L3)], "outputs/label_audit.tsv", sep = "\t")

# ---------- 3. read full expression once, subset to interGenes ----------
ig <- toupper(trimws(readLines("inputs/interGenes.txt")))
cat("\ninterGenes n:", length(ig), "\n")
expr <- fread("inputs/normalize.txt", sep = "\t", header = TRUE, check.names = FALSE)
gn <- expr[[1]]
idx <- which(toupper(gn) %in% ig)
cat("interGenes matched in rat matrix:", length(idx), "\n")
genes_matched <- gn[idx]
X <- as.matrix(expr[idx, -1, with = FALSE]); rownames(X) <- genes_matched
saveRDS(list(X = X, gsm = gsm, flab = flab, m = m), "outputs/rat_expr45.rds")

three <- toupper(c("Ccnb2", "Plk1", "Ube2c"))
idx3 <- which(toupper(genes_matched) %in% three)
cat("3 genes present:", genes_matched[idx3], "\n")
stopifnot(length(idx3) == 3)
saveRDS(X[idx3, , drop = FALSE], "outputs/rat_expr3.rds")
cat("DONE prep\n")
