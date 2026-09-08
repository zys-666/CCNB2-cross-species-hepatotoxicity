# appendix_tables.R : generate Supplementary metadata tables S1 (GSE57815) & S2 (GSE44783)
suppressMessages({library(data.table)})
setwd("D:/MLroute2")
out <- "outputs/paper/appendices"; dir.create(out, showWarnings = FALSE, recursive = TRUE)

## ---- GSE57815 (rat) ----
r <- fread("inputs/rat_samples.tsv")
r[, dnum := suppressWarnings(as.numeric(sub(" .*", "", trimws(dose))))]
r[, dose0 := !is.na(dnum) & dnum == 0]
r[, dnum := NULL]
r[, lab := ifelse(dose0, "Control", "Treat")]
sink(file.path(out, "Supplementary_Table_S1_GSE57815.txt"))
cat("Supplementary Table S1. GSE57815 (rat DrugMatrix liver) sample metadata\n")
cat("Series title      : Exposure of rat to a variety of toxicants, liver assayed by Affymetrix microarray\n")
cat("Platform          : GPL1355 [Rat230_2] Affymetrix Rat Genome 230 2.0 Array\n")
cat("Organism/Tissue   : Rattus norvegicus, liver; all male Sprague-Dawley\n")
cat("Preprocessing     : RMA log2 (GeneSpring GX 11.0.5.1, per-sample GEO VALUE)\n")
cat("Total samples     :", nrow(r), "\n")
cat("Labels (this study, audited L1 rule): Control (dose = 0)", sum(r$lab == "Control"),
    "; Exposed (dose > 0)", sum(r$lab == "Treat"), "\n")
cat("Vehicles (all; n):\n"); print(table(r$vehicle))
cat("Vehicles (dose=0 controls only):\n"); print(table(r$vehicle[r$dose0]))
cat("Time points (n):\n"); print(table(r$time))
cat("Routes (n):\n"); print(table(r$route))
cat("NOTE: compound identities are not stored in GEO sample records; the sample-to-compound\n")
cat("map is available from CEBS accession 004-00008-0000-000-2 (approx. 200 compounds per the\n")
cat("original DrugMatrix liver series). Full per-sample table: rat_samples_full.tsv\n")
sink()

## ---- GSE44783 (mouse) ----
m <- fread("inputs/mouse_samples.tsv")
veh <- c("Corn oil", "Carboxymethyl cellulose", "Water")
m[, grp := ifelse(treatment %in% veh, "Vehicle", "Compound")]
sink(file.path(out, "Supplementary_Table_S2_GSE44783.txt"))
cat("Supplementary Table S2. GSE44783 (CD-1 mouse liver, Novartis toxicogenomics) sample metadata\n")
cat("Series title      : Expression data from CD-1 mouse liver samples ... genotoxic/non-genotoxic carcinogens or non-hepatocarcinogens\n")
cat("Reference         : PMID 24040119 (Eichner et al.)\n")
cat("Platform          : GPL1261 [Mouse430_2] Affymetrix Mouse Genome 430 2.0 Array\n")
cat("Organism/Tissue   : Mus musculus (CD-1), liver; male and female\n")
cat("Preprocessing     : RMA log2 (R/Bioconductor)\n")
cat("Total samples     :", nrow(m), "\n")
cat("Groups: vehicle control", sum(m$grp == "Vehicle"), "; compound-exposed", sum(m$grp == "Compound"), "\n")
cat("Vehicle controls (n):\n"); print(table(m$treatment[m$grp == "Vehicle"]))
cat("Compounds (n; sex split M/F; time 4 d / 15 d):\n")
tab <- m[grp == "Compound", .(n = .N, M = sum(gender == "male"), F = sum(gender == "female"),
                              t4 = sum(time == "4 days"), t15 = sum(time == "15 days")), by = treatment]
print(tab[order(-n)])
cat("NOTE: per-compound doses were selected per CPDB carcinogenic doses (see PMID 24040119);\n")
cat("doses are not recorded per-sample in GEO. Full per-sample table: mouse_samples_full.tsv\n")
sink()

file.copy("inputs/rat_samples.tsv", file.path(out, "rat_samples_full.tsv"), overwrite = TRUE)
file.copy("inputs/mouse_samples.tsv", file.path(out, "mouse_samples_full.tsv"), overwrite = TRUE)
cat("appendices written to", normalizePath(out), "\n")
