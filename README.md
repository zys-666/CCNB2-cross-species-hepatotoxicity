# CCNB2-cross-species-hepatotoxicity

Reproducible R scripts for: **CCNB2 as a cross-species hepatotoxicant-responsive candidate gene**
(crucian carp NTH model → rat DrugMatrix GSE57815 → mouse GSE44783).

All results reported in the accompanying manuscript are produced by the scripts listed under
**Main pipeline** below, run in the listed order with R 4.6.0. Software versions are written to
`sessionInfo.txt` by the analysis scripts.

> **Working directory.** Every script assumes the working directory contains `inputs/` and `outputs/`
> (the `setwd()` line at the top of each script). Adjust that line — or `PROJ` in `driver_rerun.R` —
> to your own layout. Raw expression matrices are not distributed here (see *Data*).

## Cohorts & roles
| Cohort | Accession / source | Role |
|---|---|---|
| Crucian carp (NTH experiment) | own RNA-seq | **Discovery** (source of candidate genes) |
| Rat liver, DrugMatrix | **GSE57815** (GPL1355) | Screening: 9-classifier consensus ranking |
| Mouse liver, 13 model compounds | **GSE44783** (GPL1261), PMID 24040119 | **Independent external consistency** cohort |

The crucian carp transcriptome is the discovery cohort and is **not** an independent validation cohort.

## Label rule (rat)
Labels verified sample-by-sample against the GEO dose records: **numeric dose == 0 mg/kg → vehicle control
(n = 279); dose > 0 → exposed (n = 1,939, including low-dose arms of 0.01–0.65 mg/kg)**. Compound identity is
not present in the GEO sample metadata; it was taken from the **CEBS DrugMatrix annotation (CEBS accession
15670_23401)** by joining on the CEL barcode (`ARRAY_ID = <barcode>.CEL`), which matched **2,218/2,218
samples (100%)** and yielded **200 named compounds** plus 279 vehicle-only controls.

## Main pipeline (run in order; `driver_rerun.R` runs all steps)

| # | Script | What it does |
|---|---|---|
| 1 | `01_prep.R` | Build the labelled expression object (`rat_expr45.rds`) + label audit |
| 2 | `05_verify.R` | Alignment verification (explicit `match()` + `stopifnot()`, never `merge()`) and positive-control genes |
| 3 | `06_ccnb2.R` | Single-gene CCNB2, all rat strata (L1/L3, time point, vehicle) and mouse |
| 4 | `09_fish_summary.R` | Three-species CCNB2 summary table + base-R forest plot |
| 5 | **`training_v3.R`** | **9 classifiers** (RF/SVM/XGB/GLM/GBM/KNN/NNET/LASSO/DT): stratified 5×3 CV, permutation importance + ranks, SHAP, ROC, **calibration + ECE**, figures, `models_for_figures.rds` (~40 min) |
| 6 | **`nested_cv_v3.R`** | **Compound-grouped nested CV** + compound-clustered bootstrap confidence intervals |
| 7 | `youden_ccnb2.R` | Youden-optimal sensitivity/specificity (rat & mouse) |
| 8 | `loco_rerun.R`, `loco_dot.R` | Leave-one-compound-out analysis (mouse) + dot plot |
| 9 | `fig_importance_bars.R`, `fig_shap_polar.R` | Importance figure and SHAP polar/beeswarm/bar panels per model |
| 10 | **`fig_grouped_cv.R`** | **Figure 3C** (three-caliber CV comparison) and the per-fold stability figure |
| 11 | **`fig3H.R`** | **Figure 3H**: consistency of CCNB2 induction across rat strata (duration / vehicle) |
| 12 | **`fig_style2.R`** | Journal-style vector figures: **Figure 3K** forest plot (three cohorts), CV-performance and fish dose–response panels |
| 13 | `s5_sensitivity.R` | Mouse sensitivity strata (time × sex) + per-compound table |
| 14 | `appendix_tables.R`, `format_supp_tables.R` | Supplementary tables (per-sample metadata, classifier performance, ranking matrix, LOCO/Youden) |
| 15 | **`ortholog_map.R`** | Cross-species orthologue mapping (NCBI HomoloGene, offline version-pinned database) |
| 16 | `audit_all_groups.R` | Independent label/group audit of all input files |
| — | `test_caret_foldwise.R`, `verify_resample.R` | Verification checks (see *Verification scripts*); not part of the reported pipeline |

## Verification scripts

Two scripts document checks that are referred to in the Methods and in the response to review.
They do not produce any reported number or figure.

| Script | Check |
|---|---|
| `test_caret_foldwise.R` | Confirms that caret's `preProcess` is fitted **within each resampling training partition** (toy data, `preProcess = "range"`, comparison of the two scalings on the first held-out fold) |
| `verify_resample.R` | Confirms that caret's `$resample` reflects **only the final tuning configuration** (compares `$resample$ROC` with the AUC recomputed from `$pred` after filtering to `bestTune`) |
| `05_verify.R` | Alignment verification and positive-control genes (see Main pipeline) |
| `audit_all_groups.R` | Independent label/group audit of all input files (see Main pipeline) |

## Software environment

`sessionInfo()` output for every analysis step is in [`sessionInfo/`](sessionInfo/) — one file per
script group, with a summary table of the key package versions (R 4.6.0; caret 7.0-1, xgboost 3.2.1.1,
randomForest 4.7-1.2, gbm 2.3.1, kernlab 0.9-33, glmnet 5.0, nnet 7.3-20, rpart 4.1.27, pROC 1.19.0.1,
DALEX 2.5.3, homologene 1.4.68.19.3.27). Versions for components outside the R environment
(read quantification, enrichment, docking, MD) are given in the manuscript and its Supplementary Material.

## Methods summary (full details in the manuscript)

- **Model tuning.** Hyperparameters were **selected by cross-validation within the training partitions**,
  not fixed a priori: during the nested cross-validation by an inner 5-fold search inside each outer training
  partition, and otherwise by the resampling grid built into each classifier's training function. No outer
  test-fold data were used for tuning.
- **Fold construction.** Stratified sample-level folds for the 5×3 CV; **compound-grouped folds** for the
  nested CV — outer and inner folds are formed over groups (200 compound groups + 58 dosing-batch groups =
  258), so no compound is split across training and test. The scripts assert that the number of groups
  spanning more than one fold is zero.
- **Prediction footing.** Only the final hyperparameter configuration is retained, and the repeated
  predictions of each sample are averaged, so that **every sample contributes exactly one out-of-fold
  prediction**; all AUCs, intervals and calibration errors are computed on that footing.
  Pooled out-of-fold AUC **0.959 (DT) – 0.995 (XGBoost)**; expected calibration error **0.002–0.019**.
- **Internal discrimination.** Stratified 5×3 CV **0.949–0.994**; fully nested CV **0.948–0.993**
  (sample-level) and **0.948–0.995** (compound-grouped), i.e. grouping does not change the conclusions
  (each classifier within 0.006 of its sample-level per-fold mean).
- **Confidence intervals.** Because samples sharing a compound are not mutually independent,
  **compound-clustered bootstrap** intervals (2,000 resamples over the grouping variable) are reported for
  both the single-gene estimates and the pooled out-of-fold AUCs of the nine classifiers.
- **Consensus ranking.** Permutation importance (training side) + univariate AUC per gene; **CCNB2 ranked
  first in all nine classifiers and in univariate AUC (mean rank 1.0/45)**.
- **External consistency.** CCNB2 up-regulated in mouse (log₂FC = +0.944, AUC 0.741), 10/13 compounds
  induced, LOCO mean AUC 0.771 (single gene). The mouse cohort was never used for ranking or tuning.
  The LOCO analysis refits a model **within** the mouse cohort and therefore tests the directional
  single-gene signal across compounds; it is **not** transfer validation of the nine rat-trained
  classifiers, and no such cross-species transfer claim is made.
- **Winner's-curse caveat.** The rat single-gene AUC is descriptive (the rat cohort is also the ranking
  cohort); independent confirmation comes from the pre-specified mouse analysis.
- **Gene-symbol mapping.** Orthologue correspondence verified with NCBI HomoloGene; 39 of the 45 panel genes
  map 1:1:1:1 across rat, mouse, zebrafish and human (CCNB2 included).
- **Batch correction.** Expression matrices were processed **within each cohort** (between-array
  normalization, log₂ transformation, probe-to-gene annotation). **No batch-effect correction (e.g. ComBat)
  was applied**, because no pooled cross-platform matrix and no multi-gene cross-platform model was
  constructed. No held-out data informed any preprocessing, feature-selection, tuning or ranking step.
- **Random seeds** are recorded in the scripts at every stochastic step (main pipeline 123; XGBoost folds
  200 + r; outer grouped folds 300 + r; inner grouped folds 11; native XGBoost subsampling 1; nine-classifier
  cluster bootstrap 99; single-gene cluster bootstrap 7).

## Related data in this repository
`RNAseq_supplementary/` holds the NTH-exposure RNA-seq deliverable for *Carassius auratus*
(crucian carp / goldfish), prepared for a reviewer response: QC metrics, mapping statistics,
reference genome metadata, gene-ID conversion evidence, expression matrices, DEG/enrichment
outputs and expression-level QC figures. Note that reads were aligned to the *Carassius gibelio*
reference genome (carGib1.2-hapl.c, GCF_023724105.1) — see the notes in
`RNAseq_supplementary/README.md`, which also documents which DEG output the manuscript reports
and gives a full file index.

## Data
Raw expression matrices are large and must be obtained from the original sources — `normalize.txt`
(GSE57815 samples, from the DrugMatrix/GEO distribution) and `mouse_expr.csv` (GSE44783). Sample-level
metadata used for label definition are derived from the GEO series GSM records; the rat compound identity
table comes from CEBS accession 15670_23401. The crucian carp expression values reported in the manuscript
are provided in `RNAseq_supplementary/` (see above).
