# CCNB2-cross-species-hepatotoxicity

Reproducible R scripts for: **CCNB2 as a cross-species hepatotoxicant-responsive candidate gene**
(crucian carp NTH model → rat DrugMatrix GSE57815 → mouse GSE44783).

All results reported in the accompanying manuscript were produced by the scripts marked **[main]** below,
run in the listed order with R 4.6.0 (package versions in `sessionInfo.txt` of the local run).

## Cohorts & roles
| Cohort | Accession / source | Role |
|---|---|---|
| Crucian carp (NTH experiment) | own RNA-seq | Discovery (candidate genes) |
| Rat liver, DrugMatrix | **GSE57815** (GPL1355) | Screening: 9-classifier consensus ranking |
| Mouse liver, 13 model compounds | **GSE44783** (GPL1261), PMID 24040119 | External consistency (single-gene, pre-specified) |

## Label rule (rat)
Labels verified sample-by-sample against GEO dose records: **numeric dose == 0 mg/kg → vehicle
control (n = 279); dose > 0 → exposed (n = 1,939, incl. low-dose arms 0.01–0.65 mg/kg)**.
Original file labels are reproduced exactly by this rule.

## Main pipeline (run in order)
1. `01_prep.R` — build labelled expression object (`rat_expr45.rds`), label audit
2. `02_analysis.R` *(exploratory)*
3. `03_mouse_cv.R`, `03b_loco.R`, `04_rat_cv.R`, `05_verify.R` *(exploratory / verification of the merge-alignment bug)*
4. `06_ccnb2.R` — single-gene CCNB2 stats (rat L1/L3, time/vehicle strata; mouse all/time/13 compounds)
5. `09_fish_summary.R` — three-species summary table + forest plot
6. `training_v2.R` **[main]** — 9 classifiers (RF/SVM/XGB/GLM/GBM/KNN/NNET/LASSO/DT), 5×3 CV,
   permutation importance + ranks (training side), SHAP beeswarm, ROC, figures, `models_for_figures.rds`
7. `nested_cv.R` — fully nested CV (outer 5×2, inner 5-fold re-tuning)
8. `calibration.R` — calibration curves + ECE from CV out-of-fold predictions
9. `youden_ccnb2.R` — single-gene Youden-optimal sensitivity/specificity (rat & mouse)
10. `fig_importance_bars.R`, `fig_shap_polar.R` — importance figure + SHAP polar/beeswarm/bar (per model)
11. `loco_rerun.R`, `loco_dot.R` — leave-one-compound-out (mouse, control-half holdout) + figure
12. `appendix_tables.R`, `format_supp_tables.R` — Supplementary Tables S1–S6 (with `#` header notes)
13. `s5_sensitivity.R` — Supplementary Table S5 (mouse time × sex × per-compound)
14. `audit_all_groups.R` — independent group/label audit of all input files

`driver_rerun.R` runs steps 1–12 end-to-end. `audit_dose0_bug.R`, `test_imp_orig.R`,
`replot_all_figures.R` are internal diagnostics.

## Data
Raw expression matrices are large and must be downloaded from GEO / the original DrugMatrix source
(`normalize.txt` for GSE57815 samples; `mouse_expr.csv` for GSE44783); sample-level metadata tables
used for label definition are derived from GEO series GSM records.

## Methods summary (for full details see the manuscript)
- Labels: GEO numeric dose (see above); alignment by explicit `match()` + `stopifnot(identical)`
  (never `merge()`, which reorders rows).
- Screening: 9 classifiers, fixed documented hyperparameters, seed 123; internal discrimination via
  stratified 5-fold × 3 CV (AUC 0.95–0.99) and fully nested CV (AUC 0.948–0.993); calibration
  ECE ≤ 0.027.
- Consensus ranking: permutation importance (training side, n = 500) + univariate AUC per gene;
  **CCNB2 ranked first in all nine classifiers and in univariate AUC (mean rank 1.0/45)**.
- External consistency: CCNB2 up-regulated in mouse (Δ = +0.944 log2, AUC 0.741), 10/13 compounds,
  LOCO mean AUC 0.771 (single gene) — mouse cohort never used for ranking/tuning.
- Winner's-curse caveat: the rat single-gene AUC is descriptive (ranking cohort); independent
  confirmation is the pre-specified mouse analysis.
