# driver_rerun.R : run the full reported pipeline end-to-end (rat screening + cross-cohort consistency)
# NOTE: every script assumes the working directory contains inputs/ and outputs/.
#       Adjust PROJ below (or the setwd() line at the top of each script) to your own layout.
PROJ <- "D:/MLroute2"

steps <- c(
  "01_prep.R",             # rebuild labelled expression object (rat_expr45.rds) + label audit
  "05_verify.R",           # alignment verification (match/stopifnot, never merge())
  "06_ccnb2.R",            # single-gene CCNB2, rat strata + mouse
  "09_fish_summary.R",     # three-species summary + forest plot
  "training_v3.R",         # 9 classifiers: 5x3 CV, permutation importance, SHAP, ROC + ECE  (~40 min)
  "nested_cv_v3.R",        # compound-grouped nested CV + cluster bootstrap
  "youden_ccnb2.R",        # Youden-optimal sensitivity/specificity (rat & mouse)
  "loco_rerun.R", "loco_dot.R",          # leave-one-compound-out (mouse) + figure
  "fig_importance_bars.R", "fig_shap_polar.R",  # importance + SHAP polar/beeswarm/bar
  "fig_grouped_cv.R",      # Figure 3C + per-fold stability figure
  "fig3H.R",               # Figure 3H: consistency of CCNB2 induction across rat strata
  "fig_style2.R",          # Figure 3K forest plot + journal-style panels
  "s5_sensitivity.R",      # mouse sensitivity strata (time x sex) + per compound
  "appendix_tables.R", "format_supp_tables.R",  # supplementary tables
  "ortholog_map.R",        # cross-species orthologue mapping (NCBI HomoloGene)
  "audit_all_groups.R"     # independent label/group audit of all input files
)

t0 <- Sys.time()
for (s in steps) {
  cat("\n################################################\n")
  cat("STEP:", s, "| start:", format(Sys.time(), "%H:%M:%S"), "\n")
  cat("################################################\n")
  ok <- tryCatch({ source(file.path(PROJ, s), local = FALSE, echo = FALSE); TRUE },
                 error = function(e) { cat("\n!!! STEP FAILED:", s, "|", conditionMessage(e), "\n"); FALSE })
  if (!ok) { cat("ABORT at", s, "\n"); quit(status = 1) }
  cat("STEP DONE:", s, "|", format(Sys.time(), "%H:%M:%S"), "\n")
}
cat("\nALL STEPS COMPLETED in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
