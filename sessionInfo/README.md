# Software environment

`sessionInfo()` output for every analysis step of this study, as produced by the scripts in the
repository root (R 4.6.0, x86_64-w64-mingw32/x64, Windows 10 x64 build 19045).

| File | Produced by | Covers |
|---|---|---|
| `sessionInfo_training_v3.txt` | `training_v3.R` | 9-classifier screening: caret, randomForest, xgboost, gbm, kernlab (SVM), glmnet, nnet, rpart, pROC, DALEX, ggplot2, data.table |
| `sessionInfo_nested_cv_v3.txt` | `nested_cv_v3.R` | Compound-grouped nested CV + cluster bootstrap: caret, pROC, data.table |
| `sessionInfo_ortholog_map.txt` | `ortholog_map.R` | Cross-species orthologue mapping: homologene 1.4.68.19.3.27, data.table |

## Key package versions

| Package | Version | Used for |
|---|---|---|
| R | 4.6.0 (2026-04-24 ucrt) | all analyses |
| caret | 7.0-1 | resampling, tuning, `preProcess` |
| randomForest | 4.7-1.2 | RF |
| xgboost | 3.2.1.1 | XGBoost (native + caret) |
| gbm | 2.3.1 | GBM |
| kernlab | 0.9-33 | radial-basis SVM |
| glmnet | 5.0 | LASSO |
| nnet | 7.3-20 | single-layer neural network |
| rpart | 4.1.27 | CART |
| pROC | 1.19.0.1 | ROC AUC, DeLong CI |
| DALEX | 2.5.3 | permutation importance |
| homologene | 1.4.68.19.3.27 | orthologue mapping (offline, version-pinned) |
| data.table | 1.18.4 | data handling |
| ggplot2 | 4.0.3 | figures |
| doParallel / foreach / iterators | 1.0.17 / 1.5.2 / 1.0.14 | parallel resampling |

> Versions for components outside the R environment (read quantification, enrichment, docking and MD)
> are given in the manuscript and its Supplementary Material.
