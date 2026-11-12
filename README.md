# Development and internal validation of a random forest model for low or failed fertilization in conventional IVF

This repository contains the complete analysis code for a retrospective cohort study of
4,635 first conventional IVF (c-IVF) cycles. It develops and internally validates a
random forest model for low fertilization rate (LFR) or total fertilization failure
(TFF), defined as a 2PN fertilization rate below 30% of mature oocytes.

The pipeline covers cohort assembly, multiple imputation, collinearity diagnostics,
two-stage model selection across nine algorithms, internal validation, SHAP-based
interpretation, sensitivity analyses and the generation of every figure and table in
the manuscript and its supplementary material.

---

## 1. What is in this archive

```
.
├── 00_run_all.R                      Driver: runs every stage in order
├── 00_prepare_upload.R               Builds a deposit directory without the
│                                     row-level cached objects
├── config.R                          Paths, seeds, constants, variable set,
│                                     final model specification
├── R/
│   ├── utils.R                       Shared helpers (logging, folds, metrics, IO)
│   ├── 01_data_preparation.R         Analysis dataset, factor levels, 80:20 split, Table 1
│   ├── 02_multiple_imputation.R      30 PMM imputations, plus 30 RF imputations
│   ├── 03_collinearity.R             GVIF, condition index, Spearman correlation
│   ├── 04_permutation_importance.R   DALEX permutation importance, nine algorithms
│   ├── 05_variable_selection_nested_cv.R  Predictor count 21→3 by nested CV
│   ├── 06_hyperparameter_search.R    200-point random search × 30 imputations × 10 folds
│   ├── 07_final_model.R              Final random forest, apparent and test performance
│   ├── 08_model_comparison.R         Nine-algorithm comparison, DeLong, IDI/NRI
│   ├── 09_internal_validation.R      Bootstrap .632+, repeated 10-fold CV
│   ├── 10_risk_stratification.R      Two cut-offs, three risk tiers
│   ├── 11_shap.R                     SHAP values and Figure 4
│   ├── 12_sensitivity_analyses.R     RF imputation, complete case, pooled odds ratios
│   ├── 13_figures.R                  Figures 2, 3, S4, S5
│   ├── 14_session_info.R             Software versions actually used
│   └── 15_export_model_bundle.R      Builds the artefact used by the Shiny app
├── shiny/
│   └── app.R                         Interactive risk calculator (English UI)
├── data/
│   └── README.md                     Required input schema (data not distributed)
├── outputs/
│   ├── figures/                      Generated figures (PDF + PNG)
│   ├── tables/                       Generated tables (CSV)
│   ├── models/                       Fitted model bundle (written by stage 15)
│   └── cache/                        Intermediate objects shared between stages
└── docs/
    ├── ANALYSIS_PIPELINE.md          Stage-by-stage description and design notes
    ├── DATA_AVAILABILITY.md          Why the cohort is not deposited, and the
    │                                 three routes available
    └── variable_dictionary.csv       Full variable dictionary
```

---

## 2. Requirements

* **R 4.5.1**, the version used for the published analysis. R >= 4.3 should work.
* A C/C++ toolchain may be needed to compile `lightgbm`, `xgboost` or `nnet`
  if binary packages are unavailable for your platform.
* Memory: the full pipeline holds 30 imputed copies of the development cohort plus
  the tuning grids; 16 GB is comfortable.

### R packages

```r
pkgs <- c(
  # modelling framework
  "tidymodels", "rlang",
  # model engines
  "ranger", "xgboost", "lightgbm", "bonsai", "glmnet", "kernlab", "kknn",
  "nnet", "rpart",
  # data handling
  "dplyr", "tidyr", "tibble",
  # imputation
  "mice",
  # explainability
  "DALEX", "fastshap",
  # discrimination and reclassification
  "pROC", "PredictABEL", "nricens",
  # collinearity
  "car", "olsrr", "corrplot",
  # graphics
  "ggplot2", "patchwork",
  # parallelism
  "foreach", "doParallel", "doSNOW",
  # application
  "shiny"
)
install.packages(setdiff(pkgs, rownames(installed.packages())))
```

`lightgbm` may need `install.packages("lightgbm", repos = "https://dmlc.r-universe.dev")`
on some platforms. Exact versions are recorded in `outputs/session_info.txt` once
stage 14 has been run.

---

## 3. Input data

The pipeline starts from the analysis dataset: one row per cycle, the 21 predictors
listed in Table 1 of the article and the outcome. The individual-cycle file is **not**
included, because it contains identifiable clinical records. `data/README.md` documents
the exact schema the scripts expect — 22 named columns, with
`fertilization_outcome` coded `IVF` (2PN rate ≥ 30%) or `LFR/TFF` (2PN rate < 30%). The
two category names are set once in `config.R` (`EVENT_LEVEL`, `NONEVENT_LEVEL`) and every
later stage reads them from there, including the Shiny app.

Categorical columns are stored as words rather than as numeric codes, so the file can be
read without a code book; `config.R` fixes the level order used for modelling, which also
determines the reference level of each factor.

```bash
export IVF_DATA=/path/to/analysis_dataset.csv   # or set it in .Renviron
export IVF_ROOT=/path/to/this/archive           # only needed if not running from the root
Rscript 00_run_all.R
```

The scripts read `data/analysis_dataset.csv` by default. Any file that matches the schema
can be substituted through `IVF_DATA`; `01_data_preparation.R` stops with a readable
error if a column is missing or a categorical level is not one of the expected ones.

`docs/DATA_AVAILABILITY.md` sets out why the cohort itself is not deposited, the three
sharing routes that are available, and the wording of a data availability statement that
matches this archive. In short: the ethics approval covers analysis of the records, not
their redistribution, and an open deposit is also a cross-border transfer of clinical
data. A synthetic file matching the schema can be added later if a reader needs to
execute the pipeline without access to the real cohort.

---

## 4. Running the analysis

Run everything:

```bash
Rscript 00_run_all.R
```

Run selected stages (useful during development, since stages 5, 6 and 8 take hours):

```bash
Rscript 00_run_all.R 1 2 3 4      # up to permutation importance
Rscript 00_run_all.R 6 7          # hyperparameter search, final model
Rscript 00_run_all.R 13           # redraw all figures from cached results
```

From an interactive session:

```r
setwd("path/to/this/archive")
source("R/01_data_preparation.R")   # ... or source("00_run_all.R")
```

### Checking that the code runs

A full run spends days inside stages 5, 6 and 8. To confirm that the pipeline is
wired correctly before committing to that, run it once at reduced size:

```bash
IVF_SMOKE=1 Rscript 00_run_all.R
```

This shrinks the number of imputations (30 → 2), folds (10 → 3), predictors
(21→3 → 5→3), search size (200 → 6) and algorithms (9 → LR and RF), and writes
everything to `outputs_smoke/` instead of `outputs/`. The real results are not
touched. The whole chain finishes in about twelve minutes. **The numbers it produces
are meaningless** — it tests the plumbing, not the science.

Each stage reads what it needs from `outputs/cache/` and writes its own results there,
so any stage can be re-run on its own. `13_figures.R` in particular never refits a
model: it draws every panel from cached objects and takes seconds.

### What is already in `outputs/`

The deposit ships the results of the stages that finish in minutes to a couple of hours:
Table 1, the missingness and imputation diagnostics, the collinearity diagnostics and
Fig. S3, Table 3 with the per-imputation predictions, the SHAP figure and table, the
sensitivity analyses, the internal-validation table and the session record. The outputs
of stages 4, 5, 6, 8 and 10 — the importance matrix behind Fig. 2, the predictor-count
curves, the tuning grids, Table 4 with the DeLong and IDI/NRI tables, and the risk-tier
tables — are not shipped, because they take days of computation. The code that produces
them is complete and unchanged; running the pipeline writes them into `outputs/`.

---

## 5. Reproducibility

* **Every stochastic step is seeded explicitly**, with a separate seed per step
  (`config.R`, section 2). Seeds are also set inside each parallel worker, so results
  do not depend on task scheduling or on the number of cores.
* **Fold structures are generated once and reused.** The 80:20 importance split is
  drawn on an index table and applied to all 270 algorithm-by-imputation combinations;
  the outer and inner cross-validation folds are shared by every predictor count within
  an algorithm. Comparisons across algorithms and across model sizes are therefore
  paired.
* **The random grids are drawn from a single seeded call** in a fixed algorithm order
  (`ALGO_ORDER_TUNE`). Changing that order changes the grids; see the note in
  `config.R`.
* **The test set plays no part** in screening, tuning or model selection. It is used
  once, for the held-out evaluation in stage 7, and again in stage 12 for the
  recalibration and sensitivity analyses that are explicitly described as such.
* **Parallelism affects wall-clock time only, never results.**

Two cross-validated estimates appear in the output and they answer different questions:

| Estimate | Where | What it is |
|---|---|---|
| Pooled CV-AUC after tuning (e.g. RF 0.625) | `06_hyperparameter_search.R` | Performance of the selected configuration, evaluated on the development set only. The search is *not* nested inside the outer loop, so this figure is mildly optimistic. |
| Nested-CV AUC per predictor count | `05_variable_selection_nested_cv.R` | Honest comparison *between* predictor counts, because tuning happens inside each outer fold. Use it to choose a model size, not as an absolute performance claim. |
| Held-out test AUC | `07_final_model.R` | The single estimate quoted as the model's internal validation performance. |

---

## 6. Expected runtime

Runtime is dominated by three stages. The numbers below are indicative for a
16-core workstation; the fit counts are exact and let you rescale.

| Stage | Fits performed | Indicative time |
|---|---|---|
| 1–2 | Cohort assembly, 30 + 30 imputations | 5–15 min |
| 3–4 | Collinearity (30 models); permutation importance (9 × 30 fits, B = 500) | 20–60 min |
| 5 | 9 algorithms × 19 predictor counts × 30 imputations × 10 outer folds, each with inner tuning | 8–48 h |
| 6 | 8 algorithms × 200 candidates × 30 imputations × 10 folds ≈ 480 000 fits | 12–72 h |
| 7–8 | Final model (30 fits); out-of-fold predictions (9 × 30 × 10 fits) | 1–6 h |
| 9–12 | Bootstrap, repeated CV, SHAP, sensitivity analyses | 1–4 h |
| 13–15 | Figures, session info, bundle | < 5 min |

Stages 5, 6 and 8 are embarrassingly parallel and use
`min(N_CORES_MAX, detectCores() - 4)` workers. Reduce `N_CORES_MAX` in `config.R` on a
shared machine. Setting `IVF_QUICK=1` reduces the bootstrap and repeat counts in stage 9
for a smoke test, and `IVF_REDRAW=1` rebuilds Figure 4 from the cached SHAP values in
stage 11 without recomputing them (the full run takes about 45 minutes).
`IVF_DEBUG_WORKERS=1` keeps the stdout of the parallel workers, which is discarded by
default and is the only way to see an error raised inside `%dopar%`.

---

## 7. Mapping to the manuscript

| Manuscript item | Produced by |
|---|---|
| Fig. 1 (flow diagram) | Not generated from data; drawn separately |
| Fig. 2 (importance matrix and performance curves) | `13_figures.R`, from `04_permutation_importance.R` and `05_variable_selection_nested_cv.R` |
| Fig. 3 (discrimination, calibration, decision curve, confusion matrix) | `13_figures.R`, from `07_final_model.R` |
| Fig. 4 (SHAP) | `11_shap.R` |
| Fig. S1 (sample-size curve) | Not generated from data; drawn separately |
| Fig. S2 (missingness heatmap) | Not generated from data; drawn separately |
| Fig. S3 (Spearman correlation heatmap) | `03_collinearity.R` |
| Fig. S4 (complete-case ROC) | `13_figures.R`, from `12_sensitivity_analyses.R` |
| Fig. S5 (random-forest imputation ROC) | `13_figures.R`, from `12_sensitivity_analyses.R` |
| Table 1 (baseline characteristics) | `01_data_preparation.R` |
| Table 3 (final model performance) | `07_final_model.R` |
| Table 4 (nine-algorithm comparison) | `08_model_comparison.R` |
| Table S2.2 (selected hyperparameters) | `06_hyperparameter_search.R` |
| Table S (VIF/GVIF) | `03_collinearity.R` |
| Table S (DeLong, IDI/NRI) | `08_model_comparison.R` |
| Table S (internal validation) | `09_internal_validation.R` |
| Table S (risk stratification) | `10_risk_stratification.R` |
| Table S (sensitivity analyses, adjusted ORs) | `12_sensitivity_analyses.R` |
| Variable definitions (the 21 predictors) | `docs/variable_dictionary.csv` |
| Text S2 (tuning procedure) | `config.R`, `05_variable_selection_nested_cv.R`, `06_hyperparameter_search.R` |

`docs/ANALYSIS_PIPELINE.md` gives the full stage-by-stage description, including the
rationale for the two-stage design and the distinction between the two sets of
cross-validated estimates.

Table 1 lists the 21 predictors of the article plus the outcome block, which is the same
set the code models: no variable is derived, merged or dropped after the input file. The
sample sizes are 3,707 in the development set and 928 in the test set throughout, in line
with the manuscript.

---

## 8. Interactive calculator

```r
Rscript R/15_export_model_bundle.R    # once: build the model bundle
shiny::runApp("shiny")
```

The application contains no hard-coded numbers. Thresholds, tier-specific event rates,
recalibration coefficients and the reported AUC all come from
`outputs/models/model_bundle.rds`, so the calculator cannot drift out of step with the
analysis. Zone assignment uses the raw probability; the monotone recalibration changes
only the number displayed.

The model is a research tool. It has not been externally validated, and it is not a
substitute for clinical judgement.

---

## 9. Citation

If you use this code, please cite the article (see `CITATION.cff`) and this archive.

## 10. Licence

MIT. See `LICENSE`.

## 11. Contact

First author: Yiliang Zhu (3311878665@qq.com).
Corresponding authors: Peng Li, PhD, and Chuang Liu, PhD.
Author details are recorded in `CITATION.cff` and `DESCRIPTION`; only the
`repository-code` URL in `CITATION.cff` still needs to be filled in after the
repository is created on GitHub.

## 12. Preparing the deposit

`outputs/cache/01_split.rds` and `outputs/cache/02_imputation.rds` contain the
individual-cycle records in row-level form, because every later stage needs them.
They must never be uploaded, and neither must `analysis_dataset.csv`.

`00_prepare_upload.R` builds a clean copy in a sibling directory, keeping the code,
the documentation, the metadata and the generated figures, tables and model bundle,
and leaving behind every row-level object:

```bash
Rscript 00_prepare_upload.R              # -> ../<archive>_upload/
Rscript 00_prepare_upload.R "D:/tmp/deposit"
```

The script aborts if it finds a cached object or a data file in the destination.
