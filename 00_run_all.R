# =============================================================================
# 00_run_all.R -- run the complete analysis pipeline in order
#
# Usage
#   Rscript 00_run_all.R                 # run every stage
#   Rscript 00_run_all.R 1 2 3           # run selected stages only
#
# Each stage writes its outputs to outputs/ and, where a later stage depends on
# it, to outputs/cache/. Stages can therefore be re-run individually without
# redoing everything.
#
# Stage | Script                                   | Produces
# ------+------------------------------------------+--------------------------
#   1   | R/01_data_preparation.R                  | analysis set, 80:20 split, Table 1
#   2   | R/02_multiple_imputation.R               | 30 imputed datasets (PMM and RF)
#   3   | R/03_collinearity.R                      | VIF/GVIF, condition index, Spearman
#   4   | R/04_permutation_importance.R            | variable ranking, nine algorithms
#   5   | R/05_variable_selection_nested_cv.R      | predictor count per algorithm
#   6   | R/06_hyperparameter_search.R             | Table S2.2, pooled CV-AUC
#   7   | R/07_final_model.R                       | Table 3, apparent and test performance
#   8   | R/08_model_comparison.R                  | Table 4, DeLong, IDI/NRI
#   9   | R/09_internal_validation.R               | bootstrap .632+, repeated CV
#  10   | R/10_risk_stratification.R               | two cut-offs, three risk tiers
#  11   | R/11_shap.R                              | Figure 4, SHAP importance
#  12   | R/12_sensitivity_analyses.R              | Table S, adjusted odds ratios
#  13   | R/13_figures.R                           | Figures 2, 3, S4, S5
#  14   | R/14_session_info.R                      | software versions
#  15   | R/15_export_model_bundle.R               | Shiny model bundle
#
# Note on runtime. Stages 5, 6 and 8 each perform tens of thousands of model
# fits and together account for almost all of the total runtime. For a first
# pass, stages 1-4 and 7 can be run alone; each stage caches what it produces,
# so the expensive ones can be started later without repeating the cheap ones.
# =============================================================================

# -----------------------------------------------------------------------------
# Resolve project root from this script's own location
# -----------------------------------------------------------------------------
# The UTF-8 locale guard normally lives at the top of config.R, but the root has
# to be resolved -- and setwd() has to succeed -- before config.R can be sourced.
# On a Windows machine running a non-UTF-8 locale, setwd() to a directory whose
# name contains non-ASCII characters fails with "cannot change working
# directory", so the guard is repeated here.
if (!grepl("utf-?8", Sys.getlocale("LC_CTYPE"), ignore.case = TRUE)) {
  for (loc in c("Chinese (Simplified)_China.utf8", "Chinese_China.utf8",
                "English_United States.utf8", "en_US.UTF-8", "C.UTF-8")) {
    got <- suppressWarnings(tryCatch(Sys.setlocale("LC_CTYPE", loc),
                                     error = function(e) ""))
    if (nzchar(got)) break
  }
}

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
if (length(file_arg) == 1) {
  ROOT <- dirname(normalizePath(sub("^--file=", "", file_arg), mustWork = FALSE))
} else {
  ROOT <- getwd()
}
Sys.setenv(IVF_ROOT = ROOT)
setwd(ROOT)

source(file.path(ROOT, "config.R"), local = FALSE)

# -----------------------------------------------------------------------------
# Which stages to run
# -----------------------------------------------------------------------------
requested <- suppressWarnings(as.integer(commandArgs(trailingOnly = TRUE)))
if (length(requested) == 0 || anyNA(requested)) requested <- 1:15

stages <- list(
  `1`  = "R/01_data_preparation.R",
  `2`  = "R/02_multiple_imputation.R",
  `3`  = "R/03_collinearity.R",
  `4`  = "R/04_permutation_importance.R",
  `5`  = "R/05_variable_selection_nested_cv.R",
  `6`  = "R/06_hyperparameter_search.R",
  `7`  = "R/07_final_model.R",
  `8`  = "R/08_model_comparison.R",
  `9`  = "R/09_internal_validation.R",
  `10` = "R/10_risk_stratification.R",
  `11` = "R/11_shap.R",
  `12` = "R/12_sensitivity_analyses.R",
  `13` = "R/13_figures.R",
  `14` = "R/14_session_info.R",
  `15` = "R/15_export_model_bundle.R"
)

# -----------------------------------------------------------------------------
# Shared environment
# -----------------------------------------------------------------------------
# All stages share one environment so that objects created early (the imputed
# datasets, for instance) remain available downstream without being written to
# disk as large binaries.
PIPELINE_ENV <- new.env(parent = globalenv())

run_stage <- function(i) {
  script <- stages[[as.character(i)]]
  path <- file.path(ROOT, script)
  if (!file.exists(path)) {
    warning("Stage ", i, " skipped: file not found -> ", script)
    return(invisible(FALSE))
  }
  message("\n", strrep("=", 72))
  message("Stage ", i, ": ", script)
  message(strrep("=", 72))
  t0 <- Sys.time()
  ok <- tryCatch({
    sys.source(path, envir = PIPELINE_ENV, keep.source = FALSE)
    TRUE
  }, error = function(e) {
    message("Stage ", i, " failed: ", conditionMessage(e))
    FALSE
  })
  dt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  message(sprintf("Stage %d %s in %.1f min", i,
                  if (ok) "completed" else "FAILED", dt))
  invisible(ok)
}

# -----------------------------------------------------------------------------
# Execute
# -----------------------------------------------------------------------------
results <- vapply(requested, run_stage, logical(1))

message("\n", strrep("=", 72))
message("Pipeline summary")
message(strrep("=", 72))
for (k in seq_along(requested)) {
  message(sprintf("  stage %2d  %-36s  %s", requested[k],
                  stages[[as.character(requested[k])]],
                  if (results[k]) "ok" else "FAILED"))
}
if (!all(results)) {
  warning("One or more stages failed; see the messages above.")
}
