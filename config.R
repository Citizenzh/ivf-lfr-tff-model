# =============================================================================
# config.R -- central configuration for the IVF LFR/TFF prediction pipeline
#
# Every analysis script in this repository sources this file first. It defines
# directory locations, random seeds, analysis constants, the predictor set and
# the final model specification, so that no script contains hard-coded paths or
# duplicated parameter values.
#
# Project : Development and internal validation of a random forest model for
#           low or failed fertilization in conventional IVF
# Licence : MIT (see LICENSE)
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Locale guard
# -----------------------------------------------------------------------------
# If R is running under a non-UTF-8 locale (CP936 on Chinese Windows, or C in a
# bare shell), file paths containing non-ASCII characters cannot be translated
# and readRDS()/saveRDS()/read.csv() fail with "unable to translate ... to
# native encoding". This happens even though the code itself is ASCII. The
# guard switches LC_CTYPE to a UTF-8 locale when one is available and warns
# otherwise, so that the failure is explained rather than mysterious.
if (!grepl("utf-?8", Sys.getlocale("LC_CTYPE"), ignore.case = TRUE)) {
  utf8_candidates <- c(
    "Chinese (Simplified)_China.utf8",
    "Chinese_China.utf8",
    "English_United States.utf8",
    "en_US.UTF-8",
    "C.UTF-8"
  )
  for (loc in utf8_candidates) {
    got <- suppressWarnings(tryCatch(Sys.setlocale("LC_CTYPE", loc),
                                     error = function(e) ""))
    if (nzchar(got)) break
  }
  if (!grepl("utf-?8", Sys.getlocale("LC_CTYPE"), ignore.case = TRUE)) {
    warning("Could not switch to a UTF-8 locale. Scripts will still run if the ",
            "project path is ASCII-only; otherwise set LC_CTYPE manually, e.g. ",
            "Sys.setlocale('LC_CTYPE', 'Chinese (Simplified)_China.utf8').",
            call. = FALSE)
  }
}

# -----------------------------------------------------------------------------
# 1. Directory layout
# -----------------------------------------------------------------------------
# ROOT is resolved from the location of this file, so the archive can be
# unzipped anywhere. Set IVF_ROOT in .Renviron to override.
ROOT <- Sys.getenv("IVF_ROOT", unset = "")
if (!nzchar(ROOT)) {
  # Search the call stack for this file when sourced from 00_run_all.R or
  # from an interactive session.
  for (i in seq_len(sys.nframe())) {
    of <- sys.frame(i)$ofile
    if (!is.null(of) && basename(of) == "config.R") {
      ROOT <- dirname(normalizePath(of, mustWork = FALSE))
      break
    }
  }
}
if (!nzchar(ROOT)) ROOT <- getwd()

DIR_R       <- file.path(ROOT, "R")
DIR_DATA    <- file.path(ROOT, "data")
DIR_OUT     <- file.path(ROOT, "outputs")
DIR_FIG     <- file.path(DIR_OUT, "figures")
DIR_TAB     <- file.path(DIR_OUT, "tables")
DIR_MODEL   <- file.path(DIR_OUT, "models")
DIR_CACHE   <- file.path(DIR_OUT, "cache")
DIR_SHINY   <- file.path(ROOT, "shiny")

for (d in c(DIR_OUT, DIR_FIG, DIR_TAB, DIR_MODEL, DIR_CACHE)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# Input file. The deposit starts from the analysis dataset: one row per cycle,
# the 21 predictors listed in Table 1 and the outcome. The individual records
# are not distributed with this archive; data/README.md gives the schema and the
# coding of every column.
FILE_CSV <- Sys.getenv("IVF_DATA",
                       unset = file.path(DIR_DATA, "analysis_dataset.csv"))

# -----------------------------------------------------------------------------
# 2. Random seeds
# -----------------------------------------------------------------------------
# Each stochastic step has its own seed so that a change in one step does not
# cascade through the whole pipeline.
SEED_SPLIT      <- 20260823  # 80:20 development/test partition, mice
SEED_VI         <- 20260828  # 80:20 split for permutation importance,
                             #   ranger/xgboost/lightgbm engine seeds
SEED_OLDFOLDS   <- 20260808  # base; outer folds use SEED_OLDFOLDS + 1 = 20260809
SEED_TUNE       <- 20260915  # random search grids, second-stage CV folds
SEED_BOOT       <- 20260920  # bootstrap resampling
SEED_VI_WORKER  <- 1000      # worker-level seed offset for permutation importance

# -----------------------------------------------------------------------------
# 3. Analysis constants
# -----------------------------------------------------------------------------
N_IMPUTATIONS   <- 30    # number of multiply imputed datasets
MICE_MAXIT      <- 20    # chained-equations iterations
TRAIN_PROP      <- 0.80  # development share of the cohort

V_OUTER         <- 10    # outer folds
V_OUTER_FALLBACK<- 5     # used when a fold fails the balance check
V_INNER         <- 10
V_INNER_FALLBACK<- 5
FOLD_MIN_N      <- 10    # minimum assessment-fold size
FOLD_MIN_PROP   <- 0.05  # minimum event proportion in an assessment fold

VAR_COUNTS      <- 21:3  # candidate numbers of predictors, decreasing

B_VI            <- 500   # permutations for model_parts
VI_FIT_PROP     <- 0.80  # fitting share inside the importance split

TUNE_GRID_SIZE  <- 200   # random-search candidates per algorithm, stage 2
N_BOOT          <- 500   # bootstrap replicates for confidence intervals
N_REPEATED_CV   <- 10    # repeats for repeated 10-fold CV

N_CORES_MAX     <- 30    # upper bound on worker processes

# Collinearity thresholds (reported in the Methods)
GVIF_ADJ_MAX    <- 5
COND_INDEX_MAX  <- 30
VDP_MAX         <- 0.50
SPEARMAN_MAX    <- 0.70

# Risk-stratification targets (used by 10_risk_stratification.R)
SENS_TARGET     <- 0.90
SPEC_TARGET     <- 0.90

# -----------------------------------------------------------------------------
# 4. Outcome coding
# -----------------------------------------------------------------------------
OUTCOME_VAR   <- "fertilization_outcome"
# The two outcome categories, named as in the manuscript: LFR/TFF (low or
# failed fertilization, 2PN rate < 30%) is the event, conventional IVF (2PN
# rate >= 30%) is the non-event. The event must be the second level so that
# event_level = "second" selects it throughout the pipeline.
EVENT_LEVEL   <- "LFR/TFF"   # level modelled as the event
NONEVENT_LEVEL<- "IVF"       # conventional fertilization
POS_CLASS     <- paste0(".pred_", EVENT_LEVEL)  # probability column returned by tidymodels

# -----------------------------------------------------------------------------
# 5. Predictor set
# -----------------------------------------------------------------------------
# The 21 predictors of Table 1, in the order in which they are tabulated.
PREDICTORS <- c(
  # categorical
  "female_marital_status", "female_infertility_type", "male_infertility_type",
  "stimulation_protocol", "tuberculosis_history", "oligoasthenospermia",
  "teratozoospermia", "infertility_etiology",
  # continuous
  "female_age", "infertility_duration", "antral_follicle_count", "female_bmi",
  "male_bmi", "gn_starting_dose", "gn_duration", "e2_on_hcg_day",
  "lh_on_hcg_day", "progesterone_on_hcg_day", "progressive_motility",
  "abnormal_morphology_rate", "mature_oocytes"
)

CATEGORICAL_PREDICTORS <- c(
  "female_marital_status", "female_infertility_type", "male_infertility_type",
  "stimulation_protocol", "tuberculosis_history", "oligoasthenospermia",
  "teratozoospermia", "infertility_etiology"
)
NUMERIC_PREDICTORS <- setdiff(PREDICTORS, CATEGORICAL_PREDICTORS)

# Level order of every categorical variable, as it appears in Table 1. The
# dataset stores the levels as words, so this only fixes the order.
LEVELS <- list(
  female_marital_status  = c("First marriage", "Remarriage"),
  female_infertility_type= c("Primary infertility", "Secondary infertility"),
  male_infertility_type  = c("Primary infertility", "Secondary infertility"),
  stimulation_protocol   = c("Long protocol", "Short or mild stimulation",
                             "Antagonist"),
  tuberculosis_history   = c("No", "Yes"),
  oligoasthenospermia    = c("No", "Yes"),
  teratozoospermia       = c("No", "Yes"),
  infertility_etiology   = c("Tubal factor", "Diminished ovarian reserve",
                             "Other ovarian factor", "Ovulation disorder",
                             "Pelvic inflammation", "Uterine lesion",
                             "Multiple factors", "Unexplained")
)


# Variables retained by the two-stage procedure in the final random forest.
RF_VARS <- c(
  "e2_on_hcg_day", "infertility_duration", "male_bmi", "infertility_etiology",
  "abnormal_morphology_rate", "mature_oocytes", "female_infertility_type",
  "teratozoospermia", "progressive_motility", "progesterone_on_hcg_day"
)

# Hyperparameters selected for the final random forest (Table S2.2).
RF_BEST_PARAMS <- list(
  mtry            = 1L,
  trees           = 577L,
  min_n           = 40L,
  max.depth       = 11L,
  sample.fraction = 0.8583293
)

# Reference label for every variable, matching Table 1.
VAR_LABELS <- c(
  female_marital_status   = "Female marital status",
  female_infertility_type = "Female infertility type",
  male_infertility_type   = "Male infertility type",
  stimulation_protocol    = "Stimulation protocol",
  tuberculosis_history    = "Tuberculosis history",
  oligoasthenospermia     = "Oligoasthenospermia",
  teratozoospermia        = "Teratozoospermia",
  infertility_etiology    = "Infertility etiology",
  female_age              = "Female age",
  infertility_duration    = "Infertility duration",
  antral_follicle_count   = "Antral follicle count",
  female_bmi              = "Female BMI",
  male_bmi                = "Male BMI",
  gn_starting_dose        = "Gn starting dose",
  gn_duration             = "Gn duration",
  e2_on_hcg_day           = "E2 on hCG day",
  lh_on_hcg_day           = "LH on hCG day",
  progesterone_on_hcg_day = "P on hCG day",
  progressive_motility    = "Pre-wash progressive motility",
  abnormal_morphology_rate= "Abnormal morphology rate",
  mature_oocytes          = "Mature oocytes",
  fertilization_outcome   = "Fertilization outcome"
)

# Units, used when rendering tables.
VAR_UNITS <- c(
  female_age = "years", infertility_duration = "years",
  antral_follicle_count = "n", female_bmi = "kg/m2", male_bmi = "kg/m2",
  gn_starting_dose = "IU", gn_duration = "days", e2_on_hcg_day = "pg/mL",
  lh_on_hcg_day = "IU/L", progesterone_on_hcg_day = "ng/mL",
  progressive_motility = "%", abnormal_morphology_rate = "%",
  mature_oocytes = "n"
)

# -----------------------------------------------------------------------------
# 6. Algorithm registry
# -----------------------------------------------------------------------------
# Names used throughout the pipeline. LR has no tunable hyperparameters.
ALGORITHMS <- c("LR", "ENET", "DT", "RF", "XGB", "GBM", "SVM", "KNN", "MLP")
ALGORITHMS_TUNED <- setdiff(ALGORITHMS, "LR")

# Inner-loop grid size used in the nested cross-validation (Table S2.1).
INNER_GRID_SIZE <- c(
  LR = 0L, ENET = 6L, DT = 5L, RF = 27L,
  XGB = 50L, GBM = 30L, SVM = 10L, KNN = 5L, MLP = 64L
)

# Display order used in figures and tables.
ALGO_ORDER <- c("LR", "ENET", "DT", "RF", "XGB", "GBM", "SVM", "KNN", "MLP")

# Iteration order used when the second-stage random grids are drawn. All nine
# grids are generated from a single call to set.seed(), so this order is part
# of the reproducibility contract: changing it changes which candidate is
# drawn for each algorithm. It follows the order of the original analysis.
ALGO_ORDER_TUNE <- c("LR", "ENET", "RF", "XGB", "GBM", "SVM", "MLP", "KNN", "DT")

# -----------------------------------------------------------------------------
# 7. Session options
# -----------------------------------------------------------------------------
options(
  stringsAsFactors = FALSE,
  dplyr.summarise.inform = FALSE,
  tidymodels.dark = FALSE,
  scipen = 999
)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggplot2::theme_set(ggplot2::theme_bw(base_size = 11))
}


# -----------------------------------------------------------------------------
# 8. Smoke-test mode
# -----------------------------------------------------------------------------
# A full run spends days inside stages 5, 6 and 8. IVF_SMOKE=1 shrinks every
# dimension that makes them slow and redirects all output to outputs_smoke/,
# so the whole chain from stage 1 to stage 15 can be exercised in minutes
# without touching the real results:
#
#   IVF_SMOKE=1 Rscript 00_run_all.R
#
# Numbers produced in this mode are meaningless -- it exists to prove that the
# code runs end to end, not to reproduce the paper.
SMOKE <- identical(Sys.getenv("IVF_SMOKE"), "1")

if (SMOKE) {
  N_IMPUTATIONS    <- 2
  MICE_MAXIT       <- 2
  V_OUTER          <- 3
  V_OUTER_FALLBACK <- 3
  V_INNER          <- 3
  V_INNER_FALLBACK <- 3
  VAR_COUNTS       <- c(5, 3)
  B_VI             <- 5
  TUNE_GRID_SIZE   <- 6
  N_BOOT           <- 10
  N_REPEATED_CV    <- 1
  N_CORES_MAX      <- 2

  ALGORITHMS       <- c("LR", "RF")
  ALGORITHMS_TUNED <- setdiff(ALGORITHMS, "LR")
  ALGO_ORDER       <- intersect(ALGO_ORDER, ALGORITHMS)
  # Intersecting keeps the original relative order, so the grids are still
  # drawn from one seeded call in a fixed sequence.
  ALGO_ORDER_TUNE  <- intersect(ALGO_ORDER_TUNE, ALGORITHMS)

  INNER_GRID_SIZE <- c(
    LR = 0L, ENET = 3L, DT = 2L, RF = 6L,
    XGB = 3L, GBM = 3L, SVM = 3L, KNN = 2L, MLP = 4L
  )

  DIR_OUT   <- file.path(ROOT, "outputs_smoke")
  DIR_FIG   <- file.path(DIR_OUT, "figures")
  DIR_TAB   <- file.path(DIR_OUT, "tables")
  DIR_MODEL <- file.path(DIR_OUT, "models")
  DIR_CACHE <- file.path(DIR_OUT, "cache")
  for (d in c(DIR_OUT, DIR_FIG, DIR_TAB, DIR_MODEL, DIR_CACHE)) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
  }

  message("IVF_SMOKE is set: reduced dimensions and 2 algorithms, ",
          "all output redirected to outputs_smoke/.")
}

message("config.R loaded. ROOT = ", ROOT)
