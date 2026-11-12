# =============================================================================
# 00_prepare_upload.R -- build a clean directory for the code deposit
#
# The working archive accumulates intermediate objects in outputs/cache/.
# Some of them (01_split.rds, 02_imputation.rds) hold the individual-cycle
# records in row-level form, because every later stage needs them. They must
# not be uploaded: the deposit distributes code, not clinical data.
#
# This script copies the archive to a sibling directory, leaving out
#   outputs/cache/            row-level data and intermediate objects
#   data/*.xlsx, data/*.csv   the analysis dataset
#   .git/, .Rhistory, .RData  version-control and session leftovers
# and keeps everything a reader needs to re-run the analysis:
#   code, configuration, documentation, metadata,
#   outputs/figures, outputs/tables, outputs/models (the Shiny bundle)
#
# Usage
#   Rscript 00_prepare_upload.R                 # writes ../<name>_upload/
#   Rscript 00_prepare_upload.R "D:/tmp/deposit"
# =============================================================================

ROOT <- Sys.getenv("IVF_ROOT", unset = getwd())
if (!file.exists(file.path(ROOT, "config.R"))) {
  stop("config.R not found in '", ROOT, "'. Run from the archive root or set ",
       "IVF_ROOT.", call. = FALSE)
}
# config.R also installs the UTF-8 locale guard, without which file.path() on a
# non-ASCII directory fails with "unable to translate ... to native encoding".
source(file.path(ROOT, "config.R"))
ROOT <- normalizePath(ROOT)

args <- commandArgs(trailingOnly = TRUE)
DEST <- if (length(args) >= 1) args[[1]] else
  file.path(dirname(ROOT), paste0(basename(ROOT), "_upload"))
DEST <- normalizePath(DEST, mustWork = FALSE)

if (identical(DEST, ROOT)) {
  stop("Destination must differ from the archive root.", call. = FALSE)
}

message("source      : ", ROOT)
message("destination : ", DEST)

# Directories copied wholesale, and paths excluded from them.
keep_dirs <- c("R", "shiny", "docs", "data",
               file.path("outputs", "figures"),
               file.path("outputs", "tables"),
               file.path("outputs", "models"))

keep_files <- c("00_run_all.R", "00_prepare_upload.R", "config.R",
                "README.md", "DESCRIPTION", "CITATION.cff", "LICENSE",
                ".gitignore",
                file.path("outputs", "session_info.txt"),
                file.path("outputs", "figures", ".gitkeep"),
                file.path("outputs", "tables", ".gitkeep"),
                file.path("outputs", "models", ".gitkeep"))

# Files that must never travel, whatever directory they sit in.
junk_pattern <- "(^\\.[Dd][Ss]_Store$|^Thumbs\\.db$|~$|\\.swp$|^\\.Rhistory$|^\\.RData$)"

unlink(DEST, recursive = TRUE, force = TRUE)
dir.create(DEST, recursive = TRUE, showWarnings = FALSE)

copied <- character(0)

copy_one <- function(rel) {
  from <- file.path(ROOT, rel)
  if (!file.exists(from)) return(invisible(NULL))
  to <- file.path(DEST, rel)
  dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(from, to, overwrite = TRUE)) {
    warning("could not copy ", rel, call. = FALSE)
  } else {
    copied <<- c(copied, rel)
  }
  invisible(NULL)
}

for (d in keep_dirs) {
  full <- file.path(ROOT, d)
  if (!dir.exists(full)) next
  files <- list.files(full, recursive = TRUE, all.files = TRUE,
                      no.. = TRUE, full.names = FALSE)
  for (f in files) {
    rel <- file.path(d, f)
    base <- basename(f)

    # data/ ships the schema document only, never the dataset.
    if (d == "data" && !identical(base, "README.md")) next
    if (grepl(junk_pattern, base)) next
    copy_one(rel)
  }
}

for (f in keep_files) {
  if (grepl(junk_pattern, basename(f))) next
  copy_one(f)
}

# Empty directories are not preserved by file.copy(); recreate the ones the
# scripts expect to write into.
for (d in c(file.path("outputs", "cache"), file.path("outputs", "figures"),
            file.path("outputs", "tables"), file.path("outputs", "models"))) {
  dir.create(file.path(DEST, d), recursive = TRUE, showWarnings = FALSE)
}

# --- report ---------------------------------------------------------------
files <- list.files(DEST, recursive = TRUE, all.files = TRUE, no.. = TRUE)
bytes <- sum(file.info(file.path(DEST, files))$size, na.rm = TRUE)

cat("\n")
cat(sprintf("%-52s %d\n", "files copied", length(files)))
cat(sprintf("%-52s %.1f MB\n", "total size", bytes / 1024^2))
cat("\nTop-level contents:\n")
print(sort(list.files(DEST)))

risky <- c(
  list.files(file.path(DEST, "outputs", "cache"), all.files = TRUE, no.. = TRUE),
  list.files(file.path(DEST, "data"), pattern = "\\.(xlsx|xls|csv)$")
)
if (length(risky) > 0) {
  stop("the deposit still contains row-level or source data: ",
       paste(risky, collapse = ", "), call. = FALSE)
}
cat("\nNo cached row-level objects and no dataset in the deposit. Good to upload.\n")
cat("Zip the directory, or upload its contents directly to GitHub.\n")
