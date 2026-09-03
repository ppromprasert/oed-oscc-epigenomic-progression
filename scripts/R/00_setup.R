# Shared setup for the public OED-to-OSCC methylation portfolio workflow.

options(stringsAsFactors = FALSE)
set.seed(43)

pkgs <- c(
  "dplyr", "tidyr", "readr", "stringr", "purrr", "tibble",
  "ggplot2", "limma"
)

missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Install required packages: ", paste(missing, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(limma)
})

ROOT <- normalizePath(".", mustWork = FALSE)
DATA_DIR <- file.path(ROOT, "data")
RESULTS_DIR <- file.path(ROOT, "results")
TABLE_DIR <- file.path(RESULTS_DIR, "tables")
FIG_DIR <- file.path(RESULTS_DIR, "figures")

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

write_result <- function(x, filename) {
  readr::write_csv(x, file.path(TABLE_DIR, filename))
}

delta_beta <- function(beta, group) {
  mean(beta[group == "Progressor"], na.rm = TRUE) -
    mean(beta[group == "NonProgressor"], na.rm = TRUE)
}
