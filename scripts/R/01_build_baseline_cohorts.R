# Build Ref1 / Ref2 / Ref3 baseline cohorts and verify sample alignment.
#
# Public portfolio inputs use privacy-safe metadata. Restricted source metadata
# and institutional paths should never be committed.

source("scripts/R/00_setup.R")

meta_path <- file.path(DATA_DIR, "synthetic", "all_sample_metadata.csv")
if (!file.exists(meta_path)) stop("Missing synthetic metadata: ", meta_path)

meta <- read_csv(meta_path, show_col_types = FALSE)

# Expected fields can be mapped to the real project's baseline flags.
ref_fields <- c(
  Ref1 = "ref1_Lowest_Before_Thislesion",
  Ref2 = "ref2_Lowest_Ever_Thislesion",
  Ref3 = "ref3_Lowest_AnyDx_Before_Thislesion"
)

available <- ref_fields[ref_fields %in% names(meta)]

if (length(available) == 0) {
  stop("No documented baseline-definition fields found in metadata.")
}

build_ref <- function(field, label) {
  meta |>
    filter(.data[[field]] == 1) |>
    mutate(BaselineDefinition = label)
}

cohorts <- imap(available, ~ build_ref(.x, .y))

cohort_summary <- imap_dfr(cohorts, function(x, label) {
  x |>
    count(Progression, name = "N") |>
    mutate(BaselineDefinition = label)
})

write_result(cohort_summary, "baseline_cohort_counts.csv")
