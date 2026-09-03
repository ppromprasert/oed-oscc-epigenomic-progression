# Robustness checks for the small and imbalanced Ref2 cohort.
#
# Includes:
# - nonparametric Progressor vs NonProgressor comparisons
# - large-effect delta-beta summaries
# - leave-one-NonProgressor-out direction stability
#
# These analyses complement, rather than replace, the primary limma model.

source("scripts/R/00_setup.R")

beta_path <- file.path(DATA_DIR, "synthetic", "ref2_candidate_beta.csv")
meta_path <- file.path(DATA_DIR, "synthetic", "ref2_metadata.csv")

beta_tbl <- read_csv(beta_path, show_col_types = FALSE)
meta <- read_csv(meta_path, show_col_types = FALSE)

probe_col <- names(beta_tbl)[1]
beta_mat <- beta_tbl |>
  column_to_rownames(probe_col) |>
  as.matrix()

required_meta <- c("sample_id", "Progression")
missing <- setdiff(required_meta, names(meta))
if (length(missing) > 0) {
  stop("Missing metadata columns: ", paste(missing, collapse = ", "))
}

meta <- meta |>
  filter(sample_id %in% colnames(beta_mat)) |>
  arrange(match(sample_id, colnames(beta_mat)))

beta_mat <- beta_mat[, meta$sample_id, drop = FALSE]

# Nonparametric candidate-CpG sensitivity analysis.
wilcox_res <- map_dfr(rownames(beta_mat), function(probe) {
  x <- beta_mat[probe, ]
  g <- meta$Progression

  if (
    sum(g == "Progressor", na.rm = TRUE) < 2 ||
    sum(g == "NonProgressor", na.rm = TRUE) < 2
  ) {
    return(NULL)
  }

  test <- wilcox.test(
    x[g == "Progressor"],
    x[g == "NonProgressor"],
    exact = FALSE
  )

  tibble(
    Probe = probe,
    delta_beta = delta_beta(x, g),
    wilcox_p = test$p.value
  )
}) |>
  mutate(FDR = p.adjust(wilcox_p, method = "BH"))

write_result(wilcox_res, "candidate_wilcoxon_sensitivity.csv")

# CpGs with an absolute Progressor - NonProgressor beta difference >= 0.10.
large_effect <- wilcox_res |>
  mutate(abs_delta_beta = abs(delta_beta)) |>
  filter(abs_delta_beta >= 0.10) |>
  arrange(desc(abs_delta_beta))

write_result(large_effect, "candidate_delta_beta_ge_10pct.csv")

# Leave-one-NonProgressor-out sensitivity.
# With only three NonProgressors, this is used to assess directional stability,
# not to claim independent replication.
full_delta <- tibble(
  Probe = rownames(beta_mat),
  full_delta_beta = apply(
    beta_mat,
    1,
    function(x) delta_beta(x, meta$Progression)
  )
)

np_samples <- meta |>
  filter(Progression == "NonProgressor") |>
  pull(sample_id)

loo_long <- map_dfr(np_samples, function(excluded_sample) {
  keep <- meta$sample_id != excluded_sample
  meta_i <- meta[keep, , drop = FALSE]
  beta_i <- beta_mat[, meta_i$sample_id, drop = FALSE]

  tibble(
    Probe = rownames(beta_i),
    excluded_nonprogressor = excluded_sample,
    loo_delta_beta = apply(
      beta_i,
      1,
      function(x) delta_beta(x, meta_i$Progression)
    )
  )
}) |>
  left_join(full_delta, by = "Probe") |>
  mutate(
    direction_match = case_when(
      is.na(full_delta_beta) | is.na(loo_delta_beta) ~ NA,
      full_delta_beta == 0 & loo_delta_beta == 0 ~ TRUE,
      sign(full_delta_beta) == sign(loo_delta_beta) ~ TRUE,
      TRUE ~ FALSE
    )
  )

loo_summary <- loo_long |>
  group_by(Probe) |>
  summarise(
    full_delta_beta = first(full_delta_beta),
    n_leave_one_out = sum(!is.na(loo_delta_beta)),
    direction_concordance = mean(direction_match, na.rm = TRUE),
    min_loo_delta_beta = min(loo_delta_beta, na.rm = TRUE),
    max_loo_delta_beta = max(loo_delta_beta, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(direction_concordance), desc(abs(full_delta_beta)))

write_result(loo_long, "leave_one_nonprogressor_out_long.csv")
write_result(loo_summary, "leave_one_nonprogressor_out_summary.csv")
