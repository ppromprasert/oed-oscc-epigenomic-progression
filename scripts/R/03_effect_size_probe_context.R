# Join limma statistics to beta-value effect sizes and probe annotations.

source("scripts/R/00_setup.R")

beta_path <- file.path(DATA_DIR, "synthetic", "ref2_beta_matrix.csv")
meta_path <- file.path(DATA_DIR, "synthetic", "ref2_metadata.csv")
limma_path <- file.path(TABLE_DIR, "ref2_primary_limma_all_cpgs.csv")

if (!file.exists(limma_path)) stop("Run 02_fit_primary_limma.R first.")

beta_tbl <- read_csv(beta_path, show_col_types = FALSE)
meta <- read_csv(meta_path, show_col_types = FALSE)
stats <- read_csv(limma_path, show_col_types = FALSE)

probe_col <- names(beta_tbl)[1]
beta_mat <- beta_tbl |>
  column_to_rownames(probe_col) |>
  as.matrix()

meta <- meta |>
  filter(sample_id %in% colnames(beta_mat)) |>
  arrange(match(sample_id, colnames(beta_mat)))
beta_mat <- beta_mat[, meta$sample_id, drop = FALSE]

effects <- tibble(
  Probe = rownames(beta_mat),
  delta_beta = apply(
    beta_mat, 1,
    function(x) delta_beta(x, meta$Progression)
  ),
  progressor_mean_beta = apply(
    beta_mat[, meta$Progression == "Progressor", drop = FALSE],
    1, mean, na.rm = TRUE
  ),
  nonprogressor_mean_beta = apply(
    beta_mat[, meta$Progression == "NonProgressor", drop = FALSE],
    1, mean, na.rm = TRUE
  )
)

out <- stats |>
  left_join(effects, by = "Probe") |>
  arrange(adj.P.Val, P.Value)

write_result(out, "ref2_limma_with_beta_effect_sizes.csv")
