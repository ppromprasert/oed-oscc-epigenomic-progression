# Exposure-associated methylation signature scoring.
#
# Scores are standardized epigenetic patterns derived from literature-linked
# CpGs. They must not be described as direct measurements of environmental dose.

source("scripts/R/00_setup.R")

beta_path <- file.path(DATA_DIR, "public", "all_sample_beta_matrix.csv")
sig_path <- file.path(DATA_DIR, "public", "exposure_signature_cpgs.csv")

if (!file.exists(beta_path) || !file.exists(sig_path)) {
  stop(
    "Provide privacy-safe/public all-sample beta matrix and exposure signature ",
    "definition table before running this script."
  )
}

beta_tbl <- read_csv(beta_path, show_col_types = FALSE)
sig <- read_csv(sig_path, show_col_types = FALSE)

# Expected signature columns:
# Exposure, Signature, Probe, Direction
needed <- c("Exposure", "Signature", "Probe", "Direction")
missing <- setdiff(needed, names(sig))
if (length(missing) > 0) {
  stop("Missing signature columns: ", paste(missing, collapse = ", "))
}

probe_col <- names(beta_tbl)[1]
beta_mat <- beta_tbl |>
  column_to_rownames(probe_col) |>
  as.matrix()

available <- sig |>
  filter(Probe %in% rownames(beta_mat))

if (nrow(available) == 0) {
  stop("None of the exposure-signature CpGs are present in the beta matrix.")
}

# Z-score each CpG across samples, orient by published direction, then average.
# A score is calculated from the CpGs available for that sample/signature.
z_mat <- t(scale(t(beta_mat[unique(available$Probe), , drop = FALSE])))

scores <- available |>
  group_by(Exposure, Signature) |>
  group_modify(function(df, key) {
    probes <- unique(df$Probe)
    dir_lookup <- df |>
      distinct(Probe, Direction) |>
      arrange(match(Probe, probes))

    dirs <- ifelse(dir_lookup$Direction >= 0, 1, -1)
    m <- z_mat[probes, , drop = FALSE]
    oriented <- sweep(m, 1, dirs, `*`)

    tibble(
      sample_id = colnames(oriented),
      score = colMeans(oriented, na.rm = TRUE),
      n_CpGs_used = colSums(!is.na(oriented))
    )
  }) |>
  ungroup()

coverage <- available |>
  count(Exposure, Signature, name = "n_CpGs_available")

write_result(scores, "exposure_signature_scores_long.csv")
write_result(coverage, "exposure_signature_probe_coverage.csv")
