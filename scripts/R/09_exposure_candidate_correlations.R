# Correlate exposure-associated methylation scores with prioritized promoter
# methylation features across the full methylation cohort.

source("scripts/R/00_setup.R")

score_path <- file.path(TABLE_DIR, "exposure_signature_scores_long.csv")
promoter_path <- file.path(DATA_DIR, "public", "priority_promoter_methylation.csv")

scores <- read_csv(score_path, show_col_types = FALSE)
prom <- read_csv(promoter_path, show_col_types = FALSE)

# Expected promoter table: sample_id, Gene, promoter_value
needed <- c("sample_id", "Gene", "promoter_value")
missing <- setdiff(needed, names(prom))
if (length(missing) > 0) stop("Missing promoter columns: ", paste(missing, collapse = ", "))

dat <- scores |>
  inner_join(prom, by = "sample_id")

corr <- dat |>
  group_by(Exposure, Gene) |>
  group_modify(function(df, key) {
    x <- df$score
    y <- df$promoter_value
    ok <- complete.cases(x, y)
    if (sum(ok) < 10) return(tibble())

    test <- cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)

    tibble(
      spearman_rho = unname(test$estimate),
      p_value = test$p.value,
      n = sum(ok)
    )
  }) |>
  ungroup() |>
  mutate(FDR = p.adjust(p_value, method = "BH")) |>
  arrange(FDR, p_value)

write_result(corr, "exposure_priority_promoter_correlations.csv")
