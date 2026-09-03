# Ref2 Progressor vs NonProgressor comparisons of exposure-associated scores.

source("scripts/R/00_setup.R")

score_path <- file.path(TABLE_DIR, "exposure_signature_scores_long.csv")
meta_path <- file.path(DATA_DIR, "public", "all_sample_metadata.csv")

scores <- read_csv(score_path, show_col_types = FALSE)
meta <- read_csv(meta_path, show_col_types = FALSE)

ref2 <- meta |>
  filter(ref2_Lowest_Ever_Thislesion == 1) |>
  select(sample_id, Progression)

dat <- scores |>
  inner_join(ref2, by = "sample_id")

res <- dat |>
  group_by(Exposure, Signature) |>
  group_modify(function(df, key) {
    p <- df$score[df$Progression == "Progressor"]
    np <- df$score[df$Progression == "NonProgressor"]

    if (length(p) < 2 || length(np) < 2) return(tibble())

    test <- wilcox.test(p, np, exact = FALSE)

    tibble(
      n_progressor = length(p),
      n_nonprogressor = length(np),
      progressor_median = median(p, na.rm = TRUE),
      nonprogressor_median = median(np, na.rm = TRUE),
      delta_median = median(p, na.rm = TRUE) - median(np, na.rm = TRUE),
      wilcox_p = test$p.value
    )
  }) |>
  ungroup() |>
  mutate(FDR = p.adjust(wilcox_p, method = "BH")) |>
  arrange(FDR, wilcox_p)

write_result(res, "ref2_exposure_progression_tests.csv")
