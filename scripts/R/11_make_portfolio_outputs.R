# Generate compact public figures from aggregate/privacy-safe tables.

source("scripts/R/00_setup.R")

exposure_path <- file.path(TABLE_DIR, "ref2_exposure_progression_tests.csv")
corr_path <- file.path(TABLE_DIR, "exposure_priority_promoter_correlations.csv")

if (file.exists(exposure_path)) {
  x <- read_csv(exposure_path, show_col_types = FALSE)

  p <- x |>
    mutate(
      neglog10_FDR = -log10(pmax(FDR, 1e-300)),
      Exposure = reorder(Exposure, delta_median)
    ) |>
    ggplot(aes(delta_median, neglog10_FDR)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_point() +
    geom_text(aes(label = Exposure), check_overlap = TRUE, nudge_y = 0.03) +
    labs(
      title = "Exposure-associated methylation scores by OED progression",
      x = "Progressor - NonProgressor median standardized score",
      y = "-log10(BH FDR)"
    ) +
    theme_classic()

  ggsave(
    file.path(FIG_DIR, "ref2_exposure_effect_significance.png"),
    p, width = 7, height = 5, dpi = 300
  )
}

if (file.exists(corr_path)) {
  ctab <- read_csv(corr_path, show_col_types = FALSE)

  p2 <- ctab |>
    ggplot(aes(Gene, Exposure, fill = spearman_rho)) +
    geom_tile() +
    labs(
      title = "Exposure-score and priority-promoter methylation correlations",
      x = "Priority gene",
      y = "Exposure-associated score",
      fill = "Spearman rho"
    ) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(
    file.path(FIG_DIR, "exposure_priority_promoter_correlation_heatmap.png"),
    p2, width = 8, height = 6, dpi = 300
  )
}
