# Known OSCC candidate-gene analysis and CDKN2A deep dive.

source("scripts/R/00_setup.R")

stats_path <- file.path(TABLE_DIR, "ref2_limma_with_beta_effect_sizes.csv")
ann_path <- file.path(DATA_DIR, "synthetic", "known_oscc_annotation.csv")

if (!file.exists(stats_path)) stop("Run scripts 02-03 first.")

stats <- read_csv(stats_path, show_col_types = FALSE)
ann <- read_csv(ann_path, show_col_types = FALSE)

# Exact probe-to-gene annotation is preferred over substring matching.
known <- stats |>
  inner_join(ann, by = "Probe") |>
  arrange(adj.P.Val, P.Value)

gene_summary <- known |>
  group_by(Gene) |>
  summarise(
    n_CpGs = n(),
    n_nominal = sum(P.Value < 0.05, na.rm = TRUE),
    n_FDR = sum(adj.P.Val < 0.05, na.rm = TRUE),
    mean_logFC = mean(logFC, na.rm = TRUE),
    median_delta_beta = median(delta_beta, na.rm = TRUE),
    min_p = min(P.Value, na.rm = TRUE),
    min_FDR = min(adj.P.Val, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(min_FDR, min_p)

write_result(known, "known_oscc_cpg_results.csv")
write_result(gene_summary, "known_oscc_gene_summary.csv")

cdkn2a <- known |>
  filter(Gene == "CDKN2A") |>
  arrange(P.Value)

write_result(cdkn2a, "cdkn2a_probe_results.csv")
