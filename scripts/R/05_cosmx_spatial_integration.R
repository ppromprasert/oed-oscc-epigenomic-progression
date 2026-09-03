# Candidate methylation analysis for genes prioritized by spatial/CosMx work.

source("scripts/R/00_setup.R")

stats_path <- file.path(TABLE_DIR, "ref2_limma_with_beta_effect_sizes.csv")
ann_path <- file.path(DATA_DIR, "synthetic", "cosmx_annotation.csv")

stats <- read_csv(stats_path, show_col_types = FALSE)
ann <- read_csv(ann_path, show_col_types = FALSE)

cosmx <- stats |>
  inner_join(ann, by = "Probe") |>
  arrange(adj.P.Val, P.Value)

gene_summary <- cosmx |>
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

write_result(cosmx, "cosmx_candidate_cpg_results.csv")
write_result(gene_summary, "cosmx_candidate_gene_summary.csv")

# Promoter-context summaries should use explicit annotation rather than
# assuming all probes assigned to a gene are promoter probes.
if ("Context" %in% names(cosmx)) {
  promoter <- cosmx |>
    filter(str_detect(Context, "TSS200|TSS1500|5UTR|1stExon"))
  write_result(promoter, "cosmx_promoter_cpg_results.csv")
}
