# Grade-adjusted Ref2 limma model on M-values.
#
# Primary model:
#   M-value ~ Progression + Grade
#
# duplicateCorrelation is used when repeated patient samples are present.

source("scripts/R/00_setup.R")

meta_path <- file.path(DATA_DIR, "synthetic", "ref2_metadata.csv")
m_path <- file.path(DATA_DIR, "synthetic", "ref2_m_matrix.csv")

meta <- read_csv(meta_path, show_col_types = FALSE)
m_tbl <- read_csv(m_path, show_col_types = FALSE)

probe_col <- names(m_tbl)[1]
m_mat <- m_tbl |>
  column_to_rownames(probe_col) |>
  as.matrix()

required_meta <- c("sample_id", "Progression", "Grade")
missing <- setdiff(required_meta, names(meta))
if (length(missing) > 0) {
  stop("Missing metadata columns: ", paste(missing, collapse = ", "))
}

meta <- meta |>
  filter(sample_id %in% colnames(m_mat)) |>
  arrange(match(sample_id, colnames(m_mat)))

m_mat <- m_mat[, meta$sample_id, drop = FALSE]

meta <- meta |>
  mutate(
    Progression = relevel(factor(Progression), ref = "NonProgressor"),
    Grade = factor(Grade)
  )

design <- model.matrix(~ Progression + Grade, data = meta)

if ("patient_id" %in% names(meta) && anyDuplicated(meta$patient_id)) {
  corfit <- duplicateCorrelation(
    m_mat,
    design,
    block = meta$patient_id
  )
  fit <- lmFit(
    m_mat,
    design,
    block = meta$patient_id,
    correlation = corfit$consensus
  )
} else {
  fit <- lmFit(m_mat, design)
}

fit <- eBayes(fit)

coef_name <- grep("^Progression", colnames(design), value = TRUE)[1]
tt <- topTable(
  fit,
  coef = coef_name,
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
) |>
  rownames_to_column("Probe") |>
  as_tibble()

write_result(tt, "ref2_primary_limma_all_cpgs.csv")
