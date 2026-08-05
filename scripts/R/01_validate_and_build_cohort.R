source("scripts/R/00_setup.R")
meta <- read_csv(file.path(DATA,"ref2_metadata.csv"),show_col_types=FALSE)
beta <- read_csv(file.path(DATA,"ref2_beta_matrix.csv"),show_col_types=FALSE) |> column_to_rownames("probe")
M <- read_csv(file.path(DATA,"ref2_m_matrix.csv"),show_col_types=FALSE) |> column_to_rownames("probe")
stopifnot(identical(colnames(beta),colnames(M)))
stopifnot(setequal(colnames(beta),meta$sample_id))
meta <- meta |> slice(match(colnames(beta),sample_id)) |>
  mutate(progression=factor(progression,levels=c("NonProgressor","Progressor")),
         grade=factor(grade,levels=c("Mild","Moderate","Severe")))
write_csv(meta,file.path(OUT,"analysis_metadata.csv"))
