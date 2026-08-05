source("scripts/R/00_setup.R")
d <- read_csv(file.path(DATA,"trajectory_gene_methylation.csv"),show_col_types=FALSE)
pairs <- d |> arrange(patient_id,Gene,diagnosis_order) |> group_by(patient_id,Gene) |>
 summarise(precursor_beta=first(gene_mean_beta),later_beta=last(gene_mean_beta),
 delta_beta=later_beta-precursor_beta,.groups="drop")
write_csv(pairs,file.path(OUT,"synthetic_precursor_later_pairs.csv"))
