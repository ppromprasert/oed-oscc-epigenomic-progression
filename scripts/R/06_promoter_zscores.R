source("scripts/R/01_validate_and_build_cohort.R")
anno <- read_csv(file.path(DATA,"candidate_annotation.csv"),show_col_types=FALSE)
prom <- anno |> filter(promoter)
long <- beta[intersect(prom$probe,rownames(beta)),,drop=FALSE] |>
 as.data.frame() |> rownames_to_column("probe") |>
 pivot_longer(-probe,names_to="sample_id",values_to="beta") |>
 left_join(prom,by="probe")
gene_sample <- long |> group_by(gene,sample_id) |> summarise(promoter_beta=mean(beta),.groups="drop") |>
 group_by(gene) |> mutate(z=(promoter_beta-mean(promoter_beta))/sd(promoter_beta)) |> ungroup()
write_csv(gene_sample,file.path(OUT,"promoter_gene_zscores.csv"))
