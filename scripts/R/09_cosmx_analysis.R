source("scripts/R/05_gene_summary.R")
cosmx_genes <- jsonlite::fromJSON("config/gene_sets.json")$COSMX
write_csv(out |> filter(Gene %in% cosmx_genes),file.path(OUT,"cosmx_gene_summary.csv"))
