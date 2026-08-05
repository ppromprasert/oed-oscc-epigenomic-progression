source("scripts/R/05_gene_summary.R")
known_genes <- jsonlite::fromJSON("config/gene_sets.json")$KNOWN_OSCC
write_csv(out |> filter(Gene %in% known_genes),file.path(OUT,"known_oscc_gene_summary.csv"))
