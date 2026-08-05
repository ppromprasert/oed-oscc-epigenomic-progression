source("scripts/R/04_beta_effect_sizes.R")
anno <- read_csv(file.path(DATA,"known_oscc_annotation.csv"),show_col_types=FALSE)
cdk <- res |> left_join(summ,by="Probe") |> left_join(anno,by="Probe") |> filter(Gene=="CDKN2A")
write_csv(cdk,file.path(OUT,"CDKN2A_all_cpg_results.csv"))
