source("scripts/R/04_beta_effect_sizes.R")
anno <- bind_rows(read_csv(file.path(DATA,"known_oscc_annotation.csv"),show_col_types=FALSE),
                  read_csv(file.path(DATA,"cosmx_annotation.csv"),show_col_types=FALSE))
out <- res |> left_join(summ,by="Probe") |> left_join(anno,by="Probe") |> group_by(Gene) |>
 summarise(n_probes=n(),n_nominal=sum(P.Value<.05),n_fdr=sum(adj.P.Val<.05),
 median_delta_beta=median(delta_beta),prop_abs_delta_ge_0_10=mean(abs_delta_beta>=.10),
 min_p=min(P.Value),min_fdr=min(adj.P.Val),.groups="drop")
write_csv(out,file.path(OUT,"gene_level_summary.csv"))
