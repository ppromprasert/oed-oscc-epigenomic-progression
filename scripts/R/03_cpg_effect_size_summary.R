source("scripts/R/02_limma_ref2.R")
summary <- res |>
 group_by(gene) |>
 summarise(n_cpg=n(), n_abs_delta_beta_ge_0_10=sum(abs(delta_beta)>=0.10),
           proportion_ge_0_10=n_abs_delta_beta_ge_0_10/n_cpg,
           median_delta_beta=median(delta_beta), min_p=min(P.Value),
           min_fdr=min(adj.P.Val),.groups="drop")
write_csv(summary,file.path(OUT,"gene_10pct_effect_summary.csv"))
