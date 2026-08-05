source("scripts/R/04_beta_effect_sizes.R")
write_csv(summ |> select(Probe,wilcox_p,wilcox_FDR,delta_beta,effect_category),
 file.path(OUT,"wilcoxon_beta_sensitivity.csv"))
