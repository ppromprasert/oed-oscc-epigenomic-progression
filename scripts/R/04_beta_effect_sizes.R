source("scripts/R/03_fit_grade_adjusted_limma.R")
beta <- read_csv(file.path(DATA,"ref2_candidate_beta.csv"),show_col_types=FALSE) |> column_to_rownames("Probe")
summ <- map_dfr(rownames(beta),function(cpg){
 d <- tibble(matrix_sample_id=colnames(beta),beta=as.numeric(beta[cpg,])) |> left_join(ref2,by="matrix_sample_id")
 a <- d$beta[d$Progression=="Progressor"]; b <- d$beta[d$Progression=="NonProgressor"]
 wt <- wilcox.test(a,b,exact=FALSE)
 tibble(Probe=cpg,progressor_mean=mean(a),nonprogressor_mean=mean(b),
        delta_beta=mean(a)-mean(b),wilcox_p=wt$p.value)
}) |> mutate(abs_delta_beta=abs(delta_beta),
 effect_category=case_when(abs_delta_beta>=.20~"Large",abs_delta_beta>=.10~"Potentially meaningful",abs_delta_beta>=.05~"Modest",TRUE~"Small"),
 wilcox_FDR=p.adjust(wilcox_p,"BH"))
write_csv(summ,file.path(OUT,"cpg_beta_effect_sizes.csv"))
