source("scripts/R/01_build_baseline_cohorts.R")
beta <- read_csv(file.path(DATA,"ref2_candidate_beta.csv"),show_col_types=FALSE) |> column_to_rownames("Probe")
out <- map_dfr(rownames(beta),function(cpg){
 d <- tibble(matrix_sample_id=colnames(beta),beta=as.numeric(beta[cpg,])) |> mutate(z=as.numeric(scale(beta)),high=z>=1) |> left_join(ref2,by="matrix_sample_id")
 tab <- table(d$high,d$Progression)
 if(!all(dim(tab)==c(2,2))) return(tibble(Probe=cpg,odds_ratio=NA_real_,p_value=NA_real_))
 f <- fisher.test(tab); tibble(Probe=cpg,odds_ratio=unname(f$estimate),p_value=f$p.value)
}) |> mutate(FDR=p.adjust(p_value,"BH"))
write_csv(out,file.path(OUT,"fisher_high_methylation_sensitivity.csv"))
