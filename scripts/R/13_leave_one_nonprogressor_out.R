source("scripts/R/01_build_baseline_cohorts.R")
M <- read_csv(file.path(DATA,"ref2_candidate_mvalues.csv"),show_col_types=FALSE) |> column_to_rownames("Probe")
np <- ref2$matrix_sample_id[ref2$Progression=="NonProgressor"]
long <- map_dfr(np,function(drop_id){
 d <- ref2 |> filter(matrix_sample_id!=drop_id) |> droplevels()
 design <- model.matrix(~ Progression + Grade,data=d)
 fit <- eBayes(lmFit(M[,d$matrix_sample_id,drop=FALSE],design))
 topTable(fit,coef="ProgressionProgressor",number=Inf,sort.by="none") |> rownames_to_column("Probe") |>
  transmute(dropped_nonprogressor=drop_id,Probe,logFC,p_value=P.Value)
})
summary <- long |> group_by(Probe) |> summarise(n_runs=n(),min_logFC=min(logFC),max_logFC=max(logFC),
 median_logFC=median(logFC),same_direction_all_runs=all(logFC>0)|all(logFC<0),
 min_p=min(p_value),max_p=max(p_value),.groups="drop")
write_csv(long,file.path(OUT,"leave_one_out_long.csv"));write_csv(summary,file.path(OUT,"leave_one_out_summary.csv"))
