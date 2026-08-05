source("scripts/R/01_validate_and_build_cohort.R")
np_ids <- meta$sample_id[meta$progression=="NonProgressor"]
all_runs <- list()
for(drop_id in np_ids){
 keep <- meta$sample_id != drop_id
 dat <- droplevels(meta[keep,])
 M2 <- M[,keep,drop=FALSE]
 design <- model.matrix(~ progression + grade,data=dat)
 fit <- eBayes(lmFit(M2,design))
 coef_name <- grep("^progression",colnames(design),value=TRUE)[1]
 tt <- topTable(fit,coef=coef_name,number=Inf,sort.by="none") |>
   rownames_to_column("probe") |> mutate(dropped_nonprogressor=drop_id)
 all_runs[[drop_id]] <- tt
}
write_csv(bind_rows(all_runs),file.path(OUT,"leave_one_nonprogressor_out.csv"))
