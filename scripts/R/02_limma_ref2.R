source("scripts/R/01_validate_and_build_cohort.R")
design <- model.matrix(~ progression + grade,data=meta)
patient_tab <- table(meta$patient_id)
if(any(patient_tab>1)){
  corfit <- duplicateCorrelation(M,design,block=meta$patient_id)
  fit <- lmFit(M,design,block=meta$patient_id,correlation=corfit$consensus)
}else{
  fit <- lmFit(M,design)
}
fit <- eBayes(fit)
coef_name <- grep("^progression",colnames(design),value=TRUE)[1]
res <- topTable(fit,coef=coef_name,number=Inf,sort.by="P") |> rownames_to_column("probe")
means <- tibble(probe=rownames(beta),
                beta_progressor=rowMeans(beta[,meta$progression=="Progressor",drop=FALSE]),
                beta_nonprogressor=rowMeans(beta[,meta$progression=="NonProgressor",drop=FALSE])) |>
  mutate(delta_beta=beta_progressor-beta_nonprogressor)
anno <- read_csv(file.path(DATA,"candidate_annotation.csv"),show_col_types=FALSE)
res <- res |> left_join(means,by="probe") |> left_join(anno,by="probe")
write_csv(res,file.path(OUT,"ref2_limma_candidate_results.csv"))
