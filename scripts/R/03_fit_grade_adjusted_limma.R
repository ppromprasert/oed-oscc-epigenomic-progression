source("scripts/R/01_build_baseline_cohorts.R")
M <- read_csv(file.path(DATA,"ref2_candidate_mvalues.csv"),show_col_types=FALSE) |> column_to_rownames("Probe")
M <- M[,ref2$matrix_sample_id,drop=FALSE]
design <- model.matrix(~ Progression + Grade,data=ref2)
if(any(table(ref2$patient_id)>1)){
 corfit <- duplicateCorrelation(M,design,block=ref2$patient_id)
 fit <- lmFit(M,design,block=ref2$patient_id,correlation=corfit$consensus)
}else fit <- lmFit(M,design)
fit <- eBayes(fit)
res <- topTable(fit,coef="ProgressionProgressor",number=Inf,sort.by="P") |> rownames_to_column("Probe")
write_csv(res,file.path(OUT,"limma_ref2_all_candidate_cpgs.csv"))
