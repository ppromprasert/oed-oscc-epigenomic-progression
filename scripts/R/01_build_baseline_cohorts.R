source("scripts/R/00_setup.R")
ref2 <- read_csv(file.path(DATA,"ref2_metadata.csv"),show_col_types=FALSE) |>
 mutate(Progression=factor(Progression,levels=c("NonProgressor","Progressor")),
        Grade=factor(Grade,levels=c("Mild","Moderate","Severe")))
stopifnot(nrow(ref2)==28,sum(ref2$Progression=="Progressor")==25,sum(ref2$Progression=="NonProgressor")==3)
write_csv(ref2,file.path(OUT,"ref2_analysis_metadata.csv"))
