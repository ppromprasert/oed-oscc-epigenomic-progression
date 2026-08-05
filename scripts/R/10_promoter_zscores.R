source("scripts/R/00_setup.R")
d <- read_csv(file.path(DATA,"cosmx_promoter_zscores.csv"),show_col_types=FALSE)
write_csv(d,file.path(OUT,"cosmx_promoter_zscores.csv"))
