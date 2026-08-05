suppressPackageStartupMessages({library(tidyverse);library(limma)})
DATA <- "data/synthetic"
OUT <- "results/generated"
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
