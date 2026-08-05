suppressPackageStartupMessages({library(tidyverse);library(limma);library(pheatmap)})
DATA <- "data/synthetic"
OUT <- "results/generated"
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
