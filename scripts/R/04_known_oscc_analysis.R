source("scripts/R/02_limma_ref2.R")
known <- c("CDKN2A","CDKN2B","CDH1","MGMT","PAX1","ZNF583","ATG5","MAP1LC3A","DAPK1","TP63","TFPI2","SOX17","GATA4","EGFR","PTK6","AIM2","CCNA1","CCNA2","CCNB1","CCNB2")
write_csv(res |> filter(gene %in% known) |> arrange(P.Value),
          file.path(OUT,"known_oscc_cpg_results.csv"))
