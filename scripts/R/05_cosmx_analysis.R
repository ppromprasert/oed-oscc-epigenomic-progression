source("scripts/R/02_limma_ref2.R")
cosmx <- c("CXCL9","CXCL10","CXCL11","IFI6","WARS1","GBP5","CD74","STAT1","IFI44","UBE2L6","RNF213","OAS3","IFITM3","BST2","MX1","PARP14","GBP1","EPSTI1","LY6E","LCP1","UBD","XAF1","C1QC","TAP1","LAPTM5","C1QA")
write_csv(res |> filter(gene %in% cosmx) |> arrange(P.Value),
          file.path(OUT,"cosmx_cpg_results.csv"))
