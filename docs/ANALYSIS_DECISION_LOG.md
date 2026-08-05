# Analysis Decision Log

| Decision | Reason |
|---|---|
| Use M-values for limma | Better statistical properties for array methylation inference |
| Use beta values for plots and delta-beta | Directly interpretable methylation proportion |
| Adjust for dysplasia grade | Baseline grade may confound progression-associated methylation |
| Ref2 only for inference | Ref1/Ref3 non-progressor groups are too small |
| BH-FDR | CpG-level multiplicity |
| Report delta-beta alongside P/FDR | Statistical significance and biological magnitude are distinct |
| 10% absolute methylation threshold | Effect-size summary, not a significance threshold |
| Leave-one-NP-out sensitivity | Ref2 inference is highly sensitive to three non-progressors |
| Separate CpG, gene-region, and promoter summaries | Methylation direction depends on genomic context |
| Do not infer expression directly from methylation | CpG location and regulatory context matter |
