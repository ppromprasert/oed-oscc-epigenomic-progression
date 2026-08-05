# Analysis Decision Log

| Decision | Rationale |
|---|---|
| Keep Ref1/Ref2/Ref3 separate | Each baseline definition asks a different biological question |
| Use Ref2 for inference | Only definition with at least 3 non-progressors |
| Adjust for grade | Baseline dysplasia grade may confound progression-associated methylation |
| Test on M-values | Better statistical behavior for methylation inference |
| Interpret on beta scale | Direct methylation-fraction/effect-size interpretation |
| Use BH-FDR | CpG-level multiplicity |
| Report 5/10/20% delta-beta categories | Biological prioritization, not significance thresholds |
| Exact gene-symbol probe mapping | Prevent false mapping to related gene names |
| Summarize genomic regions | Promoter and gene-body methylation have different regulatory meaning |
| Run Wilcoxon/Fisher sensitivities | Test robustness to model/scale choices |
| Leave one non-progressor out | Quantify sensitivity to each of only three controls |
| Use tumor/patient, not CpGs, as biological context | Avoid pseudo-replication in interpretation |
| Treat lesion trajectories separately | Longitudinal change and baseline group difference are different questions |
| Do not overinterpret one paired lesion | n=1 trajectories are descriptive only |
| Keep relaxed QC as sensitivity only | Prevent targeted probe rescue from redefining primary evidence |
