# Analysis Architecture

This project evaluates whether baseline oral epithelial dysplasia (OED) methylation differs between patients who later progress to oral squamous cell carcinoma (OSCC) and non-progressors.

The analysis has six linked layers:
1. baseline/phenotype engineering;
2. grade-adjusted CpG-level inference;
3. beta-scale biological effect sizes and genomic-region context;
4. known OSCC and CosMx/spatial candidate biology;
5. robustness/sensitivity analysis;
6. within-patient lesion trajectories.

The project intentionally distinguishes formal FDR evidence from biological prioritization based on delta-beta magnitude, direction consistency, and sensitivity analyses.
