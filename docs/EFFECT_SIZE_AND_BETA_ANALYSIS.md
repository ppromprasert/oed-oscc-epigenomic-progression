# Beta-Scale Effect Sizes and Methylation Percentages

For each CpG, the expanded analysis calculates:
- progressor mean and median beta;
- non-progressor mean and median beta;
- mean delta beta;
- median delta beta;
- methylation percentage difference;
- Wilcoxon P/FDR;
- descriptive effect-size category.

Prespecified descriptive categories are:
- |delta-beta| >= 0.05: Modest
- |delta-beta| >= 0.10: Potentially meaningful
- |delta-beta| >= 0.20: Large

These are prioritization categories, **not universal clinical thresholds**.

A CpG is flagged for review when it is FDR-significant or has |delta-beta| >= 0.10. This explicitly separates statistical evidence from biological magnitude.
