# Sensitivity and Robustness Analyses

## Beta-scale Wilcoxon
CpG beta distributions are compared without grade adjustment. These are complementary to, not replacements for, the primary limma model.

## Fisher categorical test
Each CpG is standardized across samples and classified as high methylation at z >= 1. Fisher exact tests compare the categorical state with progression.

## Leave-one-nonprogressor-out
Ref2 contains only three non-progressors. The grade-adjusted limma model is rerun after removing each non-progressor separately. Stability summaries record direction consistency and P-value ranges.

## Diagnosis/demographic nonparametric analyses
Where cell sizes permit, CpG methylation is explored across:
- current diagnosis;
- original outcome;
- race;
- sex;
- smoking.

Kruskal-Wallis tests use BH-FDR within each variable family.

The sensitivity analyses assess robustness; they do not manufacture significance when the primary FDR screen is negative.
