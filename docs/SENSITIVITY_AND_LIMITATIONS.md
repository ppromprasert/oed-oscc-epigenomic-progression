# Sensitivity Analyses and Limitations

## Sensitivity analyses

The workflow uses:

- Ref1 / Ref2 / Ref3 baseline comparisons
- beta-value effect sizes
- nonparametric Wilcoxon tests
- large-effect delta-beta summaries
- categorical effect summaries where appropriate
- leave-one-NonProgressor-out robustness checks
- promoter / gene-context summaries
- specimen-level longitudinal trajectories

These analyses assess directional stability and biological plausibility. They do not replace the primary limma model.

## Limitations

- Ref2 is highly imbalanced: 25 Progressors vs 3 NonProgressors.
- Ref1 and Ref3 have even fewer NonProgressors.
- Candidate analyses remain vulnerable to small-group instability.
- Exposure signatures vary in CpG coverage.
- Single-CpG exposure signatures are especially exploratory.
- Exposure scores are not direct dose measures.
- Cross-platform methylation / spatial-expression associations do not establish regulation.
- The longitudinal Mild-to-Severe analysis includes only three patients.
- No clean methylation-profiled dysplasia-to-SCC specimen pair was available.

The strongest formal result is UBD cg13206902. Most other findings are best treated as prioritization signals for validation.
