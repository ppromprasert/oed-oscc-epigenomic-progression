# Analysis Architecture

The project is organized around four analytical layers:

1. **Cohort engineering**
   Resolve sample IDs, lesion chronology, and Ref1/Ref2/Ref3 baseline definitions.

2. **Progression-associated methylation**
   Fit a grade-adjusted limma model on M-values and summarize beta-value effect sizes.

3. **Biological integration**
   Evaluate known OSCC candidates, CDKN2A, spatial/CosMx genes, promoter context, and exposure-associated methylation signatures.

4. **Robustness and longitudinal interpretation**
   Use nonparametric tests, leave-one-NonProgressor-out checks, large-effect summaries, and specimen-level Mild-to-Severe trajectories.

The primary inferential cohort is Ref2 because it contains 25 Progressors and 3 NonProgressors, whereas Ref1 and Ref3 have fewer NonProgressors.
