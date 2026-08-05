# Cohort and Baseline Definitions

The real dataset contains 92 QC-passed methylation samples.

Three baseline definitions were evaluated:
- Ref1: lowest grade before the exact lesion;
- Ref2: lowest grade ever for the exact lesion;
- Ref3: lowest diagnosis before the lesion across patient lesions.

Modeled cohorts:
- Ref1: 22 progressors / 1 non-progressor — descriptive only.
- Ref2: 25 progressors / 3 non-progressors — exploratory inference.
- Ref3: 25 progressors / 2 non-progressors — descriptive only.

Each baseline definition changes the biological estimand. Ref2 is used for inference because it is the only definition with at least three non-progressors, but it remains severely imbalanced.
