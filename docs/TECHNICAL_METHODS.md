# Technical Methods

The research question is whether baseline OED methylation differs between lesions from patients who later progress to OSCC and those who do not.

## Real project architecture
- 92 QC-passed EPIC methylation samples.
- 705,131 CpGs retained after QC.
- M-values for limma inference; beta values for effect sizes and visualization.
- Baseline definitions: Ref1, Ref2, Ref3.
- Ref2 is the only exploratory inferential baseline because it contains 25 progressors and 3 non-progressors after grade filtering; Ref1 (22/1) and Ref3 (25/2) are descriptive.
- Primary candidate model: `M ~ Progression + Grade`.
- Repeated subjects are handled with `limma::duplicateCorrelation` when present.
- BH-FDR is applied across CpGs within the tested candidate set.

## Candidate layers
1. Literature-derived known OED/OSCC methylation genes.
2. CosMx-derived spatial transcriptomic genes.
3. Promoter CpG summaries and z-scores.
4. CpG-level beta-value effect sizes, including the clinically intuitive |delta-beta| >= 0.10 summary.

## Sensitivity
Because Ref2 has only three non-progressors, leave-one-non-progressor-out analyses are required for candidate stability. Paired precursor/later-lesion trajectories are complementary and require verified lesion/visit ordering.
