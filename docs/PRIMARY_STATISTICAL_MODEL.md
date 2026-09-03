# Primary Statistical Model

Differential methylation is modeled on M-values:

```text
M-value ~ Progression + Grade
```

Why:

- M-values provide suitable statistical behavior for methylation modeling.
- Grade adjustment reduces confounding by dysplasia severity.
- `limma` provides moderated CpG-level inference.
- `duplicateCorrelation` is used when repeated samples from the same patient are present.
- BH FDR controls the candidate-level multiple-testing burden.

Beta values are retained for direct methylation-fraction interpretation and delta-beta effect sizes.

FDR-significant CpGs are treated as the strongest evidence. Nominal P values and large delta-beta effects are hypothesis-generating, especially because Ref2 contains only three NonProgressors.
