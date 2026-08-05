# Primary Grade-Adjusted Differential Methylation Model

Statistical testing uses M-values.

```r
design <- model.matrix(~ Progression + Grade, data = metadata)
```

The primary coefficient is Progressor versus NonProgressor.

Grade is included to reduce confounding from Mild, Moderate, and Severe baseline dysplasia. When repeated patients are present, the project framework uses `limma::duplicateCorrelation` to model within-patient correlation.

BH-FDR is applied across tested CpGs.

Beta values are not used as the primary inferential scale; they are retained because they are directly interpretable as methylation fractions.
