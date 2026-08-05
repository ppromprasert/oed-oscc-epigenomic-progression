# OED → OSCC Epigenomic Progression: Grade-Adjusted CpG Modeling, Effect Sizes, Spatial Integration, and Lesion Trajectories

**Expanded Illumina EPIC methylation case study of oral epithelial dysplasia (OED) progression to oral squamous cell carcinoma (OSCC), integrating baseline phenotype engineering, grade-adjusted limma, CpG/gene/region effect sizes, known cancer biology, spatial-transcriptomic candidates, robustness analyses, and within-patient lesion trajectories.**

## Scientific question

Can methylation in baseline OED distinguish lesions from patients who later progress to OSCC from lesions in non-progressors, and which candidate signals remain biologically meaningful when statistical significance, effect size, genomic context, and sensitivity to the very small non-progressor group are considered together?

## Real project scale

| Measure | Value |
|---|---:|
| QC-passed EPIC samples | **92** |
| CpGs retained after QC | **705,131** |
| Known OSCC candidate CpGs | **342** |
| CosMx candidate CpGs | **355** |
| CosMx promoter CpGs | **193** |
| CosMx spatial candidate genes | **26** |
| Overlap cases for promoter/spatial follow-up | **5** |

## Full analysis architecture

```text
EPIC beta + M matrices
          |
metadata/sample alignment
          |
Ref1 / Ref2 / Ref3 baseline engineering
          |
baseline feasibility assessment
          |
     Ref2: 25 P / 3 NP
          |
  M-value limma model
  Progression + Grade
          |
   +------+------+----------------+
   |             |                |
CpG inference  beta effects   genomic region
FDR/logFC     delta-beta      promoter/body/etc
   |             |                |
   +------+------+----------------+
          |
 known OSCC biology + CosMx/spatial biology
          |
   robustness / sensitivity
 Wilcoxon | Fisher | leave-one-NP-out
 diagnosis/race/sex/smoking exploratory tests
          |
 promoter z-scores + lesion trajectories
          |
 candidate prioritization / validation
```

## 1. Baseline phenotype engineering is part of the analysis

Three baseline definitions were evaluated because the selected precursor lesion changes the biological question.

| Definition | Total baseline | Modeled | Progressors | Non-progressors | Role |
|---|---:|---:|---:|---:|---|
| Ref1 | 29 | 23 | 22 | 1 | Descriptive |
| **Ref2** | **46** | **28** | **25** | **3** | **Exploratory inference** |
| Ref3 | 45 | 27 | 25 | 2 | Descriptive |

Ref2 is not called “good powered.” It is simply the only definition with enough non-progressors to support an exploratory model.

See [`docs/COHORT_AND_BASELINE_DEFINITIONS.md`](docs/COHORT_AND_BASELINE_DEFINITIONS.md).

## 2. Primary methylation model

M-values are used for statistical testing:

```r
design <- model.matrix(
  ~ Progression + Grade,
  data = ref2_metadata
)
```

The primary coefficient is Progressor versus NonProgressor.

**Grade is an explicit covariate** because Mild, Moderate, and Severe baseline dysplasia can differ in methylation independently of later progression.

When repeated patients occur, the analysis framework uses `limma::duplicateCorrelation`.

BH-FDR is applied across tested candidate CpGs.

Beta values are retained for biological interpretation and visualization.

## 3. Why M-values and beta values are both reported

The analysis separates:

**Statistical scale**
- M-values
- limma logFC
- P value
- BH-FDR

from:

**Biological effect-size scale**
- progressor mean beta
- non-progressor mean beta
- delta beta
- percentage-point methylation difference
- genomic region
- direction consistency

This prevents a small P value from being interpreted as a large methylation shift and prevents a large effect from being mislabeled statistically significant when the cohort is underpowered.

## 4. Descriptive delta-beta thresholds

The expanded Rmd prespecifies:

| Absolute delta beta | Label |
|---:|---|
| ≥ 0.05 | Modest |
| ≥ 0.10 | Potentially meaningful |
| ≥ 0.20 | Large |

These are **descriptive prioritization categories, not universal clinical thresholds**.

Candidates may be flagged for review if they are FDR-significant **or** have an absolute delta beta ≥ 0.10.

See [`docs/EFFECT_SIZE_AND_BETA_ANALYSIS.md`](docs/EFFECT_SIZE_AND_BETA_ANALYSIS.md).

## 5. Known OED/OSCC methylation biology

The literature-derived panel covers tumor suppressors, cell-cycle genes, epithelial state, DNA repair, autophagy, growth signaling, and other OED/OSCC candidates.

### Real gene-level examples

| Gene | CpGs | Nominal CpGs | FDR CpGs | Mean logFC | Min P | Direction |
|---|---:|---:|---:|---:|---:|---|
| CDH1 | 15 | 6 | 0 | -0.46 | 7.37e-4 | Lower |
| CDKN2A | 22 | 20 | 0 | -0.75 | 0.0329 | Lower |
| ATG5 | 12 | 12 | 0 | +0.07 | 0.0133 | Higher |

No known-OSCC CpG reached FDR < 0.05.

### Strongest known-OSCC CpG

**CDH1 cg01251360**
- limma logFC = **-3.61**
- P = **7.37 × 10^-4**
- FDR = **0.252**
- direction: lower methylation in progressors

This is a strong nominal signal, not an FDR-significant progression biomarker.

## 6. CDKN2A is analyzed as a region, not one cherry-picked CpG

CDKN2A remains biologically important despite no FDR-significant individual probe.

The expanded analysis examines:
- 22 retained CpGs;
- CpG-level beta distributions;
- delta-beta magnitude;
- promoter/body context;
- directional consistency;
- leave-one-nonprogressor-out stability;
- probe retention after QC;
- within-patient lesion trajectories.

This is exactly the distinction between **candidate prioritization** and a claim of statistical validation.

See [`docs/CDKN2A_DEEP_DIVE.md`](docs/CDKN2A_DEEP_DIVE.md).

## 7. Genomic-region analysis

CpGs are summarized separately for:

`TSS200`, `TSS1500`, `First exon`, `5'UTR`, `Gene body`, `3'UTR`, `Other`.

For each gene-region combination the workflow reports:
- number of probes;
- nominal/FDR probe counts;
- median delta beta;
- median methylation-percent difference;
- percentage of positive and negative limma effects;
- direction consistency;
- minimum P/FDR.

This is important because promoter methylation and gene-body methylation cannot be interpreted interchangeably.

## 8. CosMx / spatial-transcriptomic integration

The 26-gene CosMx candidate set contains strong immune/interferon and antigen-presentation biology:

`CXCL9`, `CXCL10`, `CXCL11`, `IFI6`, `STAT1`, `RNF213`, `MX1`, `XAF1`, `TAP1`, `GBP1`, `GBP5`, `OAS3`, `C1QA`, `C1QC`, `UBD`, and others.

Two different questions are asked:

1. **Full Ref2 cohort:** do CpGs mapped to spatial-expression candidate genes differ between progressors and non-progressors?
2. **Overlap cases:** how does promoter methylation for these genes vary relative to the full methylation cohort?

### Strongest adjusted finding: UBD

**UBD cg13206902**

| Metric | Result |
|---|---:|
| limma logFC | **1.96** |
| P | **1.32 × 10^-4** |
| FDR | **0.047** |
| Progressor mean methylation | **87.75%** |
| Non-progressor mean methylation | **66.21%** |
| Delta methylation | **+21.54 percentage points** |

This is the only CosMx candidate CpG reaching FDR < 0.05 in the exploratory Ref2 analysis.

The interpretation remains region-aware: higher methylation does **not** automatically imply lower gene expression.

## 9. Promoter methylation z-scores

Promoter CpGs are defined as:

- TSS1500
- TSS200
- 5'UTR
- first exon

For each gene and sample:

```text
retained promoter CpGs
       ↓
mean promoter beta
       ↓
standardize across all 92 methylation samples
       ↓
gene-specific promoter z-score
```

`z >= 1` is relatively high promoter methylation and `z <= -1` relatively low methylation compared with the methylation cohort.

The five spatial/methylation overlap cases show substantial patient-specific heterogeneity, particularly for priority interferon-related genes such as CXCL10 and IFI6.

## 10. Nonparametric sensitivity analysis

The expanded analysis does not rely on one modeling framework.

### Wilcoxon
Beta-value distributions are compared directly between progression groups.

This intentionally removes grade adjustment and asks whether the signal is visible without parametric limma assumptions.

### Categorical methylation / Fisher exact

Each CpG is standardized across the cohort:

```r
z <- scale(beta)
high_methylation <- z >= 1
```

Fisher exact testing then evaluates categorical high-methylation state versus progression.

These are sensitivity analyses; they do not replace the primary grade-adjusted model.

## 11. Leave-one-nonprogressor-out robustness

Ref2 has only three non-progressors.

The full limma model is therefore rerun three times:

```text
remove NP1 → Progression + Grade model
remove NP2 → Progression + Grade model
remove NP3 → Progression + Grade model
```

Each CpG receives:
- minimum/maximum leave-one-out logFC;
- median logFC;
- direction consistency across runs;
- minimum/maximum leave-one-out P value.

A candidate that changes direction after removing one comparator is less robust than a candidate retaining its effect direction across all runs.

## 12. Diagnosis, outcome, race, sex, and smoking analyses

The Rmd includes exploratory Kruskal-Wallis analyses for methylation variation across:

- current diagnosis;
- original clinical outcome;
- race;
- sex;
- smoking.

Tests are skipped when groups are too small, and BH-FDR is applied within each variable family.

These analyses are exploratory and are kept separate from the primary progression model.

## 13. Within-patient lesion trajectories

The analysis also asks a longitudinal question that is distinct from Ref2.

Diagnosis is ordered:

```text
Benign / no dysplasia
        ↓
Mild
        ↓
Moderate
        ↓
Severe
        ↓
CIS
        ↓
SCC
```

Priority genes include:

`CDKN2A`, `CDH1`, `UBD`, `CXCL10`, `IFI6`, `RNF213`, `MX1`, `XAF1`.

The earliest and latest valid lesions per patient are used to calculate gene-level beta change.

Current aggregate output contains only one evaluable pair for several priority genes, so those trajectory values are **descriptive** rather than inferential.

See [`docs/LESION_TRAJECTORIES.md`](docs/LESION_TRAJECTORIES.md).

## 14. Probe-retention and QC sensitivity

The current post-QC matrix identifies which candidate probes remain.

A complete before/after probe-loss analysis requires either:
- the original pre-QC matrix; or
- a probe-level QC manifest with removal reasons.

The expanded source analysis explicitly recommends against making 25–30% missingness a primary threshold and suggests prespecified 5%/10% thresholds when pre-QC data are available.

Targeted relaxed-QC probe rescue should remain a sensitivity analysis rather than redefining the primary result.

## 15. How candidates are actually prioritized

A serious candidate should be evaluated across multiple evidence dimensions:

```text
Formal FDR
    +
delta-beta magnitude
    +
methylation percentage
    +
genomic region
    +
direction consistency
    +
leave-one-out stability
    +
prior OED/OSCC biology
    +
spatial-expression evidence
    +
longitudinal evidence where valid
```

This is intentionally more rigorous than ranking candidates by nominal P value.

## Real conclusions

1. Ref2 is the only baseline definition supporting exploratory inference, but 25:3 remains severely underpowered.
2. Known OSCC genes contain biologically plausible progression-associated methylation patterns, but no CpG survives FDR.
3. CDH1 is the strongest nominal known-gene CpG signal.
4. CDKN2A shows broad exploratory progression-associated methylation structure warranting regional follow-up.
5. UBD cg13206902 is the strongest adjusted CosMx CpG result.
6. Immune/interferon-related promoter methylation is heterogeneous across overlapping spatial cases.
7. Effect sizes and leave-one-out sensitivity are essential before candidate prioritization.
8. Within-patient lesion trajectories are complementary but currently limited by pair availability and lesion-order verification.

## Repository structure

```text
config/
data/synthetic/
docs/
results/aggregate_real_results/
scripts/R/
scripts/python/
.github/workflows/
```

## Skills demonstrated

**Epigenomics:** Illumina EPIC, beta/M values, probe annotation, methylation QC.

**Statistical genomics:** limma, empirical Bayes, duplicateCorrelation, covariate adjustment, BH-FDR.

**Effect-size interpretation:** delta beta, methylation percentages, 5/10/20% descriptive thresholds.

**Tumor biology:** CDKN2A/CDKN2B, CDH1, MGMT, autophagy/cell-death and growth-signaling candidates.

**Spatial multi-omics:** CosMx-derived candidate genes, promoter methylation z-scores, interferon/immune program interpretation.

**Robustness:** Wilcoxon, Fisher exact, Kruskal-Wallis, leave-one-control-out sensitivity.

**Longitudinal analysis:** precursor/later-lesion ordering, paired gene methylation, lesion trajectories.

**Scientific judgment:** separation of formal significance, biological magnitude, sensitivity, and validation status.

## Data governance

No real patient-level identifiers, methylation matrices, lesion/block identifiers, dates, clinical records, or restricted spatial-transcriptomic data are distributed. Real results are limited to aggregate non-patient-level summaries. Synthetic data are used for public reproducibility.
