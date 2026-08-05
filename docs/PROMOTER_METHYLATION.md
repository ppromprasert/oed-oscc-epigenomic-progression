# Promoter Methylation and Z-Score Integration

Promoter CpGs are defined using TSS1500, TSS200, 5'UTR, and first-exon annotations.

For each gene/sample:
1. select retained promoter CpGs;
2. calculate mean promoter beta;
3. standardize the gene across the 92-sample methylation cohort;
4. interpret overlapping samples relative to that cohort distribution.

A z-score >= 1 indicates relatively high promoter methylation and <= -1 relatively low promoter methylation versus the cohort.

These are relative methylation states. The repository does not assume that higher methylation always means lower expression; regulatory meaning depends on genomic context.
