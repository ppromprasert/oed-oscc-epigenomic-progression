# EPIC QC and Probe Mapping

The final analysis uses 92 QC-passed samples and 705,131 retained CpGs.

Candidate mapping uses exact gene-symbol matching. This matters because partial matching can create false probe assignments (for example, CDH1 versus PCDH-family genes or AIM2 versus FAIM2).

The clean report contains:
- 342 known-OSCC candidate CpGs;
- 355 CosMx candidate CpGs;
- 193 CosMx promoter CpGs.

The expanded workflow also tracks priority-gene probe retention after QC. A full before/after removal audit requires the pre-QC matrix or a probe-level QC manifest; relaxed-QC rescue analyses should remain separate sensitivity analyses.
