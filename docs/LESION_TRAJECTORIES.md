# Within-Patient Lesion Trajectories

The group-level Ref2 analysis asks whether baseline lesions differ between progressors and non-progressors.

The trajectory analysis asks a different question: does methylation change from an earlier precursor diagnosis to a later/higher-grade lesion within the same patient?

Diagnosis is ordered as:
Benign/no dysplasia -> Mild -> Moderate -> Severe -> CIS -> SCC.

For priority genes, mean gene beta is followed across lesions. When at least two valid patient pairs exist, paired Wilcoxon testing can be performed.

In the currently reported aggregate output, several priority genes have only one evaluable precursor/later pair, so the observed delta-beta values are descriptive and do not support inferential P values. Correct lesion and visit ordering must be verified before interpreting these as longitudinal evidence.
