# Specimen-level longitudinal lesion audit.
#
# Chronology is defined from specimen pathology and specimen date.
# Repeated same-date samples are summarized within patient/grade to avoid
# pseudoreplication. A Mild -> Severe trajectory is counted only when the
# Severe specimen occurs after the Mild specimen.

source("scripts/R/00_setup.R")

score_path <- file.path(TABLE_DIR, "exposure_signature_scores_long.csv")
meta_path <- file.path(DATA_DIR, "public", "all_sample_metadata.csv")

scores <- read_csv(score_path, show_col_types = FALSE)
meta <- read_csv(meta_path, show_col_types = FALSE)

needed <- c("sample_id", "patient_id", "pathology_grade", "specimen_date")
missing <- setdiff(needed, names(meta))
if (length(missing) > 0) {
  stop("Missing longitudinal metadata columns: ", paste(missing, collapse = ", "))
}

dat <- scores |>
  inner_join(meta, by = "sample_id") |>
  mutate(specimen_date = as.Date(specimen_date))

# Average multiple samples collected on the same date within patient, grade,
# exposure, and signature.
collapsed <- dat |>
  group_by(patient_id, specimen_date, pathology_grade, Exposure, Signature) |>
  summarise(score = mean(score, na.rm = TRUE), .groups = "drop")

# Select a chronologically valid Mild -> Severe pair for each patient/signature.
# The earliest Mild specimen is used as baseline, followed by the earliest
# Severe specimen occurring after that baseline.
paired <- collapsed |>
  filter(pathology_grade %in% c("Mild", "Severe")) |>
  group_by(patient_id, Exposure, Signature) |>
  group_modify(function(df, key) {
    mild <- df |>
      filter(pathology_grade == "Mild") |>
      arrange(specimen_date)

    if (nrow(mild) == 0) return(tibble())

    mild_row <- mild |> slice(1)

    severe <- df |>
      filter(
        pathology_grade == "Severe",
        specimen_date > mild_row$specimen_date
      ) |>
      arrange(specimen_date)

    if (nrow(severe) == 0) return(tibble())

    severe_row <- severe |> slice(1)

    tibble(
      mild_date = mild_row$specimen_date,
      severe_date = severe_row$specimen_date,
      mild_score = mild_row$score,
      severe_score = severe_row$score,
      delta_severe_minus_mild = severe_row$score - mild_row$score
    )
  }) |>
  ungroup()

write_result(paired, "mild_to_severe_exposure_trajectories.csv")
