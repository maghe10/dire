# Calculate overall mean absolute prediction error per isolate
# Run from the DIRE project root.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
})

source(file = "manuscript/manuscriptcommon.R")
source(file = "manuscript/metric_helpers_common.R")

# Same isolate-level error matrix used in errorPerGeneGroup.R
errors_tbl <- fetchPredictionErrors()

sample_mae <- errors_tbl %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample") %>%
  tidyr::pivot_longer(
    cols = -sample,
    names_to = "antibiotic",
    values_to = "error"
  ) %>%
  mutate(
    sample = as.character(sample),
    error = as.numeric(error),
    abs_error = abs(error)
  ) %>%
  group_by(sample) %>%
  summarise(
    sample_mean_abs_error = mean(abs_error, na.rm = TRUE),
    n_predictions = sum(!is.na(abs_error)),
    .groups = "drop"
  )

# Overall mean across isolates
# Each isolate contributes equally, matching errorPerGeneGroup.R.
overall_mean <- mean(sample_mae$sample_mean_abs_error, na.rm = TRUE)

cat("Number of isolates:", nrow(sample_mae), "\n")
cat("Overall mean absolute prediction error per isolate:\n")
cat(sprintf("  Raw proportion: %.8f\n", overall_mean))
cat(sprintf("  Percent: %.2f%%\n", 100 * overall_mean))
cat(sprintf("  Two significant figures: %.2g%%\n", 100 * overall_mean))

# Optional: inspect isolate-level values
print(sample_mae)
overall_mean <- mean(sample_mae$sample_mean_abs_error, na.rm = TRUE)
scales::percent(overall_mean, accuracy = 0.1)


