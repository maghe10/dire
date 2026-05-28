source(file = "common.R")
source(file = "genotype/genotypeCommon.R")
source(file = "manuscript/manuscriptcommon.R")
source(file = "manuscript/metric_helpers_common.R")
source("manuscript/generic_plot_helpers.R")

suppressPackageStartupMessages({
  library(tidyverse)
})


outdir <- file.path(manuscriptDirectory, "error_vs_genotype")
outdirPlot <- file.path(outdir, "plot")
dir.create(outdirPlot, recursive = TRUE, showWarnings = FALSE)




clean_token <- function(x) {
  x %>%
    stringr::str_trim() %>%
    na_if("") %>%
    na_if("-") %>%
    na_if("NA")
}

beta_property <- function(x) {
  dplyr::case_when(
    is.na(x) ~ NA_character_,
    stringr::str_detect(x, "^3") ~ "ESBL carba",
    x == "1" ~ "AmpC",
    x == "2be" ~ "ESBL classic",
    x %in% c("2b", "2br", "2d") ~ "Non-ESBL",
    TRUE ~ NA_character_
  )
}

quinolone_property <- function(x) {
  dplyr::case_when(
    is.na(x) ~ NA_character_,
    stringr::str_detect(x, "QRDR") ~ "QRDR",
    x %in% c("QnrB", "QnrS", "QepA") ~ "PMQR",
    TRUE ~ NA_character_
  )
}

aminoglycoside_property <- function(x) {
  dplyr::case_when(
    is.na(x) ~ NA_character_,
    x %in% c("AAC(3)", "AAC(6)") ~ "AAC",
    stringr::str_detect(x, "^APH") ~ "APH",
    x %in% c("AAD", "ANT") ~ "AAD/ANT",
    TRUE ~ NA_character_
  )
}

make_property_table <- function(genotype_tbl, col, summariser) {
  genotype_tbl %>%
    transmute(
      sample = as.character(sample),
      raw = .data[[col]]
    ) %>%
    tidyr::separate_rows(raw, sep = ";") %>%
    mutate(
      raw = clean_token(raw),
      property = summariser(raw)
    ) %>%
    filter(!is.na(property)) %>%
    distinct(sample, property)
}

# ============================================================
# Forest plot: mean MAE difference with confidence intervals
# ============================================================

plot_ttest_forest <- function(ttest_tbl,
                              p_adjust_col = "p_adjust_fdr",
                              diff_col = "difference_present_minus_absent",
                              conf_low_col = "conf_low",
                              conf_high_col = "conf_high",
                              family_col = "property",
                              significance_cutoff = 0.05,
                              base_size = 11) {
  
  plot_tbl <- ttest_tbl %>%
    mutate(
      significant = !is.na(.data[[p_adjust_col]]) &
        .data[[p_adjust_col]] < significance_cutoff,
      has_ci =
        !is.na(.data[[conf_low_col]]) &
        !is.na(.data[[conf_high_col]]),
      has_difference = !is.na(.data[[diff_col]])
    ) %>%
    filter(has_difference) %>%
    arrange(.data[[diff_col]]) %>%
    mutate(
      family_ordered = factor(
        .data[[family_col]],
        levels = unique(.data[[family_col]])
      )
    )
  
  ci_tbl <- plot_tbl %>%
    filter(has_ci)
  
  no_ci_tbl <- plot_tbl %>%
    filter(!has_ci)
  
  ggplot(
    plot_tbl,
    aes(
      x = .data[[diff_col]],
      y = family_ordered
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey50"
    ) +
    geom_errorbarh(
      data = ci_tbl,
      aes(
        xmin = .data[[conf_low_col]],
        xmax = .data[[conf_high_col]],
        color = significant
      ),
      height = 0.2,
      linewidth = 0.8
    ) +
    geom_point(
      data = ci_tbl,
      aes(color = significant),
      size = 2.8
    ) +
    geom_point(
      data = no_ci_tbl,
      aes(color = significant),
      size = 2.8,
      shape = 1,
      stroke = 1
    ) +
    scale_color_manual(
      values = c(
        "TRUE" = "red",
        "FALSE" = "black"
      ),
      guide = "none"
    ) +
    labs(
      x = "Difference in mean absolute error\n(Present − Absent)",
      y = NULL
    ) +
    theme_manuscript(base_size = base_size)
}


ALL <- function()
{
  errors_tbl <- fetchPredictionErrors()
  genotypeGroupTable <- getGenotypeGroupTable()
  
  
  sample_mae <- errors_tbl %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample") %>%
    pivot_longer(
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
      sample_median_abs_error = median(abs_error, na.rm = TRUE),
      n_predictions = sum(!is.na(abs_error)),
      .groups = "drop"
    )
  
  analysis_samples <- intersect(
    sample_mae$sample,
    as.character(genotypeGroupTable$sample)
  )
  
  sample_mae <- sample_mae %>%
    filter(sample %in% analysis_samples)
  
  n_total <- nrow(sample_mae)
  overall_mean <- mean(sample_mae$sample_mean_abs_error, na.rm = TRUE)
  overall_sum <- sum(sample_mae$sample_mean_abs_error, na.rm = TRUE)
  
  property_sample_long <- bind_rows(
    make_property_table(genotypeGroupTable, "Functional_groups", beta_property),
    make_property_table(genotypeGroupTable, "Quinolone_families", quinolone_property),
    make_property_table(genotypeGroupTable, "Aminoglycoside_families", aminoglycoside_property)
  ) %>%
    filter(sample %in% analysis_samples) %>%
    distinct(sample, property)
  
  all_properties <- property_sample_long %>%
    distinct(property) %>%
    arrange(property) %>%
    pull(property)
  
  mae_property_long <- tidyr::crossing(
    sample = sample_mae$sample,
    property = all_properties
  ) %>%
    left_join(
      property_sample_long %>% mutate(present = TRUE),
      by = c("sample", "property")
    ) %>%
    mutate(
      present = replace_na(present, FALSE),
      status = if_else(present, "Present", "Absent")
    ) %>%
    left_join(sample_mae, by = "sample")
  
  mae_property_summary <- mae_property_long %>%
    group_by(property, status, present) %>%
    summarise(
      n_samples = n(),
      mean_abs_error = mean(sample_mean_abs_error, na.rm = TRUE),
      median_abs_error = median(sample_mean_abs_error, na.rm = TRUE),
      sd_abs_error = sd(sample_mean_abs_error, na.rm = TRUE),
      min_abs_error = min(sample_mean_abs_error, na.rm = TRUE),
      max_abs_error = max(sample_mean_abs_error, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(property, desc(present))
  
  mae_property_wide <- mae_property_summary %>%
    select(
      property,
      status,
      n_samples,
      mean_abs_error,
      median_abs_error,
      sd_abs_error,
      min_abs_error,
      max_abs_error
    ) %>%
    pivot_wider(
      names_from = status,
      values_from = c(
        n_samples,
        mean_abs_error,
        median_abs_error,
        sd_abs_error,
        min_abs_error,
        max_abs_error
      )
    ) %>%
    mutate(
      mean_difference_present_minus_absent =
        mean_abs_error_Present - mean_abs_error_Absent
    )
  
  ttest_results <- mae_property_long %>%
    group_by(property) %>%
    group_modify(function(.x, .y) {
      present_values <- .x %>%
        filter(present) %>%
        pull(sample_mean_abs_error) %>%
        na.omit()
      
      absent_values <- .x %>%
        filter(!present) %>%
        pull(sample_mean_abs_error) %>%
        na.omit()
      
      n_present <- length(present_values)
      n_absent <- length(absent_values)
      
      mean_present <- mean(present_values, na.rm = TRUE)
      mean_absent <- mean(absent_values, na.rm = TRUE)
      
      if (n_present < 2 || n_absent < 2) {
        return(tibble(
          n_present = n_present,
          n_absent = n_absent,
          mean_present = mean_present,
          mean_absent = mean_absent,
          difference_present_minus_absent = mean_present - mean_absent,
          p_value = NA_real_,
          conf_low = NA_real_,
          conf_high = NA_real_,
          test_performed = FALSE,
          reason = paste0("Too few observations: present=", n_present, ", absent=", n_absent)
        ))
      }
      
      tt <- t.test(present_values, absent_values)
      
      tibble(
        n_present = n_present,
        n_absent = n_absent,
        mean_present = mean_present,
        mean_absent = mean_absent,
        difference_present_minus_absent = mean_present - mean_absent,
        p_value = tt$p.value,
        conf_low = tt$conf.int[1],
        conf_high = tt$conf.int[2],
        test_performed = TRUE,
        reason = NA_character_
      )
    }) %>%
    ungroup() %>%
    mutate(
      p_adjust_fdr = p.adjust(p_value, method = "fdr")
    )
  
  sanity_check <- mae_property_wide %>%
    transmute(
      property,
      n_total_expected = n_total,
      n_total_observed = n_samples_Present + n_samples_Absent,
      overall_mean = overall_mean,
      overall_sum = overall_sum,
      weighted_sum =
        n_samples_Present * mean_abs_error_Present +
        n_samples_Absent * mean_abs_error_Absent,
      weighted_mean = weighted_sum / n_total_observed,
      difference_from_overall_sum = weighted_sum - overall_sum,
      difference_from_overall_mean = weighted_mean - overall_mean,
      passed =
        n_total_observed == n_total_expected &
        abs(difference_from_overall_sum) < 1e-10
    )
  
  readr::write_csv2(sample_mae, file.path(outdir, "sample_mae.csv"))
  readr::write_csv2(property_sample_long, file.path(outdir, "property_per_sample_long.csv"))
  readr::write_csv2(mae_property_long, file.path(outdir, "sample_mae_by_property_long.csv"))
  readr::write_csv2(mae_property_summary, file.path(outdir, "summary_present_absent_long.csv"))
  readr::write_csv2(mae_property_wide, file.path(outdir, "summary_present_absent_wide.csv"))
  readr::write_csv2(ttest_results, file.path(outdir, "welch_ttests_present_absent.csv"))
  readr::write_csv2(sanity_check, file.path(outdir, "sanity_check_weighted_means.csv"))
  
  if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(
      list(
        sample_mae = sample_mae,
        property_per_sample = property_sample_long,
        sample_by_property = mae_property_long,
        summary_long = mae_property_summary,
        summary_wide = mae_property_wide,
        t_tests = ttest_results,
        sanity_check = sanity_check
      ),
      path = file.path(outdir, "mean_abs_error_by_property.xlsx")
    )
  }
  
  p_box <- ggplot(
    mae_property_long,
    aes(
      x = status,
      y = sample_mean_abs_error
    )
  ) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.12, alpha = 0.55, size = 1.5) +
    facet_wrap(~ property, scales = "free_x") +
    labs(
      x = NULL,
      y = "Mean absolute prediction error"
    ) +
    theme_manuscript(base_size = 11) +
    theme(
      strip.text = element_text(face = "bold")
    )
  
  p_mean <- ggplot(
    mae_property_summary,
    aes(
      x = status,
      y = mean_abs_error
    )
  ) +
    geom_col(width = 0.7) +
    facet_wrap(~ property, scales = "free_x") +
    labs(
      x = NULL,
      y = "Mean absolute prediction error"
    ) +
    theme_manuscript(base_size = 11) +
    theme(
      strip.text = element_text(face = "bold")
    )
  
  ggsave(
    file.path(outdirPlot, "boxplot_mean_abs_error_by_property.png"),
    p_box,
    width = 12,
    height = 7,
    dpi = 300
  )
  
  ggsave(
    file.path(outdirPlot, "barplot_mean_abs_error_by_property.png"),
    p_mean,
    width = 12,
    height = 7,
    dpi = 300
  )
  
  print(mae_property_wide)
  print(ttest_results)
  print(sanity_check)
  
  if (any(!sanity_check$passed)) {
    warning("At least one property failed the weighted-mean sanity check.")
  }
  
  p_forest <- plot_ttest_forest(
    ttest_results
  )
  
  ggsave(
    file.path(outdirPlot, "ttest_difference_forest_plot.png"),
    p_forest,
    width = 8,
    height = 5,
    dpi = 300
  )
  
  
  print(p_forest)
  
  print(p_box)
  print(p_mean)
}
