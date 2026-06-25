# ============================================================================
# model_genotype.R
# ----------------------------------------------------------------------------
# Purpose:
# Analyse whether AI prediction performance is associated with genotype groups,
# both overall and per antibiotic.
#
# Inputs:
#   1) readCountSampleFrameLong() available in the environment
#   2) sample_genotype_grouping_simple.csv
#
# Main outputs:
#   - sample-level performance tables
#   - sample-antibiotic-level performance tables
#   - merged analysis datasets
#   - sample-level models
#   - antibiotic-specific models
#   - plots
#   - heatmaps
#   - prediction-pattern clustering
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(forcats)
  library(broom)
  library(scales)
})

# ----------------------------------------------------------------------------
# CONFIG
# ----------------------------------------------------------------------------

source(file = "manuscript/manuscriptcommon.R")

indir <- file.path(processedRootRassembly, "genotype")

OUTDIR <- file.path(processedRootRassembly, "prediction_vs_genotype_output")
GENOTYPE_GROUP_FILE <- file.path(indir, "sample_genotype_grouping_simple.csv")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTDIR, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTDIR, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTDIR, "models"), recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------------
# HELPERS
# ----------------------------------------------------------------------------

message2 <- function(...) {
  cat(paste0(..., "\n"))
}

resolve_col <- function(df, candidates, required = TRUE) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) {
    if (required) {
      stop("Could not find any of these columns: ",
           paste(candidates, collapse = ", "))
    }
    return(NULL)
  }
  hit[[1]]
}

safe_ratio <- function(num, den) {
  ifelse(is.na(den) | den == 0, NA_real_, num / den)
}

save_csv <- function(x, filename) {
  readr::write_csv(x, file.path(OUTDIR, "tables", filename))
}

save_plot <- function(plot_obj, filename, width = 7, height = 5) {
  ggsave(
    filename = file.path(OUTDIR, "plots", filename),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 600
  )
}

sample_order_key <- function(x) {
  suppressWarnings(as.integer(normalize_sample_id(x)))
}

# ----------------------------------------------------------------------------
# STEP 1. READ AI PREDICTION COUNTS
# ----------------------------------------------------------------------------

read_prediction_counts <- function() {
  if (!exists("readCountSampleFrameLong", mode = "function")) {
    stop("Function readCountSampleFrameLong() is not available in the environment.")
  }
  
  x <- readCountSampleFrameLong() %>%
    filter_and_drop(cpmode, "normal") %>%
    filter_and_drop(mode, MODE_A) %>%
    filter_and_drop(significanceLevel, NA)
  
  sample_col <- resolve_col(x, c("sample", "sample_id"))
  ab_col <- resolve_col(x, c("antibiotic", "ab", "abx"))
  metric_col <- resolve_col(x, c("metric"))
  count_col <- resolve_col(x, c("count"))
  
  x %>%
    rename(
      sample = !!sample_col,
      antibiotic = !!ab_col,
      metric = !!metric_col,
      count = !!count_col
    ) %>%
    mutate(
      sample = normalize_sample_id(sample),
      antibiotic = as.character(antibiotic),
      metric = as.character(metric),
      count = as.numeric(count)
    )
}

# ----------------------------------------------------------------------------
# STEP 2. DERIVE PERFORMANCE TABLES
# ----------------------------------------------------------------------------

derive_prediction_tables <- function(pred_long) {
  required_metrics <- c(
    "correctS", "correctR", "falseS", "falseR",
    "zerolabelS", "zerolabelR", "twolabelS", "twolabelR",
    "S", "R", "total"
  )
  
  perf_long <- pred_long %>%
    filter(metric %in% required_metrics) %>%
    group_by(sample, antibiotic, metric) %>%
    summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      names_from = metric,
      values_from = count,
      values_fill = 0
    ) %>%
    mutate(
      correct = correctS + correctR,
      errors = falseS + falseR,
      ambiguous = zerolabelS + zerolabelR + twolabelS + twolabelR,
      called = correct + errors,
      n_true_S = S,
      n_true_R = R,
      acc_called = safe_ratio(correct, called),
      acc_total = safe_ratio(correct, total),
      error_total = safe_ratio(errors, total),
      ambiguity_total = safe_ratio(ambiguous, total),
      ME_called = safe_ratio(falseS, called),
      VME_called = safe_ratio(falseR, called),
      ME_total = safe_ratio(falseS, total),
      VME_total = safe_ratio(falseR, total),
      correctS_rate_total = safe_ratio(correctS, total),
      correctR_rate_total = safe_ratio(correctR, total)
    )
  
  perf_sample <- perf_long %>%
    group_by(sample) %>%
    summarise(
      n_antibiotics = n_distinct(antibiotic),
      total = sum(total, na.rm = TRUE),
      called = sum(called, na.rm = TRUE),
      correct = sum(correct, na.rm = TRUE),
      errors = sum(errors, na.rm = TRUE),
      ambiguous = sum(ambiguous, na.rm = TRUE),
      correctS = sum(correctS, na.rm = TRUE),
      correctR = sum(correctR, na.rm = TRUE),
      falseS = sum(falseS, na.rm = TRUE),
      falseR = sum(falseR, na.rm = TRUE),
      n_true_S = sum(n_true_S, na.rm = TRUE),
      n_true_R = sum(n_true_R, na.rm = TRUE),
      acc_called = safe_ratio(correct, called),
      acc_total = safe_ratio(correct, total),
      error_total = safe_ratio(errors, total),
      ambiguity_total = safe_ratio(ambiguous, total),
      ME_called = safe_ratio(falseS, called),
      VME_called = safe_ratio(falseR, called),
      ME_total = safe_ratio(falseS, total),
      VME_total = safe_ratio(falseR, total),
      S_fraction = safe_ratio(n_true_S, total),
      R_fraction = safe_ratio(n_true_R, total),
      phenotype_balance = abs(S_fraction - 0.5),
      .groups = "drop"
    )
  
  list(
    perf_sample_ab = perf_long,
    perf_sample = perf_sample
  )
}

# ----------------------------------------------------------------------------
# STEP 3. READ GENOTYPE GROUPING
# ----------------------------------------------------------------------------

read_genotype_groups <- function(file = GENOTYPE_GROUP_FILE) {
  if (!file.exists(file)) {
    stop("Genotype grouping file not found: ", file)
  }
  
  x <- read.csv2(
    file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  x %>%
    mutate(
      sample = normalize_sample_id(sample),
      samplegroup_functional = factor(samplegroup_functional),
      samplegroup_quinolone = factor(samplegroup_quinolone),
      samplegroup_aminoglycoside = factor(samplegroup_aminoglycoside),
      Functional_groups = as.character(Functional_groups),
      Beta_genes = as.character(Beta_genes),
      Enzyme_families = as.character(Enzyme_families),
      Quinolone_families = as.character(Quinolone_families),
      Aminoglycoside_families = as.character(Aminoglycoside_families)
    )
}

# ----------------------------------------------------------------------------
# STEP 4. MERGE DATASETS
# ----------------------------------------------------------------------------

merge_analysis_data <- function(perf, geno_groups) {
  sample_df <- perf$perf_sample %>%
    left_join(geno_groups, by = "sample")
  
  sample_ab_df <- perf$perf_sample_ab %>%
    left_join(geno_groups, by = "sample")
  
  list(
    sample_df = sample_df,
    sample_ab_df = sample_ab_df
  )
}

# ----------------------------------------------------------------------------
# STEP 5. SAMPLE-LEVEL MODELS
# ----------------------------------------------------------------------------

fit_sample_level_models <- function(sample_df) {
  model_specs <- list(
    acc_total =
      acc_total ~
      samplegroup_functional +
      samplegroup_quinolone +
      samplegroup_aminoglycoside,
    
    ambiguity_total =
      ambiguity_total ~
      samplegroup_functional +
      samplegroup_quinolone +
      samplegroup_aminoglycoside,
    
    ME_total =
      ME_total ~
      samplegroup_functional +
      samplegroup_quinolone +
      samplegroup_aminoglycoside,
    
    VME_total =
      VME_total ~
      samplegroup_functional +
      samplegroup_quinolone +
      samplegroup_aminoglycoside
  )
  
  fits <- purrr::imap(model_specs, function(form, name) {
    dat <- sample_df %>%
      filter(complete.cases(model.frame(form, data = sample_df)))
    if (nrow(dat) < 10) return(NULL)
    fit <- lm(form, data = dat)
    list(name = name, fit = fit)
  })
  
  fits <- fits[!vapply(fits, is.null, logical(1))]
  
  tidy_tbl <- purrr::map_dfr(fits, function(x) {
    broom::tidy(x$fit, conf.int = TRUE) %>%
      mutate(model = x$name, .before = 1)
  })
  
  glance_tbl <- purrr::map_dfr(fits, function(x) {
    broom::glance(x$fit) %>%
      mutate(model = x$name, .before = 1)
  })
  
  list(
    fits = fits,
    tidy = tidy_tbl,
    glance = glance_tbl
  )
}

# ----------------------------------------------------------------------------
# STEP 6. ANTIBIOTIC-SPECIFIC MODELS
# ----------------------------------------------------------------------------

fit_antibiotic_specific_models <- function(sample_ab_df, min_n = 10) {
  df <- sample_ab_df %>%
    filter(!is.na(acc_total))
  
  by_ab <- split(df, df$antibiotic)
  
  fits <- purrr::imap(by_ab, function(dat, ab) {
    dat2 <- dat %>%
      filter(complete.cases(
        acc_total,
        samplegroup_functional,
        samplegroup_quinolone,
        samplegroup_aminoglycoside
      ))
    
    if (nrow(dat2) < min_n) return(NULL)
    
    fit <- lm(
      acc_total ~
        samplegroup_functional +
        samplegroup_quinolone +
        samplegroup_aminoglycoside,
      data = dat2
    )
    
    list(antibiotic = ab, fit = fit, n = nrow(dat2))
  })
  
  fits <- fits[!vapply(fits, is.null, logical(1))]
  
  tidy_tbl <- purrr::map_dfr(fits, function(x) {
    broom::tidy(x$fit, conf.int = TRUE) %>%
      mutate(antibiotic = x$antibiotic, n = x$n, .before = 1)
  })
  
  glance_tbl <- purrr::map_dfr(fits, function(x) {
    broom::glance(x$fit) %>%
      mutate(antibiotic = x$antibiotic, n = x$n, .before = 1)
  })
  
  list(
    fits = fits,
    tidy = tidy_tbl,
    glance = glance_tbl
  )
}

# ----------------------------------------------------------------------------
# STEP 7. PLOTS
# ----------------------------------------------------------------------------

plot_sample_level <- function(sample_df) {
  p1 <- ggplot(sample_df, aes(x = samplegroup_functional, y = acc_total)) +
    geom_violin(trim = FALSE) +
    geom_jitter(width = 0.12, height = 0, alpha = 0.7) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      x = "Beta-lactamase group",
      y = "Overall accuracy",
      title = "AI prediction accuracy by beta-lactamase genotype group"
    )
  
  p2 <- ggplot(sample_df, aes(x = samplegroup_quinolone, y = acc_total)) +
    geom_violin(trim = FALSE) +
    geom_jitter(width = 0.12, height = 0, alpha = 0.7) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      x = "Quinolone genotype group",
      y = "Overall accuracy",
      title = "AI prediction accuracy by quinolone genotype group"
    )
  
  p3 <- ggplot(sample_df, aes(x = samplegroup_aminoglycoside, y = acc_total)) +
    geom_violin(trim = FALSE) +
    geom_jitter(width = 0.12, height = 0, alpha = 0.7) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      x = "Aminoglycoside genotype group",
      y = "Overall accuracy",
      title = "AI prediction accuracy by aminoglycoside genotype group"
    )
  
  p4 <- ggplot(sample_df, aes(x = samplegroup_functional, y = ambiguity_total)) +
    geom_violin(trim = FALSE) +
    geom_jitter(width = 0.12, height = 0, alpha = 0.7) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      x = "Beta-lactamase group",
      y = "Ambiguity rate",
      title = "AI ambiguity by beta-lactamase genotype group"
    )
  
  list(p1 = p1, p2 = p2, p3 = p3, p4 = p4)
}

plot_model_coefficients <- function(model_tidy, title = "Model coefficients") {
  df <- model_tidy %>%
    filter(term != "(Intercept)")
  
  ggplot(df, aes(x = estimate, y = fct_rev(term))) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_point() +
    geom_errorbar(
      aes(xmin = conf.low, xmax = conf.high),
      width = 0.2,
      orientation = "y"
    ) +
    facet_wrap(~ model, scales = "free_x") +
    theme_bw() +
    labs(
      x = "Estimate",
      y = NULL,
      title = title
    )
}

plot_ab_coefficients <- function(ab_model_tidy) {
  df <- ab_model_tidy %>%
    filter(term != "(Intercept)")
  
  ggplot(df, aes(x = antibiotic, y = estimate, ymin = conf.low, ymax = conf.high)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_pointrange() +
    facet_wrap(~ term, scales = "free_y") +
    theme_bw() +
    labs(
      x = "Antibiotic",
      y = "Coefficient estimate",
      title = "Antibiotic-specific genotype group coefficients"
    )
}

# ----------------------------------------------------------------------------
# STEP 8. HEATMAPS
# ----------------------------------------------------------------------------

plot_heatmaps <- function(sample_ab_df) {
  
  heat_beta <- sample_ab_df %>%
    group_by(antibiotic, samplegroup_functional) %>%
    summarise(
      mean_acc = mean(acc_total, na.rm = TRUE),
      n = sum(!is.na(acc_total)),
      .groups = "drop"
    )
  
  p_beta <- ggplot(heat_beta, aes(x = antibiotic, y = samplegroup_functional, fill = mean_acc)) +
    geom_tile() +
    scale_fill_viridis_c(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    theme_bw() +
    labs(
      x = "Antibiotic",
      y = "Beta-lactamase group",
      fill = "Mean accuracy",
      title = "Heatmap: accuracy by beta-lactamase group and antibiotic"
    )
  
  heat_q <- sample_ab_df %>%
    group_by(antibiotic, samplegroup_quinolone) %>%
    summarise(
      mean_acc = mean(acc_total, na.rm = TRUE),
      n = sum(!is.na(acc_total)),
      .groups = "drop"
    )
  
  p_q <- ggplot(heat_q, aes(x = antibiotic, y = samplegroup_quinolone, fill = mean_acc)) +
    geom_tile() +
    scale_fill_viridis_c(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    theme_bw() +
    labs(
      x = "Antibiotic",
      y = "Quinolone group",
      fill = "Mean accuracy",
      title = "Heatmap: accuracy by quinolone group and antibiotic"
    )
  
  heat_ag <- sample_ab_df %>%
    group_by(antibiotic, samplegroup_aminoglycoside) %>%
    summarise(
      mean_acc = mean(acc_total, na.rm = TRUE),
      n = sum(!is.na(acc_total)),
      .groups = "drop"
    )
  
  p_ag <- ggplot(heat_ag, aes(x = antibiotic, y = samplegroup_aminoglycoside, fill = mean_acc)) +
    geom_tile() +
    scale_fill_viridis_c(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    theme_bw() +
    labs(
      x = "Antibiotic",
      y = "Aminoglycoside group",
      fill = "Mean accuracy",
      title = "Heatmap: accuracy by aminoglycoside group and antibiotic"
    )
  
  heat_sample <- sample_ab_df %>%
    group_by(sample, antibiotic) %>%
    summarise(acc_total = mean(acc_total, na.rm = TRUE), .groups = "drop")
  
  p_sample <- ggplot(
    heat_sample,
    aes(
      x = antibiotic,
      y = fct_reorder(sample, acc_total, .fun = mean, na.rm = TRUE),
      fill = acc_total
    )
  ) +
    geom_tile() +
    scale_fill_viridis_c(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    theme_bw() +
    labs(
      x = "Antibiotic",
      y = "Sample",
      fill = "Accuracy",
      title = "Heatmap: sample-specific prediction accuracy"
    )
  
  list(
    heat_beta = heat_beta,
    heat_q = heat_q,
    heat_ag = heat_ag,
    heat_sample = heat_sample,
    p_beta = p_beta,
    p_q = p_q,
    p_ag = p_ag,
    p_sample = p_sample
  )
}

# ----------------------------------------------------------------------------
# STEP 9. CLUSTERING
# ----------------------------------------------------------------------------

cluster_isolates_by_prediction_pattern <- function(sample_ab_df) {
  wide <- sample_ab_df %>%
    select(sample, antibiotic, acc_total) %>%
    pivot_wider(names_from = antibiotic, values_from = acc_total)
  
  if (nrow(wide) < 3) return(NULL)
  
  mat_df <- as.data.frame(wide, stringsAsFactors = FALSE)
  rownames(mat_df) <- mat_df$sample
  mat_df$sample <- NULL
  mat <- as.matrix(mat_df)
  
  for (j in seq_len(ncol(mat))) {
    miss <- is.na(mat[, j])
    if (any(miss)) {
      col_mean <- mean(mat[, j], na.rm = TRUE)
      if (is.nan(col_mean)) col_mean <- 0
      mat[miss, j] <- col_mean
    }
  }
  
  d <- dist(mat)
  hc <- hclust(d, method = "ward.D2")
  
  cluster_tbl <- tibble(
    sample = rownames(mat),
    cluster_k3 = factor(cutree(hc, k = min(3, nrow(mat)))),
    cluster_k4 = factor(cutree(hc, k = min(4, nrow(mat))))
  )
  
  list(
    matrix = mat,
    hclust = hc,
    clusters = cluster_tbl
  )
}

plot_cluster_heatmap <- function(sample_ab_df, cluster_tbl) {
  df <- sample_ab_df %>%
    left_join(cluster_tbl, by = "sample") %>%
    group_by(sample, cluster_k3) %>%
    mutate(sample_mean_acc = mean(acc_total, na.rm = TRUE)) %>%
    ungroup()
  
  sample_order <- df %>%
    distinct(sample, cluster_k3, sample_mean_acc) %>%
    arrange(cluster_k3, sample_mean_acc) %>%
    pull(sample)
  
  df <- df %>%
    mutate(sample = factor(sample, levels = sample_order))
  
  ggplot(df, aes(x = antibiotic, y = sample, fill = acc_total)) +
    geom_tile() +
    scale_fill_viridis_c(limits = c(0, 1)) +
    facet_grid(cluster_k3 ~ ., scales = "free_y", space = "free_y") +
    theme_bw() +
    labs(
      x = "Antibiotic",
      y = "Sample",
      fill = "Accuracy",
      title = "Clustering of isolates by AI prediction pattern"
    )
}

# ----------------------------------------------------------------------------
# STEP 10. SUMMARIES
# ----------------------------------------------------------------------------

make_summary_tables <- function(sample_df, sample_ab_df) {
  overall <- sample_df %>%
    summarise(
      n_samples = n(),
      mean_acc_total = mean(acc_total, na.rm = TRUE),
      median_acc_total = median(acc_total, na.rm = TRUE),
      mean_ambiguity = mean(ambiguity_total, na.rm = TRUE),
      mean_ME_total = mean(ME_total, na.rm = TRUE),
      mean_VME_total = mean(VME_total, na.rm = TRUE)
    )
  
  by_functional <- sample_df %>%
    group_by(samplegroup_functional) %>%
    summarise(
      n_samples = n(),
      mean_acc_total = mean(acc_total, na.rm = TRUE),
      mean_ambiguity = mean(ambiguity_total, na.rm = TRUE),
      mean_ME_total = mean(ME_total, na.rm = TRUE),
      mean_VME_total = mean(VME_total, na.rm = TRUE),
      .groups = "drop"
    )
  
  by_quinolone <- sample_df %>%
    group_by(samplegroup_quinolone) %>%
    summarise(
      n_samples = n(),
      mean_acc_total = mean(acc_total, na.rm = TRUE),
      mean_ambiguity = mean(ambiguity_total, na.rm = TRUE),
      mean_ME_total = mean(ME_total, na.rm = TRUE),
      mean_VME_total = mean(VME_total, na.rm = TRUE),
      .groups = "drop"
    )
  
  by_aminoglycoside <- sample_df %>%
    group_by(samplegroup_aminoglycoside) %>%
    summarise(
      n_samples = n(),
      mean_acc_total = mean(acc_total, na.rm = TRUE),
      mean_ambiguity = mean(ambiguity_total, na.rm = TRUE),
      mean_ME_total = mean(ME_total, na.rm = TRUE),
      mean_VME_total = mean(VME_total, na.rm = TRUE),
      .groups = "drop"
    )
  
  by_antibiotic <- sample_ab_df %>%
    group_by(antibiotic) %>%
    summarise(
      n_samples = sum(!is.na(acc_total)),
      mean_acc_total = mean(acc_total, na.rm = TRUE),
      mean_ambiguity = mean(ambiguity_total, na.rm = TRUE),
      mean_ME_total = mean(ME_total, na.rm = TRUE),
      mean_VME_total = mean(VME_total, na.rm = TRUE),
      .groups = "drop"
    )
  
  list(
    overall = overall,
    by_functional = by_functional,
    by_quinolone = by_quinolone,
    by_aminoglycoside = by_aminoglycoside,
    by_antibiotic = by_antibiotic
  )
}

# ----------------------------------------------------------------------------
# STEP 11. MAIN RUNNER
# ----------------------------------------------------------------------------

run_prediction_vs_genotype_analysis <- function() {
  message2("Reading prediction counts ...")
  pred_long <- read_prediction_counts()
  
  message2("Deriving prediction performance tables ...")
  perf <- derive_prediction_tables(pred_long)
  
  message2("Reading genotype grouping ...")
  geno_groups <- read_genotype_groups(GENOTYPE_GROUP_FILE)
  
  message2("Merging prediction performance with genotype groups ...")
  merged <- merge_analysis_data(perf, geno_groups)
  
  message2("Fitting sample-level models ...")
  sample_models <- fit_sample_level_models(merged$sample_df)
  
  message2("Fitting antibiotic-specific models ...")
  ab_models <- fit_antibiotic_specific_models(merged$sample_ab_df)
  
  message2("Creating summary tables ...")
  summaries <- make_summary_tables(merged$sample_df, merged$sample_ab_df)
  
  message2("Saving tables ...")
  save_csv(pred_long, "prediction_counts_long.csv")
  save_csv(perf$perf_sample_ab, "prediction_performance_sample_antibiotic.csv")
  save_csv(perf$perf_sample, "prediction_performance_sample.csv")
  save_csv(geno_groups, "genotype_grouping_used.csv")
  save_csv(merged$sample_df, "analysis_dataset_sample.csv")
  save_csv(merged$sample_ab_df, "analysis_dataset_sample_antibiotic.csv")
  save_csv(summaries$overall, "summary_overall.csv")
  save_csv(summaries$by_functional, "summary_by_functional_group.csv")
  save_csv(summaries$by_quinolone, "summary_by_quinolone_group.csv")
  save_csv(summaries$by_aminoglycoside, "summary_by_aminoglycoside_group.csv")
  save_csv(summaries$by_antibiotic, "summary_by_antibiotic.csv")
  
  if (!is.null(sample_models$tidy)) {
    save_csv(sample_models$tidy, "sample_level_model_tidy.csv")
  }
  if (!is.null(sample_models$glance)) {
    save_csv(sample_models$glance, "sample_level_model_glance.csv")
  }
  if (!is.null(ab_models$tidy)) {
    save_csv(ab_models$tidy, "antibiotic_specific_model_tidy.csv")
  }
  if (!is.null(ab_models$glance)) {
    save_csv(ab_models$glance, "antibiotic_specific_model_glance.csv")
  }
  
  message2("Creating plots ...")
  sample_plots <- plot_sample_level(merged$sample_df)
  save_plot(sample_plots$p1, "sample_accuracy_by_functional_group.png", 8, 5)
  save_plot(sample_plots$p2, "sample_accuracy_by_quinolone_group.png", 8, 5)
  save_plot(sample_plots$p3, "sample_accuracy_by_aminoglycoside_group.png", 8, 5)
  save_plot(sample_plots$p4, "sample_ambiguity_by_functional_group.png", 8, 5)
  
  heat_obj <- plot_heatmaps(merged$sample_ab_df)
  save_csv(heat_obj$heat_beta, "heatmap_functional_group_vs_antibiotic.csv")
  save_csv(heat_obj$heat_q, "heatmap_quinolone_group_vs_antibiotic.csv")
  save_csv(heat_obj$heat_ag, "heatmap_aminoglycoside_group_vs_antibiotic.csv")
  save_csv(heat_obj$heat_sample, "heatmap_sample_vs_antibiotic.csv")
  
  save_plot(heat_obj$p_beta, "heatmap_functional_group_vs_antibiotic.png", 9, 5)
  save_plot(heat_obj$p_q, "heatmap_quinolone_group_vs_antibiotic.png", 9, 5)
  save_plot(heat_obj$p_ag, "heatmap_aminoglycoside_group_vs_antibiotic.png", 9, 5)
  save_plot(heat_obj$p_sample, "heatmap_sample_vs_antibiotic.png", 8, 9)
  
  if (!is.null(sample_models$tidy) && nrow(sample_models$tidy) > 0) {
    p_coef <- plot_model_coefficients(sample_models$tidy, "Sample-level genotype-group coefficients")
    save_plot(p_coef, "sample_level_model_coefficients.png", 10, 6)
  }
  
  if (!is.null(ab_models$tidy) && nrow(ab_models$tidy) > 0) {
    p_ab_coef <- plot_ab_coefficients(ab_models$tidy)
    save_plot(p_ab_coef, "antibiotic_specific_model_coefficients.png", 12, 8)
  }
  
  message2("Clustering isolates by prediction pattern ...")
  clust <- cluster_isolates_by_prediction_pattern(merged$sample_ab_df)
  if (!is.null(clust)) {
    save_csv(clust$clusters, "prediction_pattern_clusters.csv")
    p_cluster <- plot_cluster_heatmap(merged$sample_ab_df, clust$clusters)
    save_plot(p_cluster, "clustered_prediction_pattern_heatmap.png", 9, 9)
    
    png(file.path(OUTDIR, "plots", "prediction_pattern_dendrogram.png"),
        width = 1200, height = 800, res = 150)
    plot(clust$hclust, main = "Dendrogram: isolates clustered by prediction pattern")
    dev.off()
  }
  
  message2("Writing model summaries ...")
  if (!is.null(sample_models$fits) && length(sample_models$fits) > 0) {
    sink(file.path(OUTDIR, "models", "sample_level_models_summary.txt"))
    for (x in sample_models$fits) {
      cat("\n============================================================\n")
      cat("MODEL:", x$name, "\n")
      cat("============================================================\n")
      print(summary(x$fit))
      cat("\n")
    }
    sink()
  }
  
  if (!is.null(ab_models$fits) && length(ab_models$fits) > 0) {
    sink(file.path(OUTDIR, "models", "antibiotic_specific_models_summary.txt"))
    for (x in ab_models$fits) {
      cat("\n============================================================\n")
      cat("ANTIBIOTIC:", x$antibiotic, "\n")
      cat("N:", x$n, "\n")
      cat("============================================================\n")
      print(summary(x$fit))
      cat("\n")
    }
    sink()
  }
  
  message2("Done.")
  invisible(list(
    prediction_long = pred_long,
    perf = perf,
    geno_groups = geno_groups,
    merged = merged,
    sample_models = sample_models,
    ab_models = ab_models,
    summaries = summaries
  ))
}

# ----------------------------------------------------------------------------
# OPTIONAL QUICK INTERPRETATION TABLES
# ----------------------------------------------------------------------------

make_quick_interpretation_tables <- function(sample_df, sample_model_tidy = NULL, ab_model_tidy = NULL) {
  hardest_functional <- sample_df %>%
    group_by(samplegroup_functional) %>%
    summarise(
      n = n(),
      mean_acc = mean(acc_total, na.rm = TRUE),
      mean_ambiguity = mean(ambiguity_total, na.rm = TRUE),
      mean_VME = mean(VME_total, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(mean_acc)
  
  hardest_quinolone <- sample_df %>%
    group_by(samplegroup_quinolone) %>%
    summarise(
      n = n(),
      mean_acc = mean(acc_total, na.rm = TRUE),
      mean_ambiguity = mean(ambiguity_total, na.rm = TRUE),
      mean_VME = mean(VME_total, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(mean_acc)
  
  hardest_aminoglycoside <- sample_df %>%
    group_by(samplegroup_aminoglycoside) %>%
    summarise(
      n = n(),
      mean_acc = mean(acc_total, na.rm = TRUE),
      mean_ambiguity = mean(ambiguity_total, na.rm = TRUE),
      mean_VME = mean(VME_total, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(mean_acc)
  
  save_csv(hardest_functional, "quick_interpretation_hardest_functional_groups.csv")
  save_csv(hardest_quinolone, "quick_interpretation_hardest_quinolone_groups.csv")
  save_csv(hardest_aminoglycoside, "quick_interpretation_hardest_aminoglycoside_groups.csv")
  
  if (!is.null(sample_model_tidy)) {
    top_effects <- sample_model_tidy %>%
      filter(term != "(Intercept)") %>%
      arrange(p.value)
    
    save_csv(top_effects, "quick_interpretation_sample_model_effects.csv")
  }
  
  if (!is.null(ab_model_tidy)) {
    top_ab_effects <- ab_model_tidy %>%
      filter(term != "(Intercept)") %>%
      arrange(p.value)
    
    save_csv(top_ab_effects, "quick_interpretation_antibiotic_model_effects.csv")
  }
  
  invisible(TRUE)
}

# ----------------------------------------------------------------------------
# RUN
# ----------------------------------------------------------------------------
results <- run_prediction_vs_genotype_analysis()

make_quick_interpretation_tables(
  results$merged$sample_df,
  results$sample_models$tidy,
  results$ab_models$tidy
)