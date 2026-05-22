# ============================================================================
# model_genotype_all_classes_mcc_me_vme.R
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
#   - antibiotic-level aggregated performance tables
#   - merged analysis datasets
#   - sample-level models
#   - antibiotic-specific models
#   - plots:
#       * MCC
#       * ME rate (= falseS / S = FNR)
#       * VME rate (= falseR / R = FPR)
#       * accuracy as secondary reference only
#   - heatmaps
#   - beta-lactamase dichotomy plots
#
# Main outcomes:
#   - MCC
#   - ME rate
#   - VME rate
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
OUTDIR <- file.path(processedRootRassembly, "model_vs_genotype")
GENOTYPE_GROUP_FILE <- file.path(indir, "sample_genotype_grouping_simple.csv")

ANTIBIOTICS <- c(
  "AMP", "AMC", "PIP", "TZP",
  "CAZ", "CRO", "CTX", "FEP",
  "CIP", "OFX", "LVX", "MFX",
  "GEN", "TOB"
)

REQUIRED_PREDICTION_METRICS <- c(
  "correctS", "correctR", "falseS", "falseR",
  "zerolabelS", "zerolabelR", "twolabelS", "twolabelR",
  "S", "R", "total"
)

BETA_DICHOTOMIES <- c(
  "has_betalactamase",
  "has_ESBL",
  "has_AmpC_group",
  "has_OXA_group",
  "has_TEM_SHV_like"
)

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

safe_div <- function(num, den) {
  ifelse(is.na(den) | den == 0, NA_real_, num / den)
}


# MCC defined for all confusion matrices
safe_mcc <- function(tp, tn, fp, fn) {
  tp <- as.numeric(tp)
  tn <- as.numeric(tn)
  fp <- as.numeric(fp)
  fn <- as.numeric(fn)
  
  n <- length(tp)
  out <- rep(NA_real_, n)
  
  for (i in seq_len(n)) {
    TPi <- tp[i]
    TNi <- tn[i]
    FPi <- fp[i]
    FNi <- fn[i]
    
    vals <- c(TPi, TNi, FPi, FNi)
    
    if (any(is.na(vals))) {
      out[i] <- NA_real_
      next
    }
    
    num <- TPi * TNi - FPi * FNi
    den <- sqrt((TPi + FPi) * (TPi + FNi) * (TNi + FPi) * (TNi + FNi))
    
    if (!is.na(den) && den > 0) {
      out[i] <- num / den
      next
    }
    
    nonzero_n <- sum(vals > 0)
    
    if (nonzero_n == 1) {
      if (TPi > 0 || TNi > 0) {
        out[i] <- 1
      } else if (FPi > 0 || FNi > 0) {
        out[i] <- -1
      } else {
        out[i] <- 0
      }
      next
    }
    
    out[i] <- 0
  }
  
  out
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
    dpi = 300
  )
}

empty_to_label <- function(x, label) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- label
  x
}

# ----------------------------------------------------------------------------
# STEP 1. READ AI PREDICTION COUNTS ONCE
# ----------------------------------------------------------------------------

read_prediction_counts <- function() {
  if (!exists("readCountSampleFrameLong", mode = "function")) {
    stop("Function readCountSampleFrameLong() is not available in the environment.")
  }
  
  x <- readCountSampleFrameLong() %>%
    filter_and_drop(cpmode, "normal") %>%
    filter_and_drop(mode, MODE_A) %>%
    filter_and_drop(significanceLevel, NA) %>%
    filter_and_drop(noinputab,6)
  
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
      antibiotic = factor(as.character(antibiotic), levels = ANTIBIOTICS),
      metric = as.character(metric),
      count = as.numeric(count)
    ) %>%
    filter(!is.na(antibiotic))
}

# ----------------------------------------------------------------------------
# STEP 2. DERIVE PERFORMANCE TABLES
# ----------------------------------------------------------------------------

derive_prediction_tables <- function(pred_long) {
  pred_use <- pred_long %>%
    filter(metric %in% REQUIRED_PREDICTION_METRICS)
  
  perf_sample_ab <- pred_use %>%
    group_by(sample, antibiotic, metric) %>%
    summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      names_from = metric,
      values_from = count,
      values_fill = 0
    ) %>%
    mutate(
      TP = correctR,
      TN = correctS,
      FP = falseR,
      FN = falseS,
      
      correct = TP + TN,
      errors = FP + FN,
      ambiguous = zerolabelS + zerolabelR + twolabelS + twolabelR,
      called = correct + errors,
      
      accuracy = safe_div(TP + TN, TP + TN + FP + FN),
      
      ME  = safe_div(FP, TN + FP),
      VME = safe_div(FN, TP + FN),
      
      mcc = safe_mcc(TP, TN, FP, FN)
    ) %>%
    mutate(
      antibiotic = factor(as.character(antibiotic), levels = ANTIBIOTICS)
    )
  
  perf_sample <- pred_use %>%
    group_by(sample, metric) %>%
    summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      names_from = metric,
      values_from = count,
      values_fill = 0
    ) %>% mutate(
      TP = correctR,
      TN = correctS,
      FP = falseR,
      FN = falseS,
      
      correct = TP + TN,
      errors = FP + FN,
      ambiguous = zerolabelS + zerolabelR + twolabelS + twolabelR,
      called = correct + errors,
      
      accuracy = safe_div(TP + TN, TP + TN + FP + FN),
      
      ME  = safe_div(FP, TN + FP),
      VME = safe_div(FN, TP + FN),
      
      mcc = safe_mcc(TP, TN, FP, FN)
    )
  
  n_ab_tbl <- pred_use %>%
    distinct(sample, antibiotic) %>%
    count(sample, name = "n_antibiotics")
  
  perf_sample <- perf_sample %>%
    left_join(n_ab_tbl, by = "sample")
  
  perf_ab <- pred_use %>%
    group_by(antibiotic, metric) %>%
    summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      names_from = metric,
      values_from = count,
      values_fill = 0
    ) %>%
    mutate(
      TP = correctR,
      TN = correctS,
      FP = falseR,
      FN = falseS,
      
      correct = TP + TN,
      errors = FP + FN,
      ambiguous = zerolabelS + zerolabelR + twolabelS + twolabelR,
      called = correct + errors,
      
      accuracy = safe_div(TP + TN, TP + TN + FP + FN),
      
      ME  = safe_div(FP, TN + FP),
      VME = safe_div(FN, TP + FN),
      
      mcc = safe_mcc(TP, TN, FP, FN)
    ) %>%
    mutate(
      antibiotic = factor(as.character(antibiotic), levels = ANTIBIOTICS)
    ) %>%
    arrange(antibiotic)
  
  list(
    perf_sample_ab = perf_sample_ab,
    perf_sample = perf_sample,
    perf_ab = perf_ab
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
      samplegroup_functional = factor(empty_to_label(samplegroup_functional, "FG00")),
      samplegroup_quinolone = factor(empty_to_label(samplegroup_quinolone, "QG00")),
      samplegroup_aminoglycoside = factor(empty_to_label(samplegroup_aminoglycoside, "AG00")),
      Functional_groups = empty_to_label(Functional_groups, ""),
      Beta_genes = empty_to_label(Beta_genes, ""),
      Enzyme_families = empty_to_label(Enzyme_families, ""),
      Quinolone_families = empty_to_label(Quinolone_families, ""),
      Aminoglycoside_families = empty_to_label(Aminoglycoside_families, "")
    ) %>%
    mutate(
      has_betalactamase = Functional_groups != "",
      has_ESBL = str_detect(Functional_groups, fixed("2be")),
      has_OXA_group = str_detect(Functional_groups, fixed("2d")),
      has_AmpC_group = str_detect(Functional_groups, fixed("1")),
      has_TEM_SHV_like = str_detect(Functional_groups, fixed("2b")),
      
      has_QRDR = str_detect(Quinolone_families, fixed("QRDR")),
      has_qnr = str_detect(Quinolone_families, regex("Qnr", ignore_case = TRUE)),
      has_oqx = str_detect(Quinolone_families, fixed("Oqx")),
      has_qepA = str_detect(Quinolone_families, fixed("QepA")),
      
      has_AAC = str_detect(Aminoglycoside_families, fixed("AAC")),
      has_APH = str_detect(Aminoglycoside_families, fixed("APH")),
      has_AAD = str_detect(Aminoglycoside_families, fixed("AAD")),
      has_ANT = str_detect(Aminoglycoside_families, fixed("ANT")),
      has_16S_rmt = str_detect(Aminoglycoside_families, regex("ArmA|Rmt", ignore_case = TRUE))
    )
}

# ----------------------------------------------------------------------------
# STEP 4. MERGE DATASETS
# ----------------------------------------------------------------------------

merge_analysis_data <- function(perf, geno_groups) {
  sample_df <- perf$perf_sample %>%
    left_join(geno_groups, by = "sample")
  
  sample_ab_df <- perf$perf_sample_ab %>%
    left_join(geno_groups, by = "sample") %>%
    mutate(
      antibiotic = factor(as.character(antibiotic), levels = ANTIBIOTICS)
    )
  
  list(
    sample_df = sample_df,
    sample_ab_df = sample_ab_df
  )
}


# ----------------------------------------------------------------------------
# STEP 6. PLOTS
# ----------------------------------------------------------------------------


build_group_antibiotic_metrics <- function(sample_ab_df, group_col) {
  group_sym <- rlang::sym(group_col)
  
  sample_ab_df %>%
    group_by(antibiotic, !!group_sym) %>%
    summarise(
      TP = sum(TP, na.rm = TRUE),
      TN = sum(TN, na.rm = TRUE),
      FP = sum(FP, na.rm = TRUE),
      FN = sum(FN, na.rm = TRUE),
      correctS = sum(correctS, na.rm = TRUE),
      correctR = sum(correctR, na.rm = TRUE),
      falseS = sum(falseS, na.rm = TRUE),
      falseR = sum(falseR, na.rm = TRUE),
      S = sum(S, na.rm = TRUE),
      R = sum(R, na.rm = TRUE),
      total = sum(total, na.rm = TRUE),
      ambiguous = sum(ambiguous, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      mcc = safe_mcc(TP, TN, FP, FN),
      ME  = safe_div(FP, TN + FP),
      VME = safe_div(FN, TP + FN),
      accuracy = safe_div(TP + TN, TP + TN + FP + FN),
      antibiotic = factor(as.character(antibiotic), levels = ANTIBIOTICS)
    )
}

plot_heatmaps <- function(sample_ab_df) {
  heat_beta <- build_group_antibiotic_metrics(sample_ab_df, "samplegroup_functional") %>%
    mutate(samplegroup_functional = factor(as.character(samplegroup_functional),
                                           levels = sort(unique(as.character(samplegroup_functional)))))
  
  p_beta_mcc <- ggplot(heat_beta, aes(x = antibiotic, y = samplegroup_functional, fill = mcc)) +
    geom_tile() + scale_fill_viridis_c(limits = c(-1, 1)) + theme_bw() +
    labs(x = "Antibiotic", y = "Beta-lactamase group", fill = "MCC",
         title = "Heatmap: MCC by beta-lactamase group and antibiotic")
  
  p_beta_me <- ggplot(heat_beta, aes(x = antibiotic, y = samplegroup_functional, fill = ME)) +
    geom_tile() + scale_fill_viridis_c(limits = c(0, 1), na.value = "grey90") + theme_bw() +
    labs(x = "Antibiotic", y = "Beta-lactamase group", fill = "ME rate",
         title = "Heatmap: ME rate by beta-lactamase group and antibiotic")
  
  p_beta_vme <- ggplot(heat_beta, aes(x = antibiotic, y = samplegroup_functional, fill = VME)) +
    geom_tile() + scale_fill_viridis_c(limits = c(0, 1), na.value = "grey90") + theme_bw() +
    labs(x = "Antibiotic", y = "Beta-lactamase group", fill = "VME rate",
         title = "Heatmap: VME rate by beta-lactamase group and antibiotic")
  
  list(
    heat_beta = heat_beta,
    p_beta_mcc = p_beta_mcc,
    p_beta_me = p_beta_me,
    p_beta_vme = p_beta_vme
  )
}

plot_overall_reference <- function(sample_df, sample_ab_df, perf_ab) {
  p1 <- ggplot(sample_df, aes(x = mcc)) +
    geom_histogram(binwidth = 0.1, boundary = 0, closed = "left") +
    theme_bw() +
    labs(x = "Sample-level MCC", y = "Number of samples",
         title = "Distribution of sample-level MCC across all samples")
  
  ranked_df <- sample_df %>%
    arrange(mcc, sample) %>%
    mutate(sample_rank = row_number())
  
  p2 <- ggplot(ranked_df, aes(x = sample_rank, y = mcc)) +
    geom_point() +
    theme_bw() +
    labs(x = "Samples ranked by MCC", y = "Sample-level MCC",
         title = "Sample-level MCC across all samples")
  
  p3 <- ggplot(perf_ab, aes(x = antibiotic, y = accuracy)) +
    geom_col() +
    theme_bw() +
    labs(x = "Antibiotic", y = "Aggregated accuracy",
         title = "Aggregated accuracy per antibiotic")
  
  p4 <- ggplot(perf_ab, aes(x = antibiotic, y = mcc)) +
    geom_col() +
    theme_bw() +
    labs(x = "Antibiotic", y = "Aggregated MCC",
         title = "Aggregated MCC per antibiotic")
  
  p_me <- ggplot(perf_ab, aes(x = antibiotic, y = ME)) +
    geom_col() + theme_bw() +
    labs(x = "Antibiotic", y = "Aggregated ME",
         title = "Aggregated ME per antibiotic")
  
  p_vme <- ggplot(perf_ab, aes(x = antibiotic, y = VME)) +
    geom_col() + theme_bw() +
    labs(x = "Antibiotic", y = "Aggregated VME",
         title = "Aggregated VME per antibiotic")
  
  list(
    ranked_df = ranked_df,
    perf_ab = perf_ab,
    p1 = p1,
    p2 = p2,
    p3 = p3,
    p4 = p4,
    p_me = p_me,
    p_vme = p_vme
  )
}


# ----------------------------------------------------------------------------
# STEP 8. SUMMARIES
# ----------------------------------------------------------------------------

make_summary_tables <- function(sample_df, perf_ab) {
  overall <- sample_df %>%
    summarise(
      n_samples = n(),
      mean_sample_mcc = mean(mcc, na.rm = TRUE),
      median_sample_mcc = median(mcc, na.rm = TRUE),
      mean_sample_me = mean(ME, na.rm = TRUE),
      mean_sample_vme = mean(VME, na.rm = TRUE),
      mean_accuracy = mean(accuracy, na.rm = TRUE)
    )
  
  by_functional <- sample_df %>%
    group_by(samplegroup_functional) %>%
    summarise(
      n_samples = n(),
      mean_sample_mcc = mean(mcc, na.rm = TRUE),
      mean_sample_me = mean(ME, na.rm = TRUE),
      mean_sample_vme = mean(VME, na.rm = TRUE),
      mean_accuracy = mean(accuracy, na.rm = TRUE),
      .groups = "drop"
    )
  
  by_antibiotic <- perf_ab %>%
    select(antibiotic, TP, TN, FP, FN, S, R, total, mcc, ME, VME, accuracy)
  
  list(
    overall = overall,
    by_functional = by_functional,
    by_antibiotic = by_antibiotic
  )
}

# ----------------------------------------------------------------------------
# STEP 9. BETA-LACTAMASE DICHOTOMY PLOTS
# Correct aggregation: pred_long -> split on sample -> aggregate -> metric
# ----------------------------------------------------------------------------

build_antibiotic_metrics_by_dichotomy_from_predlong <- function(pred_long, geno_groups, dichotomy_var) {
  dich_sym <- rlang::sym(dichotomy_var)
  
  pred_long %>%
    filter(metric %in% REQUIRED_PREDICTION_METRICS) %>%
    left_join(
      geno_groups %>% select(sample, !!dich_sym),
      by = "sample"
    ) %>%
    mutate(
      dichotomy = dplyr::if_else(!!dich_sym, "Yes", "No"),
      antibiotic = factor(as.character(antibiotic), levels = ANTIBIOTICS)
    ) %>%
    group_by(antibiotic, dichotomy, metric) %>%
    summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      names_from = metric,
      values_from = count,
      values_fill = 0
    ) %>%
    mutate(
      TP = correctR,
      TN = correctS,
      FP = falseR,
      FN = falseS,
      
      correct = TP + TN,
      errors = FP + FN,
      ambiguous = zerolabelS + zerolabelR + twolabelS + twolabelR,
      called = correct + errors,
      
      accuracy = safe_div(TP + TN, TP + TN + FP + FN),
      ME  = safe_div(FP, TN + FP),
      VME = safe_div(FN, TP + FN),
      mcc = safe_mcc(TP, TN, FP, FN),
      
      dichotomy_var = dichotomy_var,
      dichotomy = factor(dichotomy, levels = c("No", "Yes"))
    ) %>%
    arrange(antibiotic, dichotomy)
}

plot_antibiotic_metric_by_dichotomy <- function(pred_long,
                                                geno_groups,
                                                dichotomy_var,
                                                metric = c("mcc", "ME", "VME"),
                                                pretty_name = dichotomy_var) {
  metric <- match.arg(metric)
  
  plot_df <- build_antibiotic_metrics_by_dichotomy_from_predlong(pred_long, geno_groups, dichotomy_var)
  count_df <- build_sample_counts_by_dichotomy(geno_groups, dichotomy_var)
  
  n_yes <- count_df %>% filter(dichotomy == "Yes") %>% pull(n)
  n_no  <- count_df %>% filter(dichotomy == "No") %>% pull(n)
  
  y_lab <- dplyr::case_when(
    metric == "mcc" ~ "Aggregated MCC",
    metric == "ME"  ~ "ME = falseR / (correctS + falseR)",
    metric == "VME" ~ "VME = falseS / (correctR + falseS)"
  )
  
  title_prefix <- dplyr::case_when(
    metric == "mcc" ~ "Aggregated MCC by antibiotic and",
    metric == "ME"  ~ "ME by antibiotic and",
    metric == "VME" ~ "VME by antibiotic and"
  )
  
  p <- ggplot(plot_df, aes(x = antibiotic, y = .data[[metric]], fill = dichotomy)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.72) +
    theme_bw() +
    labs(
      x = "Antibiotic",
      y = y_lab,
      fill = pretty_name,
      title = sprintf("%s %s (yes=%d, no=%d)", title_prefix, pretty_name, n_yes, n_no)
    )
  
  list(
    data = plot_df,
    counts = count_df,
    plot = p
  )
}

build_sample_counts_by_dichotomy <- function(geno_groups, dichotomy_var) {
  dich_sym <- rlang::sym(dichotomy_var)
  
  geno_groups %>%
    distinct(sample, !!dich_sym) %>%
    mutate(dichotomy = dplyr::if_else(!!dich_sym, "Yes", "No")) %>%
    count(dichotomy, name = "n") %>%
    tidyr::complete(dichotomy = c("Yes", "No"), fill = list(n = 0)) %>%
    mutate(dichotomy = factor(dichotomy, levels = c("No", "Yes")))
}


plot_all_beta_dichotomy_metric_bars <- function(pred_long,
                                                geno_groups,
                                                metric = c("mcc", "ME", "VME")) {
  metric <- match.arg(metric)
  
  dichotomy_labels <- c(
    has_betalactamase = "Any beta-lactamase",
    has_ESBL = "ESBL",
    has_AmpC_group = "AmpC",
    has_OXA_group = "OXA",
    has_TEM_SHV_like = "TEM/SHV-like"
  )
  
  plot_tbl <- purrr::map_dfr(
    names(dichotomy_labels),
    function(v) {
      tmp_plot <- build_antibiotic_metrics_by_dichotomy_from_predlong(pred_long, geno_groups, v)
      tmp_counts <- build_sample_counts_by_dichotomy(geno_groups, v)
      
      n_yes <- tmp_counts %>% filter(dichotomy == "Yes") %>% pull(n)
      n_no  <- tmp_counts %>% filter(dichotomy == "No") %>% pull(n)
      
      tmp_plot %>%
        mutate(
          dichotomy_label = dichotomy_labels[[v]],
          title_label = sprintf("%s (yes=%d, no=%d)", dichotomy_labels[[v]], n_yes, n_no)
        )
    }
  ) %>%
    mutate(
      title_label = factor(title_label, levels = unique(title_label)),
      antibiotic = factor(as.character(antibiotic), levels = ANTIBIOTICS),
      dichotomy = factor(dichotomy, levels = c("No", "Yes"))
    )
  
  y_lab <- dplyr::case_when(
    metric == "mcc" ~ "MCC",
    metric == "ME"  ~ "ME",
    metric == "VME" ~ "VME"
  )
  
  title_txt <- dplyr::case_when(
    metric == "mcc" ~ "MCC by antibiotic and beta-lactamase dichotomy",
    metric == "ME"  ~ "ME by antibiotic and beta-lactamase dichotomy",
    metric == "VME" ~ "VME by antibiotic and beta-lactamase dichotomy"
  )
  subtitle_txt <- dplyr::case_when(
    metric == "mcc" ~ NA_character_,
    metric == "ME"  ~ "Bars are omitted where no true susceptible isolates were present",
    metric == "VME" ~ "Bars are omitted where no true resistant isolates were present"
  )
  
  p <- ggplot(plot_tbl, aes(x = antibiotic, y = .data[[metric]], fill = dichotomy)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.72) +
    facet_wrap(~ title_label, ncol = 1) +
    theme_bw() +
    labs(
      x = "Antibiotic",
      y = y_lab,
      fill = NULL,
      title = title_txt,
      subtitle = subtitle_txt
    )
  
  list(
    data = plot_tbl,
    plot = p
  )
}


# ----------------------------------------------------------------------------
# STEP 10. MAIN RUNNER
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

  message2("Creating summary tables ...")
  summaries <- make_summary_tables(merged$sample_df, perf$perf_ab)
  
  message2("Saving tables ...")
  save_csv(pred_long, "prediction_counts_long.csv")
  save_csv(perf$perf_sample_ab, "prediction_performance_sample_antibiotic.csv")
  save_csv(perf$perf_sample, "prediction_performance_sample.csv")
  save_csv(perf$perf_ab, "prediction_performance_antibiotic_aggregated.csv")
  save_csv(geno_groups, "genotype_grouping_used.csv")
  save_csv(merged$sample_df, "analysis_dataset_sample.csv")
  save_csv(merged$sample_ab_df, "analysis_dataset_sample_antibiotic.csv")
  save_csv(summaries$overall, "summary_overall.csv")
  save_csv(summaries$by_functional, "summary_by_functional_group.csv")
  save_csv(summaries$by_antibiotic, "summary_by_antibiotic.csv")
  
  message2("Creating sample-level plots ...")

  message2("Creating heatmaps ...")
  heat_obj <- plot_heatmaps(merged$sample_ab_df)
  save_csv(heat_obj$heat_beta, "heatmap_functional_group_vs_antibiotic.csv")
  save_plot(heat_obj$p_beta_mcc, "heatmap_functional_group_vs_antibiotic_mcc.png", 9, 5)
  save_plot(heat_obj$p_beta_me, "heatmap_functional_group_vs_antibiotic_me.png", 9, 5)
  save_plot(heat_obj$p_beta_vme, "heatmap_functional_group_vs_antibiotic_vme.png", 9, 5)
  
  message2("Creating overall reference plots ...")
  overall_ref <- plot_overall_reference(merged$sample_df, merged$sample_ab_df, perf$perf_ab)
  #save_csv(overall_ref$ranked_df, "overall_sample_mcc_ranked.csv")
  save_csv(overall_ref$perf_ab, "overall_antibiotic_aggregated_metrics.csv")
  #save_plot(overall_ref$p1, "overall_sample_mcc_histogram.png", 7, 5)
  #save_plot(overall_ref$p2, "overall_sample_mcc_ranked.png", 8, 5)
  save_plot(overall_ref$p3, "overall_antibiotic_aggregated_accuracy.png", 10, 5)
  save_plot(overall_ref$p4, "overall_antibiotic_aggregated_mcc.png", 10, 5)
  save_plot(overall_ref$p_me,  "overall_antibiotic_aggregated_ME.png", 10, 5)
  save_plot(overall_ref$p_vme, "overall_antibiotic_aggregated_VME.png", 10, 5)
  
  message2("Creating beta-lactamase dichotomy plots ...")
  
  # MCC
  beta_dichotomy_mcc_all <- plot_all_beta_dichotomy_metric_bars(pred_long, geno_groups, metric = "mcc")
  save_csv(beta_dichotomy_mcc_all$data, "beta_dichotomy_antibiotic_aggregated_mcc.csv")
  save_plot(beta_dichotomy_mcc_all$plot, "beta_dichotomy_antibiotic_aggregated_mcc.png", 11, 12)
  
  # ME
  beta_dichotomy_me_all <- plot_all_beta_dichotomy_metric_bars(pred_long, geno_groups, metric = "ME")
  save_csv(beta_dichotomy_me_all$data, "beta_dichotomy_antibiotic_me_rate.csv")
  save_plot(beta_dichotomy_me_all$plot, "beta_dichotomy_antibiotic_me_rate.png", 11, 12)
  
  # VME
  beta_dichotomy_vme_all <- plot_all_beta_dichotomy_metric_bars(pred_long, geno_groups, metric = "VME")
  save_csv(beta_dichotomy_vme_all$data, "beta_dichotomy_antibiotic_vme_rate.csv")
  save_plot(beta_dichotomy_vme_all$plot, "beta_dichotomy_antibiotic_vme_rate.png", 11, 12)
  
  # # Individual plots
  # dichotomy_labels <- c(
  #   has_betalactamase = "Any beta-lactamase",
  #   has_ESBL = "ESBL",
  #   has_AmpC_group = "AmpC",
  #   has_OXA_group = "OXA",
  #   has_TEM_SHV_like = "TEM/SHV-like"
  # )
  # 
  # for (v in names(dichotomy_labels)) {
  #   pretty <- dichotomy_labels[[v]]
  #   
  #   p_mcc <- plot_antibiotic_metric_by_dichotomy(pred_long, geno_groups, v, metric = "mcc", pretty_name = pretty)
  #   save_csv(p_mcc$data, paste0("beta_dichotomy_", v, "_mcc.csv"))
  #   save_csv(p_mcc$counts, paste0("beta_dichotomy_", v, "_counts.csv"))
  #   save_plot(p_mcc$plot, paste0("beta_dichotomy_", v, "_mcc.png"), 10, 5)
  #   
  #   p_me <- plot_antibiotic_metric_by_dichotomy(pred_long, geno_groups, v, metric = "ME", pretty_name = pretty)
  #   save_csv(p_me$data, paste0("beta_dichotomy_", v, "_me_rate.csv"))
  #   save_plot(p_me$plot, paste0("beta_dichotomy_", v, "_me_rate.png"), 10, 5)
  #   
  #   p_vme <- plot_antibiotic_metric_by_dichotomy(pred_long, geno_groups, v, metric = "VME", pretty_name = pretty)
  #   save_csv(p_vme$data, paste0("beta_dichotomy_", v, "_vme_rate.csv"))
  #   save_plot(p_vme$plot, paste0("beta_dichotomy_", v, "_vme_rate.png"), 10, 5)
  # }
  
  
  message2("Done.")
  invisible(list(
    prediction_long = pred_long,
    perf = perf,
    geno_groups = geno_groups,
    merged = merged,
    summaries = summaries
  ))
}

# ----------------------------------------------------------------------------
# OPTIONAL QUICK INTERPRETATION TABLES
# ----------------------------------------------------------------------------

make_quick_interpretation_tables <- function(sample_df) {
  hardest_functional <- sample_df %>%
    group_by(samplegroup_functional) %>%
    summarise(
      n = n(),
      mean_sample_mcc = mean(mcc, na.rm = TRUE),
      mean_sample_me = mean(ME, na.rm = TRUE),
      mean_sample_vme = mean(VME, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(mean_sample_mcc)
  
  save_csv(hardest_functional, "quick_interpretation_hardest_functional_groups.csv")
  
  invisible(TRUE)
}

# ----------------------------------------------------------------------------
# RUN
# ----------------------------------------------------------------------------

RUN <- function() {
  results <- run_prediction_vs_genotype_analysis()
  
  make_quick_interpretation_tables(
    results$merged$sample_df)
  
  invisible(results)
}