library(dplyr)
library(rlang)
source('common.R')


BASE_GROUP_COLS <- c("noinputab", "antibiotic", "significanceLevel", "sample", "mode","cpmode")
FRAC_COLS <- c("correct","ME","VME")
RATE_COLS <- c("ME","VME")

SIGNIFICANCE_LEVELS <- c("STD", "10%", "5%", "2.5%")
CONFIDENCE_LEVELS <- c("STD", "90%", "95%", "97.5%")

ANTIBIOTICS_LEVELS <- c("AMP","AMC","PIP","TZP","CAZ","CRO","CTX","FEP","CIP","OFX","LVX","MFX","GEN","TOB")
NO_INPUT_LEVELS <- c("4","5","6","7","8")
METRIC_LEVELS = c("correct", "ME", "VME")
AB_GROUPS_LEVELS = AB_GROUPS





makeWide <- function(metrics_long, which_metric) {
  library(dplyr)
  library(tidyr)
  
  which_metric <- rlang::as_string(which_metric)
  
  # First column name is the ID column (e.g. "antibiotic")
  id_name <- names(metrics_long)[1]
  
  metrics_long %>%
    dplyr::filter(.data$metric == which_metric) %>%
    tidyr::pivot_wider(
      names_from  = index,
      values_from = value
    ) %>%
    dplyr::select(-metric) %>%
    # ensure id column stays first
    dplyr::relocate(all_of(id_name))
}


groupColumns <- function(df)
{
  intersect(BASE_GROUP_COLS, names(df))
}

fracAndGroupCoulmns <- function(df)
{
  c(intersect(BASE_GROUP_COLS, names(df)),FRAC_COLS)
}

rateAndGroupCoulmns <- function(df)
{
  c(intersect(BASE_GROUP_COLS, names(df)),RATE_COLS)
}



summarise_over_and_drop <- function(df, col) {
  col_quo  <- enquo(col)
  col_name <- as_name(col_quo)
  
  # "önskade" gruppkolumner
  
  
  # de gruppkolumner som faktiskt finns i df
  GROUP_COLS <- intersect(BASE_GROUP_COLS, names(df))
  
  # om kolumnen vi ska summera över inte finns där, varna lite lätt
  if (!(col_name %in% GROUP_COLS)) {
    warning(sprintf("Kolumnen '%s' finns inte i grouping-kolumnerna. Ingen kolumn tas bort.", col_name))
  }
  
  group_cols_use <- setdiff(GROUP_COLS, col_name)
  count_cols     <- setdiff(names(df), GROUP_COLS)
  
  df %>%
    group_by(across(all_of(group_cols_use))) %>%   # NA behålls automatiskt som egen nivå
    summarise(
      across(all_of(count_cols), sum, na.rm = TRUE),
      .groups = "drop"
    )
}


countsToFracCompareWithPEK <- function(df_count=readStatisticsExcel("abErrorStatisticsCount"))
{
  cols <- fracAndGroupCoulmns(df_count)
  result <- df_count %>% mutate(correct=correct/total,ME=(falseR+notpredictedS)/total,VME=(falseS+notpredictedR)/total) %>% select(all_of(cols))
  
  # # Sanity check
  # df_Frac <- readStatisticsExcel("abErrorStatisticsFrac")
  # toCompareA <- df_Frac %>% select(all_of(FRAC_COLS))
  # toCompareB <- result %>% select(all_of(FRAC_COLS))
  # diff = toCompareA-toCompareB
  # stopifnot(all(abs(diff)<0.001))
  
  result
}


countsToErrorErrorRates <-function(df_count)
{
  cols <- rateAndGroupCoulmns(df_count)
  result <- df_count %>% mutate(ME=(falseR+zerolabelS)/S,VME=(falseS+zerolabelR)/R) %>% select(all_of(cols))
  result
}

countsToFracClinical <- function(df_count=readStatisticsExcel("abErrorStatisticsCount"))
{
  cols <- fracAndGroupCoulmns(df_count)
  result <- df_count %>% mutate(correct=correct/total,ME=falseR/total,VME=falseS/total) %>% select(all_of(cols))
  
  # # Sanity check
  # df_Frac <- readStatisticsExcel("abErrorStatisticsFrac")
  # toCompareA <- df_Frac %>% filter_and_drop(col = significanceLevel,values = NA) %>% select(all_of(FRAC_COLS))
  # toCompareB <- result %>% filter_and_drop(col = significanceLevel,values = NA) %>% select(all_of(FRAC_COLS))
  # diff = toCompareA-toCompareB
  # stopifnot(all(abs(diff)<0.001))
  
  result
}

countsToFracClinicalSumToOne <- function(df_count=readStatisticsExcel("abErrorStatisticsCount"))
{
  cols <- fracAndGroupCoulmns(df_count)
  
  result <- df_count %>% mutate(correct=correct/(falseS+falseR+correctS+correctR),ME=falseR/(falseS+falseR+correctS+correctR),VME=falseS/(falseS+falseR+correctS+correctR)) %>% select(all_of(cols))
  
  # # Sanity check
  # df_Frac <- readStatisticsExcel("abErrorStatisticsFrac")
  # toCompareA <- df_Frac %>% filter_and_drop(col = significanceLevel,values = NA) %>% select(all_of(FRAC_COLS))
  # toCompareB <- result %>% filter_and_drop(col = significanceLevel,values = NA) %>% select(all_of(FRAC_COLS))
  # diff = toCompareA-toCompareB
  # stopifnot(all(abs(diff)<0.001))
  # 
  # # Another sanity check
  # aSum <- result %>% mutate(sum = correct+ME+VME) %>% select(sum)
  # diff <- aSum - 1
  # stopifnot(all(abs(diff)<0.001))
  
  result
}

# ---- Utility: safe_div -------------------------------------------------------
# Safe division: returns NA when denominator is NA or 0
safe_div <- function(num, den) {
  ifelse(is.na(den) | den == 0, NA_real_, num / den)
}

# ---- Compute metrics from count data (TP/FP/TN/FN) ---------------------------
# Expects df to contain count columns specified in 'mapping':
#   TP = correctR, FP = falseR, TN = correctS, FN = falseS  (defaults)
#
# Returns a long table with one row per group and metric:
#   group_vars..., metric, value
compute_metrics_long_counts <- function(
    df,
    group_vars = c("noinputab", "antibiotic", "significanceLevel", "mode"),
    mapping    = c(TP = "correctR",
                   FP = "falseR",
                   TN = "correctS",
                   FN = "falseS")
) {
  missing_cols <- setdiff(unname(mapping), names(df))
  if (length(missing_cols) > 0) {
    stop(
      "These count columns from 'mapping' are missing in df: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  group_vars <- intersect(group_vars, names(df))
  
  summarised <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(unname(mapping)), ~ sum(.x, na.rm = TRUE)),
      .groups = "drop"
    )
  
  TP_col <- mapping["TP"]
  FP_col <- mapping["FP"]
  TN_col <- mapping["TN"]
  FN_col <- mapping["FN"]
  
  metrics_wide <- summarised %>%
    dplyr::mutate(
      precision = safe_div(.data[[TP_col]], .data[[TP_col]] + .data[[FP_col]]),
      recall    = safe_div(.data[[TP_col]], .data[[TP_col]] + .data[[FN_col]]),
      F1        = safe_div(
        2 * .data[[TP_col]],
        2 * .data[[TP_col]] + .data[[FP_col]] + .data[[FN_col]]
      ),
      ME        = safe_div(.data[[FP_col]], .data[[TN_col]] + .data[[FP_col]]),
      VME       = safe_div(.data[[FN_col]], .data[[TP_col]] + .data[[FN_col]]),
      correct   = safe_div(
        .data[[TP_col]] + .data[[TN_col]],
        .data[[TP_col]] + .data[[FP_col]] + .data[[TN_col]] + .data[[FN_col]]
      ),
      MCC = {
        TP <- as.double(.data[[TP_col]])
        FP <- as.double(.data[[FP_col]])
        TN <- as.double(.data[[TN_col]])
        FN <- as.double(.data[[FN_col]])
        
        denom <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
        
        safe_div(TP * TN - FP * FN, denom)
      }
    )
  
  metrics_wide %>%
    tidyr::pivot_longer(
      cols      = c(precision, recall, F1, MCC, ME, VME, correct),
      names_to  = "metric",
      values_to = "value"
    )
}

# ---- Compute unambiguous metrics (totals-aware, long format) -----------------
#
# Expects df to contain:
#   TP = correctR
#   FP = falseR
#   TN = correctS
#   FN = falseS
#   totalP = total truly resistant (positive class)
#   totalN = total truly susceptible (negative class)
#
# Metrics returned:
#   unambiguous_correct = (TP + TN) / (totalP + totalN)
#   ME  = FP / totalN
#   VME = FN / totalP
#
# Returns long table:
#   group_vars..., metric, value
#
compute_metrics_unambiguous_long_counts <- function(
    df,
    group_vars = c("noinputab", "antibiotic", "significanceLevel", "mode"),
    mapping    = c(
      TP      = "correctR",
      FP      = "falseR",
      TN      = "correctS",
      FN      = "falseS",
      total  = "total"
    ),
    check_totals = TRUE,
    tol = 0
) {
  
  required <- unname(mapping)
  missing_cols <- setdiff(required, names(df))
  if (length(missing_cols) > 0) {
    stop(
      "These required columns from 'mapping' are missing in df: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  group_vars <- intersect(group_vars, names(df))
  
  summarised <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(required), ~ sum(.x, na.rm = TRUE)),
      .groups = "drop"
    )
  
  TP_col    <- mapping[["TP"]]
  FP_col    <- mapping[["FP"]]
  TN_col    <- mapping[["TN"]]
  FN_col    <- mapping[["FN"]]
  total_col <- mapping[["total"]]

  # Optional consistency check
  if (isTRUE(check_totals)) {
    bad_total <- summarised[[total_col]] + tol <
      (summarised[[TP_col]] + summarised[[FN_col]] +
      summarised[[TN_col]] + summarised[[FP_col]])
    
    if (any(bad_total, na.rm = TRUE)) {
      stop(
        "Totals check failed: expected total >= TP+FN + TN+FP (within tol).",
        call. = FALSE
      )
    }
  }
  
  metrics_wide <- summarised %>%
    dplyr::mutate(
      correct = safe_div(
        .data[[TP_col]] + .data[[TN_col]],
        .data[[total_col]] 
      ),
      ME  = safe_div(.data[[FP_col]], .data[[total_col]]),
      VME = safe_div(.data[[FN_col]], .data[[total_col]])
    )
  
  metrics_wide %>%
    tidyr::pivot_longer(
      cols      = c(correct, ME, VME),
      names_to  = "metric",
      values_to = "value"
    )
}





convertSignificanceLevelForPresentation <- function(df) {
  df %>%
    dplyr::mutate(
      significanceLevel = dplyr::case_when(
        is.na(significanceLevel) ~ "STD",
        as.numeric(significanceLevel) == 0.1   ~ "10%",
        as.numeric(significanceLevel) == 0.05  ~ "5%",
        as.numeric(significanceLevel) == 0.025 ~ "2.5%",
        TRUE ~ as.character(significanceLevel)
      ),
      significanceLevel = factor(
        significanceLevel,
        levels = SIGNIFICANCE_LEVELS
      )
    )
}



ENSURE_READS <- function()
{
  readCountFrameWide() %>% defaults()
  readCountFrameLong() %>% defaults()
  readCountWordFrameLong() %>% defaults()
  readCountSampleFrameLong() %>% defaults()
  
  A <- readCountFrameWide() %>% makeMetricLong()
  B <- readCountFrameLong()
  stopifnot(all(A==B,na.rm = TRUE))
  
  X <- readCountFrameLong() %>% makeMetricWide()
  Y <- readCountFrameWide()
  stopifnot(all(X==Y,na.rm = TRUE))
}

ALL_GROUPS_VAR <- c("noinputab", "antibiotic", "significanceLevel", "mode","cpmode")
AB_GROUP_GROUPS_VAR <- c("noinputab", "ab_group", "significanceLevel", "mode","cpmode")
SAMPLE_GROUPS_VAR <- c("noinputab","sample", "antibiotic", "significanceLevel", "mode","cpmode")
WORD_GROUPS_VAR <- c("noinputab","word", "significanceLevel", "mode","cpmode")


PERFORMANCE_METRICS_DEFS <- tibble::tribble(
  ~metric,   ~tag, ~ylab,~y_as_percent, 
  "MCC", "A",  "Matthews correlation coefficient",FALSE,
  "F1",      "B",  "F1 score",FALSE,
  "VME",     "C",  "Very major error rate",TRUE,
  "ME",      "D",  "Major error rate",TRUE
)

VME_ME_METRICS_DEFS <- tibble::tribble(
  ~metric,   ~tag, ~ylab,~y_as_percent, 
  "VME",     "A",  "Very major error rate",TRUE,
  "ME",      "B",  "Major error rate",TRUE
)


contingencyCountsUnambiguous <- function(wideCountsFrame)
{
  wideCountsFrame %>% mutate( TP = correctR,
                              FP = falseR,
                              TN = correctS,
                              FN = falseS,
                              total = R + S) %>% select(-c(S, R,correctS, correctR, falseS, falseR, zerolabelS, zerolabelR, twolabelS, twolabelR))
}

contingencyCounts <- function(wideCountsFrame)
{
  wideCountsFrame %>% mutate( TP = correctR + twolabelR,
                              FP = falseR + zerolabelS,
                              TN = correctS + twolabelS,
                              FN = falseS + zerolabelR) %>% select(-c(S, R,correctS, correctR, falseS, falseR, zerolabelS, zerolabelR, twolabelS, twolabelR))
}



derivedMetricsFrameUnambiguous <-function(prefix=NA,countFrameWide = readCountFrameWide(prefix), vars = ALL_GROUPS_VAR)
{
  # Compute derived metrics
  compute_metrics_unambiguous_long_counts(
    contingencyCountsUnambiguous(countFrameWide),
    vars,
    mapping    = c(TP = "TP",
                   FP = "FP",
                   TN = "TN",
                   FN = "FN",
                   total = "total")
  ) %>% select(-c(TP,FP,TN,FN,total)) %>% convertSignificanceLevelForPresentation()
}

derivedMetricsFrame <-function(prefix=NA,countFrameWide = readCountFrameWide(prefix), vars = ALL_GROUPS_VAR)
{
  # Compute derived metrics
  compute_metrics_long_counts(
    contingencyCounts(countFrameWide),
    vars,
    mapping    = c(TP = "TP",
                   FP = "FP",
                   TN = "TN",
                   FN = "FN")
  ) %>% select(-c(TP,FP,TN,FN)) %>% convertSignificanceLevelForPresentation()
}

EXAMPLE_UANMBIGOUS <- function()
{
  aaa <- compute_metrics_unambiguous_long_counts(readCountFrameWideSum()) %>% filter_and_drop(mode,"Mode-A")

  #correct  
  aaa %>% filter_and_drop(significanceLevel,NA) %>% filter_and_drop(metric,c("correct"))

  #ME
  aaa %>% filter_and_drop(significanceLevel,NA) %>% filter_and_drop(metric,c("ME"))

  #VME
  aaa %>% filter_and_drop(significanceLevel,NA) %>% filter_and_drop(metric,c("VME"))

}

fetchPredictionErrors <- function()
{
  
  metrics_long_ab <- derivedMetricsFrame(countFrameWide = readCountSampleFrameWide(),vars = SAMPLE_GROUPS_VAR)
  
  errors <- metrics_long_ab %>% 
    filter_and_drop(mode,MODE_A) %>% 
    filter_and_drop(cpmode,"normal") %>%
    filter_and_drop(significanceLevel,"STD") %>%
    filter_and_drop(noinputab,6) %>%
    filter(metric %in% c("ME","VME"))
  
  error_matrix <- errors %>%
    pivot_wider(
      names_from = metric,
      values_from = value
    ) %>%
    mutate(
      n_na = rowSums(is.na(across(c(ME, VME))))
    ) %>%
    {
      if (any(.$n_na != 1)) {
        stop("Sanity check failed: each sample-antibiotic pair must have exactly one NA among ME and VME")
      }
      .
    } %>%
    mutate(
      combined = if_else(!is.na(ME), ME, -VME)
    ) %>%
    select(sample, antibiotic, combined) %>%
    mutate(
      antibiotic = factor(antibiotic, levels = ANTIBIOTICS)
    ) %>%
    complete(sample, antibiotic = factor(ANTIBIOTICS, levels = ANTIBIOTICS)) %>%
    pivot_wider(
      names_from = antibiotic,
      values_from = combined
    ) %>%
    arrange(sample) %>%
    as.data.frame()
  
  rownames(error_matrix) <- error_matrix$sample
  error_matrix$sample <- NULL
  
  error_matrix <- error_matrix[, ANTIBIOTICS]
  
  stopifnot(all(colnames(error_matrix) == ANTIBIOTICS))
  error_matrix
}



