source(file = 'manuscript/manuscriptcommon.R')
source(file = 'manuscript/selectExcelStatisticsTables.R')
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(forcats)
library(readxl)
source(file = 'manuscript/plotOverallABChatGPT.R')
source(file = 'manuscript/plotAntibioticABChatGPT.R')
source(file = 'manuscript/metric_plot_helpers.R')
source(file = 'manuscript/metric_helpers.R')

selectedMode <- MODE_A


computeIndexMetricsLong <- function(
    mode,
    fetch_fun = fetchfractionMetricAntibioticsVsIndex,
    mapping = c(TP = "correctR",
                FP = "falseR",
                TN = "correctS",
                FN = "falseS")
) {
  library(dplyr)
  library(tidyr)
  library(purrr)
  
  safe_div <- function(num, den) ifelse(den == 0, NA_real_, num / den)
  
  # 1) Fetch frames (each: first column = id, rest = indices)
  raw_frames <- imap(mapping, ~ fetch_fun(mode, .x))
  
  # Name of the first (ID) column, e.g. "antibiotic"
  id_name <- names(raw_frames[[1]])[1]
  
  # 2) Stack them long with the original id column name preserved
  counts_long <- imap_dfr(
    raw_frames,
    function(df, metric_name) {
      df %>%
        tidyr::pivot_longer(
          cols = -all_of(id_name),
          names_to  = "index",
          values_to = "count"
        ) %>%
        mutate(metric_component = metric_name)
    }
  )
  
  # 3) Go wide on TP/FP/TN/FN, compute metrics
  metrics_wide <- counts_long %>%
    tidyr::pivot_wider(
      names_from  = metric_component,
      values_from = count,
      values_fill = 0
    ) %>%
    mutate(
      precision = safe_div(TP, TP + FP),
      recall    = safe_div(TP, TP + FN),
      F1        = safe_div(2 * TP, 2 * TP + FP + FN),
      ME        = safe_div(FP, TN + FP),
      VME       = safe_div(FN, TP + FN),
      correct   = safe_div(TP + TN, TP + FP + TN + FN),
    )
  
  # 4) Final tidy-long: keep original id_name
  metrics_long <- metrics_wide %>%
    tidyr::pivot_longer(
      cols      = c(precision, recall, F1, ME, VME,correct),
      names_to  = "metric",
      values_to = "value"
    ) %>%
    dplyr::relocate(all_of(id_name), index, metric, value)
  
  metrics_long %>% select(- c(TP,    FP,    TN,    FN))
}




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



fetchCorrectSignificanceLevelVsIndex <- function(mode  = MODE_A)
{
  frame <- readStatisticsExcel("errorStatisticsFrac")
  subFrame <- frame
  subFrame <- subFrame %>% 
    filterMode(mode) %>% 
    filterIndexRange(c(4:8)) %>% 
    filterMetricRange("correct") %>% 
    filterSignificanceRange(SELECTED_SIGNIS)
  subFrame <- subFrame %>% tidyr::pivot_wider(names_from = noinputab ,values_from = correct)
  
  
  
  subFrame <- subFrame %>% dplyr::slice(4:1)
  subFrame$significanceLevel <- unlist(lapply(subFrame$significanceLevel,function(x) {signiToPercent(x)}))
  subFrame
}


fetchMetricVsSignificanceLevel <- function(mode  = MODE_A)
{
  indexRange <- c(4,5,6,7,8)
  frame <- readStatisticsExcel("errorStatisticsFrac")
  subFrame <- frame
  subFrame <- subFrame %>% 
    filterMode(mode) %>% 
    filterIndexRange(indexRange) %>% 
    filterSignificanceRange(SELECTED_SIGNIS)
  subFrame <- subFrame %>% tidyr::pivot_longer(cols = all_of(METRICS_COLS),names_to = "metric" )
  subFrame$significanceLevel <- unlist(lapply(subFrame$significanceLevel,function(x) {signiToPercent(x)}))
  subFrame <- subFrame %>% tidyr::pivot_wider(names_from = significanceLevel ,values_from = value)
  subFrame <- subFrame %>% dplyr::select(noinputab, metric ,`non-conformal`,`10%`,`5%`,`2.5%`)  
  aList <- lapply(indexRange, function(x) { subFrame %>% filterIndexRange(x)})
  names(aList) <-  as.character(indexRange)

  aList  
}


fetchCorrectAntibioticsVsIndex <- function(mode  = "Mode-A")
{
  frame <- readStatisticsExcel("abErrorStatisticsFrac") 
  subFrame <- frame
  
  indexRange <- c(4:8)
  subFrame <- subFrame %>% filterMode(mode) %>% filterSignificanceRange(NA) %>% filterIndexRange(indexRange) %>% filterMetricRange("correct")
  subFrame <- subFrame %>% tidyr::pivot_wider(names_from = noinputab,values_from = correct)
  
  subFrame
}


fetchfractionMetricAntibioticsGroupVsIndex <- function(mode  = "Mode-A",fractionMetric)
{
  frame <- countByAntibioticsGroup()
  subFrame <- frame
  #FIXME downstream cant handle ab_group
  subFrame <- subFrame %>% rename(antibiotic = ab_group)
  indexRange <- c(4:8)
  subFrame <- subFrame %>% filterMode(mode) %>% filterSignificanceRange(NA) %>% filterIndexRange(indexRange) %>% filterMetricRange(fractionMetric)
  subFrame <- subFrame %>% tidyr::pivot_wider(names_from = noinputab,values_from = dplyr::all_of(fractionMetric))
  subFrame
}


fetchfractionMetricAntibioticsVsIndex <- function(mode  = "Mode-A",fractionMetric)
{
  frame <- readStatisticsExcel("abErrorStatisticsCount") 
  subFrame <- frame
  indexRange <- c(4:8)
  subFrame <- subFrame %>% filterMode(mode) %>% filterSignificanceRange(NA) %>% filterIndexRange(indexRange) %>% filterMetricRange(fractionMetric)
  subFrame <- subFrame %>% tidyr::pivot_wider(names_from = noinputab,values_from = dplyr::all_of(fractionMetric))
  subFrame
}


fetchfractionMetricAntibioticsVsMode <- function(index  = 6,fractionMetric = "correctS")
{
  frame <- readStatisticsExcel("abErrorStatisticsCount") 
  subFrame <- frame
  subFrame <- subFrame %>% filterSignificanceRange(NA) %>% filterIndexRange(index) %>% filterMetricRange(fractionMetric)
  subFrame <- subFrame %>% tidyr::pivot_wider(names_from = mode,values_from = dplyr::all_of(fractionMetric))
  subFrame
}


fetchfractionMetricAntibioticsVsSignicanceLevel <- function(mode  = "Mode-A",fractionMetric = "correctR")
{
  frame <- readStatisticsExcel("abErrorStatisticsCount") 
  subFrame <- frame
  indexRange <- c(6)
#  significanceRange <- c(NA,signis)
  subFrame <- subFrame %>% filterMode(mode)  %>% filterIndexRange(indexRange) %>% filterSignificanceRange(SELECTED_SIGNIS) %>% filterMetricRange(fractionMetric)
  subFrame$significanceLevel <- unlist(lapply(subFrame$significanceLevel,function(x) {signiToPercent(x)}))

  subFrame <- subFrame %>% tidyr::pivot_wider(names_from = significanceLevel ,values_from = dplyr::all_of(fractionMetric))
  
  
  subFrame <- subFrame %>% dplyr::select(antibiotic ,`non-conformal`,`10%`,`5%`,`2.5%`)  
  
  subFrame
}



fetchAntibioticsVsMetric <- function(mode  = "Mode-A")
{
  frame <- readStatisticsExcel("abErrorStatisticsFrac") 
  subFrame <- frame
  subFrame <- subFrame %>% filterMode(mode) %>% filterSignificanceRange(NA) %>% reorderMetrics()
  
  indexRange <- c(4,6,8)
  aList <- lapply(indexRange, function(x) { subFrame %>% filterIndexRange(x)})
  names(aList) <-  as.character(indexRange)
  aList
}


countByAntibioticsGroup <- function(frame = readStatisticsExcel("abErrorStatisticsCount"))
{
  ab_groups <- tribble(
    ~antibiotic, ~ab_group,
    # Penicillins / β-lactam + inhibitor
    "AMP", "Penicillins",
    "AMC", "Penicillins",
    "PIP", "Penicillins",
    "TZP", "Penicillins",
    
    # Cephalosporins
    "CAZ", "Cephalosporins",
    "CRO", "Cephalosporins",
    "CTX", "Cephalosporins",
    "FEP", "Cephalosporins",
    
    # Fluoroquinolones
    "CIP", "Fluoroquinolones",
    "OFX", "Fluoroquinolones",
    "LVX", "Fluoroquinolones",
    "MFX", "Fluoroquinolones",
    
    # Aminoglycosides
    "GEN", "Aminoglycosides",
    "TOB", "Aminoglycosides"
  )

  frame_by_group <- frame %>%
    left_join(ab_groups, by = "antibiotic") %>%
    
    # build dynamic cols_to_sum
    { 
      cols_to_exclude <- intersect(
        c("noinputab", "significanceLevel", "antibiotic", "mode", "ab_group"),
        names(.)
      )
      cols_to_sum <- setdiff(names(.), cols_to_exclude)
      
      group_by(., noinputab, significanceLevel, ab_group, mode) %>%
        summarise(
          across(all_of(cols_to_sum), ~ sum(.x, na.rm = TRUE)),
          .groups = "drop"
        )
    }
  frame_by_group 
}






EXCEL_PLOTS <- function ()
{
  
  frameCorrectAntibioticsVsIndex <- fetchCorrectAntibioticsVsIndex(selectedMode)
  frameAntibioticsVsMetric <- fetchAntibioticsVsMetric(selectedMode)
  
  
  frameCorrectSignificanceLevelVsIndexA <- fetchCorrectSignificanceLevelVsIndex(selectedMode)
  p5A <- plot_siglevel_grouped(frameCorrectSignificanceLevelVsIndexA)
  p5A
  
  frameMetricVsSignificanceLevel <- fetchMetricVsSignificanceLevel(selectedMode)
  p5B <- plot_prediction_metrics_stacked(frameMetricVsSignificanceLevel)
  p5B
  
  p5_combined <- patchwork::wrap_plots(p5A, p5B, ncol = 2)
  
  ggplot2::ggsave(file.path(manuscriptPlotDirectory, "figure_5AB.png"),
                  p5_combined, width = 13, height = 5.2, dpi = 600, bg = "white")
  ggplot2::ggsave(file.path(manuscriptPlotDirectory, "figure_5AB.svg"),
                  p5_combined, width = 13, height = 5.2, bg = "white")
  p5_combined
  
  
  frameCorrectAntibioticsVsIndex <- fetchCorrectAntibioticsVsIndex(selectedMode)
  p7A <- plot_fig7A(frameCorrectAntibioticsVsIndex)
  p7A
  
  frameAntibioticsVsMetric <- fetchAntibioticsVsMetric(selectedMode)
  p7B <- plot_fig7B(frameAntibioticsVsMetric)
  p7B
  
  p7_combined <- patchwork::wrap_plots(p7A, p7B, ncol = 2)
  ggplot2::ggsave(file.path(manuscriptPlotDirectory, "figure_7AB.png"),
                  p7_combined, width = 13, height = 5.2, dpi = 600, bg = "white")
  ggplot2::ggsave(file.path(manuscriptPlotDirectory, "figure_7AB.svg"),
                  p7_combined, width = 13, height = 5.2, bg = "white")
  p7_combined

  
  
  
  metric_defs_S4 <- tibble::tribble(
    ~metric,     ~tag, ~ylab,
    "correct",   "A",  "correct",
    "F1",        "B",  "F1 score",
    "VME",       "C",  "VME",
    "ME",        "D",  "ME"
  )
  
  pS4 <- plot_figure_panels_with_callback(
    selectedMode,
    metric_defs = metric_defs_S4,
    fetch_fun   = fetchfractionMetricAntibioticsVsIndex,
    fig_stub    = "figure_S4",
    nrow        = 2
  )
  
  pS4
  
  pS4prim <- plot_figure_panels_with_callback(
    selectedMode,
    metric_defs = metric_defs_S4,
    fetch_fun   = fetchfractionMetricAntibioticsGroupVsIndex,
    fig_stub    = "figure_S4prim",
    nrow        = 2
  )
  
  pS4prim

    
  
  
  metric_defs_S5 <- tibble::tribble(
    ~metric,     ~tag, ~ylab,
    "correct",   "A",  "correct",
    "F1",        "B",  "F1 score",
    "VME",       "C",  "VME",
    "ME",        "D",  "ME"
  )
  
  pS5 <- plot_figure_panels_with_callback(
    selectedMode,
    metric_defs = metric_defs_S5,
    fetch_fun   = fetchfractionMetricAntibioticsVsSignicanceLevel,
    fig_stub    = "figure_S5",
    nrow        = 3
  )

  pS5
  
 
}


REF_PLOTS <- function()
{
  library(readxl)
  source("manuscript/metric_plot_helpers.R")
  
  # 1) Read the attached file
  df_counts <- readStatisticsExcel("abErrorStatisticsCount") 
  
  # 2) (Optionally) collapse antibiotics to groups
  df_group <- collapse_to_ab_groups(df_counts)
  
  # 3) Compute metrics – here with antibiotic as rows
  metrics_long_ab <- compute_metrics_long_counts(
    df_counts,
    group_vars = c("noinputab", "antibiotic", "significanceLevel", "mode")
  )
  
  # or: grouped to ab_group
  metrics_long_group <- compute_metrics_long_counts(
    df_group,
    group_vars = c("noinputab", "antibiotic", "significanceLevel", "mode")
  )
metrics_long_group_ref <- metrics_long_group

metrics_long_group_ref <- metrics_long_group_ref %>%
  tibble::add_row(
    noinputab = 6,
    antibiotic = "Cephalosporins",
    significanceLevel = NA,
    correctR = NA,
    falseR   = NA,
    correctS = NA,
    falseS   = NA,
    mode     = "ref",
    metric   = "F1",
    value    = 0.91
  ) %>%
  tibble::add_row(
    noinputab = 6,
    antibiotic = "Fluoroquinolones",
    significanceLevel = NA,
    correctR = NA,
    falseR   = NA,
    correctS = NA,
    falseS   = NA,
    mode     = "ref",
    metric   = "F1",
    value    = 0.86
  )  %>%
  tibble::add_row(
    noinputab = 6,
    antibiotic = "Penicillins",
    significanceLevel = NA,
    correctR = NA,
    falseR   = NA,
    correctS = NA,
    falseS   = NA,
    mode     = "ref",
    metric   = "F1",
    value    = 0.74
  ) %>%
  tibble::add_row(
    noinputab = 6,
    antibiotic = "Aminoglycosides",
    significanceLevel = NA,
    correctR = NA,
    falseR   = NA,
    correctS = NA,
    falseS   = NA,
    mode     = "ref",
    metric   = "F1",
    value    = 0.65
  ) %>%
  tibble::add_row(
    noinputab = 6,
    antibiotic = "Cephalosporins",
    significanceLevel = NA,
    correctR = NA,
    falseR   = NA,
    correctS = NA,
    falseS   = NA,
    mode     = "ref",
    metric   = "ME",
    value    = 0.014
  ) %>%
  tibble::add_row(
    noinputab = 6,
    antibiotic = "Fluoroquinolones",
    significanceLevel = NA,
    correctR = NA,
    falseR   = NA,
    correctS = NA,
    falseS   = NA,
    mode     = "ref",
    metric   = "ME",
    value    = 0.022
  )  %>%
  tibble::add_row(
    noinputab = 6,
    antibiotic = "Penicillins",
    significanceLevel = NA,
    correctR = NA,
    falseR   = NA,
    correctS = NA,
    falseS   = NA,
    mode     = "ref",
    metric   = "ME",
    value    = 0.125
  ) %>%
  tibble::add_row(
    noinputab = 6,
    antibiotic = "Aminoglycosides",
    significanceLevel = NA,
    correctR = NA,
    falseR   = NA,
    correctS = NA,
    falseS   = NA,
    mode     = "ref",
    metric   = "ME",
    value    = 0.073
  ) %>%
  tibble::add_row(
    noinputab = 6,
    antibiotic = "Cephalosporins",
    significanceLevel = NA,
    correctR = NA,
    falseR   = NA,
    correctS = NA,
    falseS   = NA,
    mode     = "ref",
    metric   = "VME",
    value    = 0.067
  ) %>%
  tibble::add_row(
    noinputab = 6,
    antibiotic = "Fluoroquinolones",
    significanceLevel = NA,
    correctR = NA,
    falseR   = NA,
    correctS = NA,
    falseS   = NA,
    mode     = "ref",
    metric   = "VME",
    value    = 0.192
  )  %>%
  tibble::add_row(
    noinputab = 6,
    antibiotic = "Penicillins",
    significanceLevel = NA,
    correctR = NA,
    falseR   = NA,
    correctS = NA,
    falseS   = NA,
    mode     = "ref",
    metric   = "VME",
    value    = 0.23
  ) %>%
  tibble::add_row(
    noinputab = 6,
    antibiotic = "Aminoglycosides",
    significanceLevel = NA,
    correctR = NA,
    falseR   = NA,
    correctS = NA,
    falseS   = NA,
    mode     = "ref",
    metric   = "VME",
    value    = 0.163
  ) 

metric_defs_S7 <- tibble::tribble(
  ~metric,   ~tag, ~ylab,
  "VME",     "A",  "VME",
  "ME",      "B",  "ME",
  "F1",      "C",  "F1 score"
)

pS7 <- plot_metric_panels(
  metrics_long = metrics_long_group_ref,
  metric_defs  = metric_defs_S7,
  id_col       = antibiotic,
  column_var   = mode,
  filter_expr  = is.na(significanceLevel) & noinputab == 6,
  fig_stub     = "figure_S7",
  nrow = 2,
  export       = TRUE
)
pS7
}


SAMPLE_LEVEL_PLOTS <- function()
{
  df_countsSample <- readStatisticsExcel("abSampleErrorStatisticsCount") 
  metrics_long_ab_sample <- compute_metrics_long_counts(
    df_countsSample,
    group_vars = c("noinputab", "antibiotic", "significanceLevel", "mode","sample")
  )
  
  df <- df_countsSample %>% 
    filter_and_drop(mode, "Mode-A") %>% 
    filter_and_drop(significanceLevel, NA) %>%    
    filter_and_drop(noinputab, 6) %>%
    summarise_over_and_drop(antibiotic)
  #VME
  sum(df$falseS)/(sum(df$correctR) + sum(df$falseS))
  #ME
  sum(df$falseR)/(sum(df$correctS) + sum(df$falseR))
  # fraction correct
  (sum(df$correctS)+sum(df$correctR))/(sum(df$correctS) + sum(df$falseR)+ sum(df$correctR) + sum(df$falseS))

  metrics_long_sample <- compute_metrics_long_counts(
    df,
    group_vars = c("sample")
  )


  p_S11 <- plot_metric_bars_simple(
    metrics_long = metrics_long_sample,
    metric       = "correct",
    id_col       = sample,   # x-axel = sample
    ylab         = "correct",
    export       = TRUE
  )
  p_S11
  
}

METRIC_PLOTS <- function()
{
  library(readxl)
  source("manuscript/metric_plot_helpers.R")

  # 1) Read the attached file
  df_counts <- readStatisticsExcel("abErrorStatisticsCount") 
  

  # 2) (Optionally) collapse antibiotics to groups
  df_group <- collapse_to_ab_groups(df_counts)
  
  # 3) Compute metrics – here with antibiotic as rows
  metrics_long_ab <- compute_metrics_long_counts(
    df_counts,
    group_vars = c("noinputab", "antibiotic", "significanceLevel", "mode")
  )

  
  
  
  # or: grouped to ab_group
  metrics_long_group <- compute_metrics_long_counts(
    df_group,
    group_vars = c("noinputab", "antibiotic", "significanceLevel", "mode")
  )
 

  # 5) Multi-panel figure (correct, F1, VME, ME), similar to your S4
  metric_defs_S4 <- tibble::tribble(
    ~metric,   ~tag, ~ylab,
    "correct", "A",  "Correct",
    "F1",      "B",  "F1 score",
    "VME",     "C",  "VME",
    "ME",      "D",  "ME"
  )
  
  pS4 <- plot_metric_panels(
    metrics_long = metrics_long_ab,
    metric_defs  = metric_defs_S4,
    id_col       = antibiotic,
    column_var   = noinputab,
    filter_expr  = is.na(significanceLevel) & mode == "Mode-A",
    fig_stub     = "figure_S4",
    nrow = 2,
    export       = TRUE
  )
  pS4
}

METRIC_PLOTS_BEST <-function()
{
  metric_defs_selected <- tibble::tribble(
    ~metric,   ~tag, ~ylab,
    "correct", "A",  "Correct",
    "F1",      "B",  "F1 score",
    "VME",     "C",  "VME",
    "ME",      "D",  "ME"
  )
  for(prefix in names(BEST_SELECTION)){
    df_counts <- readStatisticsExcel(paste(prefix,"abErrorStatisticsCount",sep="-"))
    metrics_long_ab <- compute_metrics_long_counts(
      df_counts,
    group_vars = c("noinputab", "antibiotic", "significanceLevel", "mode")
  )
  
  
  p <- plot_metric_panels(
    metrics_long = metrics_long_ab,
    metric_defs  = metric_defs_selected,
    id_col       = antibiotic,
    column_var   = noinputab,
    filter_expr  = is.na(significanceLevel) & mode == "Mode-A",
    fig_stub     = prefix,
    nrow = 2,
    export       = TRUE
  )
  p
  
  }
}




METRIC_PLOTS_ATLEASTONE <-function()
{
  metric_defs_selected <- tibble::tribble(
    ~metric,   ~tag, ~ylab,
    "correct", "A",  "Correct",
    "F1",      "B",  "F1 score",
    "VME",     "C",  "VME",
    "ME",      "D",  "ME"
  )
  df_counts_oneperabgroup <- readStatisticsExcel("fourbest-abErrorStatisticsCount")
  df_counts_oneperabgroup <- readStatisticsExcel("oneperabgroup-abErrorStatisticsCount")
  
  df_group_oneperabgroup <- collapse_to_ab_groups(df_counts_oneperabgroup)
  
  
  metrics_long_ab_oneperabgroup <- compute_metrics_long_counts(
    df_counts_oneperabgroup,
    group_vars = c("noinputab", "antibiotic", "significanceLevel", "mode")
  )
  
  metrics_long_group_oneperabgroup <- compute_metrics_long_counts(
    df_group_oneperabgroup,
    group_vars = c("noinputab", "antibiotic", "significanceLevel", "mode")
  )
  
  pS10 <- plot_metric_panels(
    metrics_long = metrics_long_ab_oneperabgroup,
    metric_defs  = metric_defs_selected,
    id_col       = antibiotic,
    column_var   = noinputab,
    filter_expr  = is.na(significanceLevel) & mode == "Mode-A",
    fig_stub     = "figure_S10",
    nrow = 2,
    export       = TRUE
  )
  pS10
}





CORRECT_VME_ME_PLOTS <- function()
{
  
  df_Frac <- countsToFracClinical(readStatisticsExcel("abErrorStatisticsCount"))
  metrics_long_ab_frac <- df_Frac  %>% 
    tidyr::pivot_longer(cols <- c("correct","ME","VME"),values_to = "value",names_to = "metric" )
  
  
  pS4A <- plot_stacked_errors(
    metrics_long = metrics_long_ab_frac,
    id_col       = antibiotic,
    column_var   = noinputab,
    value_col    = value,
    filter_expr  = is.na(significanceLevel) & mode == "Mode-A",
    export       = TRUE,
    file_stub    = "figure_S4A_stacked"
  )
  
  
  pS4A
}


PLOT_ESCMID_CONF <- function()
{
#  mFilterExpr <- quo(is.na(.data$significanceLevel) & .data$mode == "Mode-A")
#  mColumnVar <- quo(noinputab)
  mFilterExpr <- quo(!is.na(.data$significanceLevel) & .data$noinputab == 6 & .data$mode == "Mode-A")
  mColumnVar <- quo(significanceLevel)
  
  
  #Figure 1 is setup
  
  
  antibioticsDescription <-"\
  AMC-amoxicillin/clavulanic acid, \
  AMP-ampicillin, \
  CAZ-ceftazidime, \
  CIP-ciprofloxacin, \
  CRO-ceftriaxone, \
  CTX-cefotaxime, \
  FEP-cefepime, \
  GEN-gentamicin, \
  LVX-levofloxacin, \
  MFX-moxifloxacin, \
  OFX-ofloxacin \
  PIP-piperacillin, \
  TOB-tobramycin, \
  TZP-piperacillin-tazobactam."
  
  ME_VME_description <-"ME-major error, VME-very major error."
  
  fig2captionDescription <- str_wrap(
    "Figure 2 - Gross percentages of correctly and incorrectly predicted isolates using all combinations of 6 input antibiotics Using conformal predictions",
    width=80)
  
  fig2subcaptionDescription <- str_wrap(paste(
    "Proportion of the 99 isolates correctly predicted or misclassified as major errors or very major errors using conformal predictions.",
    antibioticsDescription,
    ME_VME_description,
    sep=" "),90)
  
  fig3captionDescription <- str_wrap(
    "Figure 3 - Antibiotic-specific major error and very major error rates across combinations of 6 input antibiotics using conformal predictions",
    width=80)
  
  fig3subcaptionDescription <- str_wrap(paste(
    "Prediction error rates for each antibiotic, showing false-resistant (major error) and false-susceptible (very major error) outcomes.",
    antibioticsDescription,
    ME_VME_description,
    sep=" "),90)
  
  
  # Figure 2 antibiotics
  df_Frac <- countsToFracClinical(readStatisticsExcel("abErrorStatisticsCount"))
  metrics_long_ab_frac <- df_Frac  %>% 
    tidyr::pivot_longer(cols <- c("correct","ME","VME"),values_to = "value",names_to = "metric" )
  
  # Figure 2 

  signiToPercent_vec <- Vectorize(signiToPercent)
  
  metrics_long_ab_frac_percent <- 
    metrics_long_ab_frac %>% 
    filter(!is.na(significanceLevel)) %>% 
    mutate(significanceLevel = signiToPercent_vec(significanceLevel))
  
  p2 <- plot_stacked_errors(
    metrics_long   = metrics_long_ab_frac_percent,
    id_col         = antibiotic,
    column_var     = !!mColumnVar, #"noinputab",  # bare, not quoted, since you use ensym()
    value_col      = value,
    filter_expr    = !!mFilterExpr,   # unquote the quosure here
    export         = TRUE,
    file_stub      = "ESCMID_figure_2prim",
    width          = 1000/150,
    height         = 800/150,
    dpi            = 150,
    fig_title_main = fig2captionDescription,
    fig_title_sub  = fig2subcaptionDescription
  )
  
  p2
  
  
  # Figure S2 is only used as support. Aggregation of all antibiotics 
  df_Frac_All <- countsToFracClinical(readStatisticsExcel("errorStatisticsCount"))
  metrics_long_frac <- df_Frac_All  %>%
    tidyr::pivot_longer(cols <- c("correct","ME","VME"),values_to = "value",names_to = "metric" )
  
  metrics_long_frac_percent <- 
    metrics_long_frac %>% 
    filter(!is.na(significanceLevel)) %>% 
    mutate(significanceLevel = signiToPercent_vec(significanceLevel))
  
  pS2 <- plot_stacked_errors_single(
    metrics_long = metrics_long_frac_percent,
    id_col       = !!mColumnVar,
    value_col    = value,
    filter_expr  = !!mFilterExpr,
    export       = TRUE,
    width        = 6,
    height       = 4,
    dpi          = 150,
    file_stub    = "ESCMID_figure_S2prim"
  )
  
  pS2
  
  
  #Figure 3  
  
  df_counts <- readStatisticsExcel("abErrorStatisticsCount") 
  metrics_long_ab <- compute_metrics_long_counts(
    df_counts,
    group_vars = c("noinputab", "antibiotic", "significanceLevel", "mode")
  )
  metric_defs_3 <- tibble::tribble(
    ~metric,   ~tag, ~ylab,
    "VME",     "A",  "VME rate",
    "ME",      "B",  "ME rate"
  )
  
  metrics_long_ab_percent <- 
    metrics_long_ab %>% 
    filter(!is.na(significanceLevel)) %>% 
    mutate(significanceLevel = signiToPercent_vec(significanceLevel))
  
  p3 <- plot_metric_panels(
    metrics_long = metrics_long_ab_percent,
    metric_defs  = metric_defs_3,
    id_col       = antibiotic,
    column_var   = !!mColumnVar,
    filter_expr  = !!mFilterExpr,
    fig_stub     = "ESCMID_figure_3prim",
    nrow = 2,
    export       = TRUE,
    width        = 1000/150,
    height       = 1000/150,
    dpi          = 150,
    fig_title_main = fig3captionDescription,
    fig_title_sub = fig3subcaptionDescription
  )
  p3
  
  # Figure S4 is when a selection of six well 
  #  
  df_FracBest <- countsToFracClinical(readStatisticsExcel(paste("sixbest","abErrorStatisticsCount",sep="-")))
  metrics_long_ab_frac_best <- df_FracBest  %>% 
    tidyr::pivot_longer(cols <- c("correct","ME","VME"),values_to = "value",names_to = "metric" )
  metrics_long_ab_frac_best_010 <- metrics_long_ab_frac_best %>% filter(significanceLevel == 0.1) 
    
  # Figure S4 is only used as support. Aggregation of all antibiotics.
  pS4 <- plot_stacked_errors_single(
    metrics_long = metrics_long_ab_frac_best_010,
    id_col       = antibiotic,
    value_col    = value,
    filter_expr  = !!mFilterExpr,
    export       = TRUE,
    file_stub    = "ESCMID_figure_S4prim",
  )
  pS4
  
}


PLOT_ESCMID <- function()
{
  #Figure 1 is setup
  
  
  antibioticsDescription <-"\
  AMC-amoxicillin/clavulanic acid, \
  AMP-ampicillin, \
  CAZ-ceftazidime, \
  CIP-ciprofloxacin, \
  CRO-ceftriaxone, \
  CTX-cefotaxime, \
  FEP-cefepime, \
  GEN-gentamicin, \
  LVX-levofloxacin, \
  MFX-moxifloxacin, \
  OFX-ofloxacin \
  PIP-piperacillin, \
  TOB-tobramycin, \
  TZP-piperacillin-tazobactam."

  ME_VME_description <-"ME-major error, VME-very major error."
  
  fig2captionDescription <- str_wrap(
    "Figure 2 - Gross percentages of correctly and incorrectly predicted isolates using all combinations of 4 to 8 input antibiotics",
    width=80)
  
  fig2subcaptionDescription <- str_wrap(paste(
    "Proportion of the 99 isolates correctly predicted or misclassified as major errors or very major errors.",
    antibioticsDescription,
    ME_VME_description,
    sep=" "),90)
  
  fig3captionDescription <- str_wrap(
    "Figure 3 - Antibiotic-specific major error and very major error rates across combinations of 4 to 8 input antibiotics",
    width=80)
  
  fig3subcaptionDescription <- str_wrap(paste(
    "Prediction error rates for each antibiotic, showing false-resistant (major error) and false-susceptible (very major error) outcomes.",
    antibioticsDescription,
    ME_VME_description,
    sep=" "),90)
    
  
  # Figure 2 antibiotics
  df_Frac <- countsToFracClinical(readStatisticsExcel("abErrorStatisticsCount"))
  metrics_long_ab_frac <- df_Frac  %>% 
    tidyr::pivot_longer(cols <- c("correct","ME","VME"),values_to = "value",names_to = "metric" )

  # Figure 2 
  p2 <- plot_stacked_errors(
    metrics_long = metrics_long_ab_frac,
    id_col       = antibiotic,
    column_var   = noinputab,
    value_col    = value,
    filter_expr  = is.na(significanceLevel) & mode == "Mode-A",
    export       = TRUE,
    file_stub    = "ESCMID_figure_2",
    width        = 1000/150,
    height       = 800/150,
    dpi          = 150,
    fig_title_main = fig2captionDescription,
    fig_title_sub = fig2subcaptionDescription
  )

  p2
  

  # Figure S2 is only used as support. Aggregation of all antibiotics 
  df_Frac_All <- countsToFracClinical(readStatisticsExcel("errorStatisticsCount"))
  metrics_long_frac <- df_Frac_All  %>%
    tidyr::pivot_longer(cols <- c("correct","ME","VME"),values_to = "value",names_to = "metric" )


  pS2 <- plot_stacked_errors_single(
    metrics_long = metrics_long_frac,
    id_col       = noinputab,
    value_col    = value,
    filter_expr  = is.na(significanceLevel) & mode == "Mode-A",
    export       = TRUE,
    width        = 6,
    height       = 4,
    dpi          = 150,
    file_stub    = "ESCMID_figure_S2"
  )

  pS2


  #Figure 3  
  
  df_counts <- readStatisticsExcel("abErrorStatisticsCount") 
  metrics_long_ab <- compute_metrics_long_counts(
    df_counts,
    group_vars = c("noinputab", "antibiotic", "significanceLevel", "mode")
  )
  metric_defs_3 <- tibble::tribble(
    ~metric,   ~tag, ~ylab,
    "VME",     "A",  "VME rate",
    "ME",      "B",  "ME rate"
  )
  
  p3 <- plot_metric_panels(
    metrics_long = metrics_long_ab,
    metric_defs  = metric_defs_3,
    id_col       = antibiotic,
    column_var   = noinputab,
    filter_expr  = is.na(significanceLevel) & mode == "Mode-A",
    fig_stub     = "ESCMID_figure_3",
    nrow = 2,
    export       = TRUE,
    width        = 1000/150,
    height       = 1000/150,
    dpi          = 150,
    fig_title_main = fig3captionDescription,
    fig_title_sub = fig3subcaptionDescription
  )
  p3

  # Figure S4 is when a selection of six well 
  #  
  df_FracBest <- countsToFracClinical(readStatisticsExcel(paste("sixbest","abErrorStatisticsCount",sep="-")))
  metrics_long_ab_frac_best <- df_FracBest  %>% 
    tidyr::pivot_longer(cols <- c("correct","ME","VME"),values_to = "value",names_to = "metric" )
  
  # Figure S4 is only used as support. Aggregation of all antibiotics.
  pS4 <- plot_stacked_errors_single(
    metrics_long = metrics_long_ab_frac_best,
    id_col       = antibiotic,
    value_col    = value,
    filter_expr  = is.na(significanceLevel) & mode == "Mode-A",
    export       = TRUE,
    file_stub    = "ESCMID_figure_S4",
  )
  pS4

}

ALL <- function()
{
  EXCEL_PLOTS()
  METRIC_PLOTS()
  METRIC_PLOTS_ATLEASTONE()
  METRIC_PLOTS_BEST()
  REF_PLOTS()
  PLOT_ESCMID()
  PLOT_ESCMID_CONF()
}


