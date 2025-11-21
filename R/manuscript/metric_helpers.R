library(dplyr)
library(rlang)

#GROUP_COLS <- c("noinputab", "antibiotic", "significanceLevel", "sample", "mode")

BASE_GROUP_COLS <- c("noinputab", "antibiotic", "significanceLevel", "sample", "mode")
FRAC_COLS <- c("correct","ME","VME")


library(dplyr)
library(rlang)

groupColumns <- function(df)
{
  intersect(BASE_GROUP_COLS, names(df))
}

fracAndGroupCoulmns <- function(df)
{
  c(intersect(BASE_GROUP_COLS, names(df)),FRAC_COLS)
}

filter_and_drop <- function(df, col, values) {
  col_quo <- enquo(col)
  
  # dela upp värden i NA och icke-NA
  values_no_na <- values[!is.na(values)]
  want_na      <- any(is.na(values))
  
  df %>%
    filter(
      (!!col_quo %in% values_no_na) | (want_na & is.na(!!col_quo))
    ) %>%
    select(-!!col_quo)
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


