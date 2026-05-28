source(file = 'manuscript/manuscriptcommon.R')
source(file = 'manuscript/metric_helpers_common.R')
library(dplyr)
library(tidyr)
library(readxl)


vmemePerformanceDataNonconformal <- function()
{
  # vme and me per antibiotic
  metrics_long_ab <- derivedMetricsFrame()
  df_base <- metrics_long_ab %>% 
    filter_and_drop(mode,MODE_A) %>% 
    filter_and_drop(cpmode,"normal") %>%
    filter_and_drop(significanceLevel,"STD") %>%
    filter_and_drop(noinputab,6) %>%
    filter(metric %in% c("ME","VME"))
  
  MEs <-  df_base %>% filter_and_drop(metric,"ME") %>% arrange(value)
  VMEs <-  df_base %>% filter_and_drop(metric,"VME") %>% arrange(value)
  print(MEs)
  print(VMEs)

  # overall vme and me
  metrics_long_ab <- derivedMetricsFrame(countFrameWide = readCountFrameWideSum())
  df_base <- metrics_long_ab %>% 
    filter_and_drop(mode,MODE_A) %>% 
    filter_and_drop(cpmode,"normal") %>%
    filter_and_drop(significanceLevel,"STD") %>%
    filter_and_drop(noinputab,6) %>%
    filter(metric %in% c("ME","VME"))
  
  print(df_base)

  # vme and me per v.s. input antibiotic
  metrics_long_ab <- derivedMetricsFrame(countFrameWide = readCountFrameWideSum())
  df_base <- metrics_long_ab %>% 
    filter_and_drop(mode,MODE_A) %>% 
    filter_and_drop(cpmode,"normal") %>%
    filter_and_drop(significanceLevel,"STD") %>%
    filter(metric %in% c("ME","VME"))
  
  print(df_base)

  # vme and me per v.s. input antibiotic
  metrics_long_ab <- derivedMetricsFrame(countFrameWide = readCountFrameWideSum())
  df_base <- metrics_long_ab %>% 
    filter_and_drop(mode,MODE_A) %>% 
    filter_and_drop(cpmode,"normal") %>%
    filter_and_drop(significanceLevel,"STD") %>%
    filter(metric %in% c("ME","VME"))
  
  print(df_base)
  
}

correctPerformanceDataNonconformal <- function()
{

  metrics_long_ab <- derivedMetricsFrame(countFrameWide = readCountFrameWideSum())
  df_base <- metrics_long_ab %>% 
    filter_and_drop(mode,MODE_A) %>% 
    filter_and_drop(cpmode,"normal") %>%
    filter_and_drop(significanceLevel,"STD") %>%
    filter(metric %in% c("correct"))
  
  print(df_base)
}

f1MMCPerformanceDataNonconformal <- function()
{
  metrics_long_ab <- derivedMetricsFrame()
  df_base <- metrics_long_ab %>% 
    filter_and_drop(mode,MODE_A) %>% 
    filter_and_drop(cpmode,"normal") %>%
    filter_and_drop(significanceLevel,"STD") %>%
    filter_and_drop(noinputab,6) %>%
    filter(metric %in% c("F1","MCC"))

  wide <- df_base %>% pivot_wider(values_from = value,names_from = metric)
  print(wide)
  
  metrics_long_ab <- derivedMetricsFrame(
    countFrameWide = readCountFrameWideSumByGroup(),
    vars = AB_GROUP_GROUPS_VAR
  )
  df_base <- metrics_long_ab %>% 
    filter_and_drop(mode,MODE_A) %>% 
    filter_and_drop(cpmode,"normal") %>%
    filter_and_drop(significanceLevel,"STD") %>%
    filter_and_drop(noinputab,6) %>%
    filter(metric %in% c("F1","MCC"))
  
  wide <- df_base %>% pivot_wider(values_from = value,names_from = metric)
  print(wide)
  
}

f1MMCPerformanceData90percent <- function()
{
  metrics_long_ab <- derivedMetricsFrame()
  df_base <- metrics_long_ab %>% 
    filter_and_drop(mode,MODE_A) %>% 
    filter_and_drop(cpmode,"normal") %>%
    filter_and_drop(significanceLevel,"10%") %>%
    filter_and_drop(noinputab,6) %>%
    filter(metric %in% c("MCC"))
  
  print(df_base)

  metrics_long_ab <- derivedMetricsFrame(countFrameWide = readCountFrameWideSum())
  df_base <- metrics_long_ab %>% 
    filter_and_drop(mode,MODE_A) %>% 
    filter_and_drop(cpmode,"normal") %>%
#    filter_and_drop(significanceLevel,"10%") %>%
    filter_and_drop(noinputab,6) %>%
    filter(metric %in% c("MCC"))
  
  print(df_base)
  
}


f1MMCPerformanceDataOnePerGroup <- function()
{
  metrics_long_ab_oneperabgroup <- derivedMetricsFrame(prefix = "oneperabgroup" )
  df_base_oneperabgroup <- metrics_long_ab_oneperabgroup %>% 
    filter_and_drop(mode,MODE_A) %>% 
    filter_and_drop(cpmode,"normal") %>%
    filter_and_drop(significanceLevel,"STD") %>%
    filter_and_drop(noinputab,6) %>%
    filter(metric %in% c("MCC"))
  
  metrics_long_ab <- derivedMetricsFrame()
  df_base <- metrics_long_ab %>% 
    filter_and_drop(mode,MODE_A) %>% 
    filter_and_drop(cpmode,"normal") %>%
    filter_and_drop(significanceLevel,"STD") %>%
    filter_and_drop(noinputab,6) %>%
    filter(metric %in% c("MCC"))

  df_combined <- df_base %>%
    rename(MCC_allcombinations = value) %>%
    inner_join(
      df_base_oneperabgroup %>% rename(MCC_oneperabgroup = value),
      by = c("antibiotic", "metric")
    ) %>%
    mutate(diff = MCC_oneperabgroup - MCC_allcombinations)
  
  
  print(df_combined)
  
  ############## Per antibiotics group ##############################
  
  metrics_long_ab <- derivedMetricsFrame(
    countFrameWide = readCountFrameWideSumByGroup(),
    vars = AB_GROUP_GROUPS_VAR
  )
  metrics_long_ab_oneperabgroup <- derivedMetricsFrame(prefix = "oneperabgroup",
    countFrameWide = readCountFrameWideSumByGroup(prefix = "oneperabgroup"),
    vars = AB_GROUP_GROUPS_VAR
  )

  df_base <- metrics_long_ab %>% 
    filter_and_drop(mode,MODE_A) %>% 
    filter_and_drop(cpmode,"normal") %>%
    filter_and_drop(significanceLevel,"STD") %>%
    filter_and_drop(noinputab,6) %>%
    filter(metric %in% c("MCC","VME","ME"))
  
  df_base_oneperabgroup <- metrics_long_ab_oneperabgroup %>% 
    filter_and_drop(mode,MODE_A) %>% 
    filter_and_drop(cpmode,"normal") %>%
    filter_and_drop(significanceLevel,"STD") %>%
    filter_and_drop(noinputab,6) %>%
    filter(metric %in% c("MCC","VME","ME"))

  df_combined <- df_base %>%
    rename(allcombinations = value) %>%
    inner_join(
      df_base_oneperabgroup %>% rename(oneperabgroup = value),
      by = c("ab_group", "metric")
    ) %>%
    mutate(diff = oneperabgroup - allcombinations)
  
  print(df_combined)
  ################## Overall ##############################################
  
  metrics_long_ab <- derivedMetricsFrame(
    countFrameWide = readCountFrameWideSum(),
    vars = AB_GROUP_GROUPS_VAR
  )
  metrics_long_ab_oneperabgroup <- derivedMetricsFrame(prefix = "oneperabgroup",
                                                       countFrameWide = readCountFrameWideSum(prefix = "oneperabgroup"),
                                                       vars = AB_GROUP_GROUPS_VAR
  )
  
  df_base <- metrics_long_ab %>% 
    filter_and_drop(mode,MODE_A) %>% 
    filter_and_drop(cpmode,"normal") %>%
    filter_and_drop(significanceLevel,"STD") %>%
    filter_and_drop(noinputab,6) %>%
    filter(metric %in% c("MCC","VME","ME"))
  
  df_base_oneperabgroup <- metrics_long_ab_oneperabgroup %>% 
    filter_and_drop(mode,MODE_A) %>% 
    filter_and_drop(cpmode,"normal") %>%
    filter_and_drop(significanceLevel,"STD") %>%
    filter_and_drop(noinputab,6) %>%
    filter(metric %in% c("MCC","VME","ME"))
  
  df_combined <- df_base %>%
    rename(allcombinations = value) %>%
    inner_join(
      df_base_oneperabgroup %>% rename(oneperabgroup = value),
      by = c("metric")
    ) %>%
    mutate(diff = oneperabgroup - allcombinations)
  
  print(df_combined)
}

optimalMCC <- function(){
  metrics_long_ab <- derivedMetricsFrame(
    countFrameWide = readCountWordFrameWide_AggregateOnAntibiotic(),
    vars = WORD_GROUPS_VAR
  )
  df_base <- metrics_long_ab %>% 
    filter_and_drop(mode,MODE_A) %>% 
    filter_and_drop(cpmode,"normal") %>%
    filter_and_drop(significanceLevel,"STD") %>%
    filter(metric %in% c("MCC")) 
  
  
  words <- list()
  for(i in c(4,6,8)){
    sorted <- df_base %>%
      filter_and_drop(noinputab,i) %>%
      arrange(desc(value))
    
    words[[as.character(i)]] <- sorted$word[1]
    print(sprintf("%s: %f" ,sorted$word[1], sorted$value[1]))
  }
  
  df_base <- metrics_long_ab %>% 
    filter_and_drop(mode,MODE_A) %>% 
    filter_and_drop(cpmode,"normal") %>%
    filter_and_drop(significanceLevel,"STD") %>%
    filter_and_drop(noinputab,6) %>%
    filter_and_drop(word,words$`6`) %>%
    filter(metric %in% c("MCC","ME","VME")) 
  print(df_base)
}

translationModePerformanceDataNonconformal <- function()
{

  metrics_long_ab <- derivedMetricsFrame(
    countFrameWide = readCountFrameWideSum(),
    vars = ALL_GROUPS_VAR
  )
  df_base <- metrics_long_ab %>% 
    filter_and_drop(cpmode,"normal") %>%
    filter_and_drop(significanceLevel,"STD") %>%
    filter_and_drop(noinputab,6) %>%
    filter(metric %in% c("MCC","ME","VME")) 
  
  
  wide <- df_base %>% pivot_wider(values_from = value,names_from = metric)
  print(wide)
  
}

noPredictions <- function()
{
  
  metrics_long_ab <- derivedMetricsFrameUnambiguous(countFrameWide = readCountFrameWideSum())
  df_base <- metrics_long_ab %>% 
    filter_and_drop(mode,MODE_A) %>% 
    filter_and_drop(cpmode,"normal") %>%
    filter_and_drop(noinputab,6) %>%
    filter(!significanceLevel %in% "STD")
    
  nopredictions <- df_base %>% pivot_wider(values_from =  value,names_from = metric) %>% mutate (notPredicted = 1-correct-ME-VME)  
    
  print(nopredictions)
}








ALL <- function()
{
  performance_log_calls <- list(
    vmemePerformanceDataNonconformal = function() {
      vmemePerformanceDataNonconformal()
    },
    
    correctPerformanceDataNonconformal = function() {
      correctPerformanceDataNonconformal()
    },
    
    f1MMCPerformanceDataNonconformal = function() {
      f1MMCPerformanceDataNonconformal()
    },
    
    f1MMCPerformanceData90percent = function() {
      f1MMCPerformanceData90percent()
    },
    
    f1MMCPerformanceDataOnePerGroup = function() {
      f1MMCPerformanceDataOnePerGroup()
    },
    
    optimalMCC = function() {
      optimalMCC()
    },
    
    translationModePerformanceDataNonconformal = function() {
      translationModePerformanceDataNonconformal()
    },
    
    noPredictions = function() {
      noPredictions()
    }
  )
  
  run_with_log(
    calls = performance_log_calls,
    log_file = paste(manuscriptDirectory,"manuscript_data_log.txt",sep = "/")
  )
}
