source(file = 'manuscript/manuscriptcommon.R')
source(file = 'model/modelcommon.R')

# FIXME remove
#MODES <- list(MODE_A,MODE_B,MODE_C)
#MODES <- list(MODE_A)

SELECTED_SIGNIS <- c(NA,0.1,0.05,0.025)

mergeRiskStratifiedStatistics <-function(prefix=NA,subfolder = "")
{
  names <- c("RStatistics",
             "SStatistics",
             "RStatisticsAb",
             "SStatisticsAb",
             "RStatisticsAbSample",
             "SStatisticsAbSample")
  for(name in names){
    orgName <- name
    if(!is.na(prefix)){
      name <- paste(prefix,name,sep="-")
    }
    framesAllModes <- lapply(MODES,function(mode){readStatisticsFrameCSV(name,subfolder,mode)})
    merged <- NULL
    for(i in 1:length(MODES)){
      mode <- MODES[[i]]
      f <- framesAllModes[[i]]
      addedMode <- f %>% dplyr::mutate(mode = mode)
      if(is.null(merged)){
        merged <- addedMode
      } else {
        merged <- rbind(merged,addedMode)
      }
    }
    writeStatisticsExcel(merged,name,orgName)
  }  
}




mergeModesErrorStatistics <- function(prefix=NA,subfolder = "")
{
  #prefix <- "AMP_CRO_CIP_TOB"
  names <- c("abErrorStatisticsFrac",
             "errorStatisticsFrac","abSampleErrorStatisticsFrac",
             "errorStatisticsCount","abErrorStatisticsCount",
             "abSampleErrorStatisticsCount")
  
  for(name in names){
    orgName <- name
    if(!is.na(prefix)){
      name <- paste(prefix,name,sep="-")
    }
    framesAllModes <- lapply(MODES,function(mode){readStatisticsFrameCSV(name,subfolder,mode)})
    merged <- NULL
    for(i in 1:length(MODES)){
      mode <- MODES[[i]]
      f <- framesAllModes[[i]]
      addedMode <- f %>% dplyr::mutate(mode = mode)
      if(is.null(merged)){
        merged <- addedMode
      } else {
        merged <- rbind(merged,addedMode)
      }
    }
    writeStatisticsExcel(merged,name,orgName,subfolder=subfolder)
  }  
}



filterSignificanceRange <-  function(frame,significanceRange=NA)
{
  myIsNa <- function(x){!length(x)>1 && is.na(x)}
  
  #frame <- readStatisticsExcel(name)
  subFrame <- frame 
  #significanceRange <- c(NA,0.1)
  if(!myIsNa(significanceRange)){
    subFrame <- subFrame %>% dplyr::filter(significanceLevel %in% significanceRange)
  } else {
    subFrame <- subFrame %>% dplyr::filter(is.na(significanceLevel))
  }
  if(length(unique(subFrame$significanceLevel))==1){
    subFrame <- subFrame %>% dplyr::select(!significanceLevel) 
  }
  subFrame
}

filterIndexRange <- function(frame,indexRange){
  #frame <- readStatisticsExcel(name)
  subFrame <- frame 
  #indexRange <- c(4:8)
  subFrame <- subFrame %>% dplyr::filter(noinputab %in% indexRange)
  if(length(unique(subFrame$noinputab))==1){
    subFrame <- subFrame %>% dplyr::select(!noinputab) 
  }
  subFrame
}




filterMetricRangeNorm <- function(frame, metricRange) {
  filterMetricRange(frame, metricRange)
}

filterMetricRange <- function(frame, metricRange) {
  library(dplyr)
  
  # Normalize metricRange to a character vector
  if (is.null(metricRange)) {
    stop("metricRange must be provided (scalar or list/vector).")
  }
  if (is.list(metricRange)) {
    metricRange <- unlist(metricRange, recursive = TRUE, use.names = FALSE)
  }
  metricRange <- as.character(metricRange)
  metricRange <- unique(metricRange)
  
  if (length(metricRange) == 0) {
    stop("metricRange is empty after normalization.")
  }
  
  # Optional identifier columns (keep if present)
  first_list <- c("noinputab", "significanceLevel", "antibiotic","ab_group", "mode")
  id_cols <- intersect(first_list, names(frame))
  
  # All metric columns must be present
  missing_metrics <- setdiff(metricRange, names(frame))
  if (length(missing_metrics) > 0) {
    stop(
      paste("Missing required metric columns:",
            paste(missing_metrics, collapse = ", "))
    )
  }
  # Return same rows, keeping (optional) IDs + required metrics
  frame %>%
    dplyr::select(all_of(c(id_cols, metricRange)))
}  






filterMode <- function(frame,modeRange){
  #frame <- readStatisticsExcel(name)
  subFrame <- frame 
  #modeRange <- c("ME")
  subFrame <- subFrame %>% dplyr::filter(mode %in% modeRange)
  if(length(unique(subFrame$mode))==1){
    subFrame <- subFrame %>% dplyr::select(!mode) 
  }
  subFrame
}

antibioticsVsMode <- function(index  = 6)
{
  frame <- readStatisticsExcel("abErrorStatisticsFrac")
  subframe <- frame
  subframe <- subframe %>% filterIndexRange(index) %>% filterSignificanceRange(NA)
  correct <- subframe %>% filterMetricRange("correct") %>% tidyr::pivot_wider(names_from = mode,values_from = correct)
  VME <- subframe %>% filterMetricRange("VME") %>% tidyr::pivot_wider(names_from = mode,values_from = VME)
  ME <- subframe %>% filterMetricRange("ME") %>% tidyr::pivot_wider(names_from = mode,values_from = ME)
  aList <- list(correct=correct,VME=VME,ME=ME)  
  writeStatisticsExcel(aList,paste("antibioticsVsMode",index,sep="-"),paste("abMetricMode",index,sep="-"))
}


antibioticsVsMetric <- function(mode  = "Mode-A")
{
  frame <- readStatisticsExcel("abErrorStatisticsFrac") 
  subFrame <- frame
  subFrame <- subFrame %>% filterMode(mode) %>% filterSignificanceRange(NA) %>% reorderMetrics()
  
  indexRange <- c(4,6,8)
  aList <- lapply(indexRange, function(x) { subFrame %>% filterIndexRange(x)})
  names(aList) <-  as.character(indexRange)
  writeStatisticsExcel(aList,paste("antibioticsVsMetric",mode,sep="-"),paste("abVsMetric",substring(mode,nchar(mode),nchar(mode)),sep="-"))
}


correctAntibioticsVsIndex <- function(mode  = "Mode-A")
{
  frame <- readStatisticsExcel("abErrorStatisticsFrac") 
  subFrame <- frame
  indexRange <- c(4:8)
  subFrame <- subFrame %>% filterMode(mode) %>% filterSignificanceRange(NA) %>% filterIndexRange(indexRange) %>% filterMetricRange("correct")
  subFrame <- subFrame %>% tidyr::pivot_wider(names_from = noinputab,values_from = correct)

  writeStatisticsExcel(subFrame,paste("correctAntibioticsVsIndex",mode,sep="-"),paste("correctAntibioticsVsIndex",substring(mode,nchar(mode),nchar(mode)),sep="-"))
}

reorderMetrics <-function(frame)
{
  longSubframe <- frame %>% tidyr::pivot_longer(cols = METRICS_COLS, names_to = "metric")
  longSubframe %>% tidyr::pivot_wider(names_from = metric ,values_from = value)
}



signiToPercent <- function(signi)
{
  if(is.null(signi) || is.na(signi)){
    rv <- "non-conformal"
  } else {
    rv <- format_percent(signi)
  }
  rv
}


correctSignificanceLevelVsIndex <- function(mode  = MODE_A)
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
  writeStatisticsExcel(subFrame,paste("correctSignificanceLevelVsIndex",mode,sep="-"),paste("correctSignilVsIndex",substring(mode,nchar(mode),nchar(mode)),sep="-"))
}

correctModeVsSignificanceLevel <- function(index  = 6)
{
  frame <- readStatisticsExcel("errorStatisticsFrac")
  subFrame <- frame
  subFrame <- subFrame %>% 
    filterIndexRange(index) %>% 
    filterMetricRange("correct") %>% 
    filterSignificanceRange(SELECTED_SIGNIS)
  subFrame$significanceLevel <- unlist(lapply(subFrame$significanceLevel,function(x) {signiToPercent(x)}))
  
  subFrame <- subFrame %>% tidyr::pivot_wider(names_from = significanceLevel ,values_from = correct)

  subFrame <- subFrame %>% dplyr::select(mode ,`10%`,`5%`,`2.5%`)  

  writeStatisticsExcel(subFrame,paste("correctModeVsSignificanceLevel",index,sep="-"),paste("correctModeVsSigni",index,sep="-"))
}



meMetricVsAntibiotic <- function(index = 6,significanceLevel = NA)
{
  name <-  "SStatisticsAb"
  
  frame <- readStatisticsExcel(name)
  subFrame <- frame
  subFrame <- subFrame %>% 
    filterIndexRange(index) %>% 
    filterSignificanceRange(significanceLevel) %>%
    filterMode(modeRange = MODE_A)
  #Swap
  subFrame <- subFrame %>% tidyr::pivot_longer(cols = c("correct","ME"),names_to = "metric")  %>%
        tidyr::pivot_wider(names_from = antibiotic, values_from = value)
  name <- paste("meMetricVsAntibiotic",index,significanceLevel,sep="-")
  writeStatisticsExcel(subFrame,name,name)
  
}

vmeRateMetricVsAntibiotic <- function(index = 6,significanceLevel = NA)
{
  name <-  "RStatisticsAb"
  
  frame <- readStatisticsExcel(name)
  subFrame <- frame
  subFrame <- subFrame %>% 
    filterIndexRange(index) %>% 
    filterSignificanceRange(significanceLevel) %>%
    filterMode(modeRange = MODE_A)
  #Swap
  subFrame <- subFrame %>% tidyr::pivot_longer(cols = c("correct","VME"),names_to = "metric")  %>%
    tidyr::pivot_wider(names_from = antibiotic, values_from = value)
  name <- paste("vmeRateMetricVsAntibiotic",index,significanceLevel,sep="-")
  writeStatisticsExcel(subFrame,name,name)
}


correctModeVsIndex <- function(significanceLevel  = 0.1)
{
  frame <- readStatisticsExcel("errorStatisticsFrac")
  subFrame <- frame
  subFrame <- subFrame %>% 
    filterIndexRange(c(4:8)) %>%
    filterMetricRange("correct") %>% 
    filterSignificanceRange(significanceLevel)

  subFrame <- subFrame %>% tidyr::pivot_wider(names_from = noinputab ,values_from = correct)
  writeStatisticsExcel(subFrame,paste("correctModeVsIndex",significanceLevel,sep="-"),paste("correctModeVsIndex",significanceLevel,sep="-"))
}



metricVsSignificanceLevel <- function(mode  = MODE_A)
{
  indexRange <- c(4,6,8)
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
  
  writeStatisticsExcel(aList,paste("metricVsSignificanceLevel",mode,sep="-"),paste("metricVsSigni",substring(mode,nchar(mode),nchar(mode)),sep="-"))
}




vmeRateModeVsSignificanceLevel <- function(index  = 6,prefix=NA,subfolder="")
{
  name <-  "RStatistics"
  if(!is.na(prefix)){
    name <- paste(prefix,name,sep="-")
  }
  frame <- readStatisticsExcel(name,subfolder=subfolder)
  subFrame <- frame
  subFrame <- subFrame %>% 
    filterIndexRange(index) %>% 
    dplyr::select(-correct) %>%
    filterSignificanceRange(SELECTED_SIGNIS)
  subFrame$significanceLevel <- unlist(lapply(subFrame$significanceLevel,function(x) {signiToPercent(x)}))
  
  subFrame <- subFrame %>% tidyr::pivot_wider(names_from = significanceLevel ,values_from = VME)
  
  subFrame <- subFrame %>% dplyr::select(mode ,`10%`,`5%`,`2.5%`)  
  
  outname <- "vmeRateModeVsSignificanceLevel"
  if(!is.na(prefix)){
    outname <- paste(prefix,outname,sep="-")
  } else {
    outname <- paste(outname,index,sep="-")
  }
  writeStatisticsExcel(subFrame,outname,paste("vmeRateModeVsSigni",index,sep="-"),subfolder=subfolder)
}

meRateModeVsSignificanceLevel <- function(index  = 6,prefix=NA,subfolder="")
{
  name <-  "SStatistics"
  if(!is.na(prefix)){
    name <- paste(prefix,name,sep="-")
  }
  frame <- readStatisticsExcel(name,subfolder=subfolder)
  subFrame <- frame
  subFrame <- subFrame %>% 
    filterIndexRange(index) %>% 
    dplyr::select(-correct) %>%
    filterSignificanceRange(SELECTED_SIGNIS)
  subFrame$significanceLevel <- unlist(lapply(subFrame$significanceLevel,function(x) {signiToPercent(x)}))
  
  subFrame <- subFrame %>% tidyr::pivot_wider(names_from = significanceLevel ,values_from = ME)
  
  subFrame <- subFrame %>% dplyr::select(mode ,`10%`,`5%`,`2.5%`)  
  
  outname <- "meRateModeVsSignificanceLevel"
  if(!is.na(prefix)){
    outname <- paste(prefix,outname,sep="-")
  } else {
    outname <- paste(outname,index,sep="-")
  }
  writeStatisticsExcel(subFrame,outname,paste("meRateModeVsSigni",index,sep="-"),subfolder=subfolder)
}

meRateAntibioticsVsSignificanceLevel <- function(index  = 6,prefix=NA,mode = MODE_A,subfolder = "")
{
  name <-  "SStatisticsAb"
  if(!is.na(prefix)){
    name <- paste(prefix,name,sep="-")
  }
  frame <- readStatisticsExcel(name)
  subFrame <- frame
  subFrame <- subFrame %>% 
    filterMode(mode) %>% 
    filterIndexRange(index) %>% 
    dplyr::select(-correct) %>% 
    filterSignificanceRange(SELECTED_SIGNIS)
  subFrame$significanceLevel <- unlist(lapply(subFrame$significanceLevel,function(x) {signiToPercent(x)}))
  
  subFrame <- subFrame %>% tidyr::pivot_wider(names_from = significanceLevel ,values_from = ME)
  subFrame <- subFrame %>% dplyr::select(antibiotic,`non-conformal` ,`10%`,`5%`,`2.5%`)  
  outname <- "meRateAntibioticsVsSignificanceLevel"
  if(!is.na(prefix)){
    outname <- paste(prefix,outname,index,sep="-")
  } else {
    outname <- paste(outname,index,sep="-")
  }
  
  writeStatisticsExcel(subFrame,outname,paste("meRateAntibioticsVsSigni",index,sep="-"),subfolder=subfolder)
}





meRateAntibioticsVsIndex <- function(prefix=NA,subfolder = "",indexRange = c(4,6,8),mode = MODE_A)
{
  name <-  "SStatisticsAb"
  if(!is.na(prefix)){
    name <- paste(prefix,name,sep="-")
  }
  frame <- readStatisticsExcel(name)
  subFrame <- frame
  subFrame <- subFrame %>% 
    filterMode(mode) %>% 
    dplyr::select(-correct) %>% 
    filterSignificanceRange(NA) %>% 
    filterIndexRange(indexRange = indexRange)

  subFrame <- subFrame %>% tidyr::pivot_wider(names_from = noinputab ,values_from = ME)
  outname <- "meRateAntibioticsVsIndex"
  if(!is.na(prefix)){
    outname <- paste(prefix,outname,paste(indexRange,collapse = "_"),mode,sep="-")
  } else {
    outname <- paste(outname,paste(indexRange,collapse = "_"),mode,sep="-")
  }
  
  writeStatisticsExcel(subFrame,outname,paste("meRateAntibioticsVsIndex",sep="-"),subfolder=subfolder)
}

vmeRateAntibioticsVsIndex <- function(prefix=NA,subfolder = "",indexRange = c(4,6,8),mode = MODE_A)
{
  name <-  "RStatisticsAb"
  if(!is.na(prefix)){
    name <- paste(prefix,name,sep="-")
  }
  frame <- readStatisticsExcel(name)
  subFrame <- frame
  subFrame <- subFrame %>% 
    filterMode(mode) %>% 
    dplyr::select(-correct) %>% 
    filterSignificanceRange(NA) %>% 
    filterIndexRange(indexRange = indexRange)
  
  subFrame <- subFrame %>% tidyr::pivot_wider(names_from = noinputab ,values_from = VME)
  outname <- "vmeRateAntibioticsVsIndex"
  if(!is.na(prefix)){
    outname <- paste(prefix,outname,paste(indexRange,collapse = "_"),mode,sep="-")
  } else {
    outname <- paste(outname,paste(indexRange,collapse = "_"),mode,sep="-")
  }
  
  writeStatisticsExcel(subFrame,outname,paste("vmeRateAntibioticsVsIndex",sep="-"),subfolder=subfolder)
}



vmeRateAntibioticsVsSignificanceLevel <- function(index  = 6,prefix=NA,mode = MODE_A,subfolder = "")
{
  name <-  "RStatisticsAb"
  if(!is.na(prefix)){
    name <- paste(prefix,name,sep="-")
  }
  frame <- readStatisticsExcel(name)
  subFrame <- frame
  subFrame <- subFrame %>% 
    filterMode(mode) %>% 
    filterIndexRange(index) %>% 
    dplyr::select(-correct) %>% 
    filterSignificanceRange(SELECTED_SIGNIS)
  subFrame$significanceLevel <- unlist(lapply(subFrame$significanceLevel,function(x) {signiToPercent(x)}))
  
  subFrame <- subFrame %>% tidyr::pivot_wider(names_from = significanceLevel ,values_from = VME)
  subFrame <- subFrame %>% dplyr::select(antibiotic,`non-conformal` ,`10%`,`5%`,`2.5%`)  
  outname <- "vmeRateAntibioticsVsSignificanceLevel"
  if(!is.na(prefix)){
    outname <- paste(prefix,outname,index,sep="-")
  } else {
    outname <- paste(outname,index,sep="-")
  }
  
  writeStatisticsExcel(subFrame,outname,paste("vmeRateAntibioticsVsSigni",index,sep="-"),subfolder=subfolder)
}



correctRAntibioticsVsSignificanceLevel <- function(index  = 6,prefix=NA,subfolder = "")
{
  name <-  "RStatisticsAb"
  if(!is.na(prefix)){
    name <- paste(prefix,name,sep="-")
  }
  frame <- readStatisticsExcel(name)
  
  subFrame <- frame
  subFrame <- subFrame %>% 
    filterMode(MODE_A) %>% 
    filterIndexRange(index) %>% 
    dplyr::select(-VME) %>% 
    filterSignificanceRange(SELECTED_SIGNIS)
  subFrame$significanceLevel <- unlist(lapply(subFrame$significanceLevel,function(x) {signiToPercent(x)}))
  
  subFrame <- subFrame %>% tidyr::pivot_wider(names_from = significanceLevel ,values_from = correct)
  subFrame <- subFrame %>% dplyr::select(antibiotic,`non-conformal` ,`10%`,`5%`,`2.5%`)  
  outname <- "correctRAntibioticsVsSignificanceLevel"
  if(!is.na(prefix)){
    outname <- paste(prefix,outname,index,sep="-")
  } else {
    outname <- paste(outname,index,sep="-")
  }
  
  writeStatisticsExcel(subFrame,outname,paste("correctRAntibioticsVsSigni",index,sep="-"),subfolder=subfolder)
}

correctSAntibioticsVsSignificanceLevel <- function(index  = 6,prefix=NA,subfolder = "")
{
  name <-  "SStatisticsAb"
  if(!is.na(prefix)){
    name <- paste(prefix,name,sep="-")
  }
  frame <- readStatisticsExcel(name)
  
  subFrame <- frame
  subFrame <- subFrame %>% 
    filterMode(MODE_A) %>% 
    filterIndexRange(index) %>% 
    dplyr::select(-ME) %>% 
    filterSignificanceRange(SELECTED_SIGNIS)
  subFrame$significanceLevel <- unlist(lapply(subFrame$significanceLevel,function(x) {signiToPercent(x)}))
  
  subFrame <- subFrame %>% tidyr::pivot_wider(names_from = significanceLevel ,values_from = correct)
  subFrame <- subFrame %>% dplyr::select(antibiotic,`non-conformal` ,`10%`,`5%`,`2.5%`)  
  outname <- "correctSAntibioticsVsSignificanceLevel"
  if(!is.na(prefix)){
    outname <- paste(prefix,outname,index,sep="-")
  } else {
    outname <- paste(outname,index,sep="-")
  }
  
  writeStatisticsExcel(subFrame,outname,paste("correctSAntibioticsVsSigni",index,sep="-"),subfolder=subfolder)
}

correctSAntibioticsVsMode <- function(index  = 6,prefix=NA,subfolder = "")
{
  name <-  "SStatisticsAb"
  if(!is.na(prefix)){
    name <- paste(prefix,name,sep="-")
  }
  frame <- readStatisticsExcel(name)
  
  subFrame <- frame
  subFrame <- subFrame %>% 
    filterSignificanceRange() %>% 
    filterIndexRange(index) %>% 
    dplyr::select(-ME)

  subFrame <- subFrame %>% tidyr::pivot_wider(names_from = mode ,values_from = correct)
  outname <- "correctSAntibioticsVsMode"
  if(!is.na(prefix)){
    outname <- paste(prefix,outname,index,sep="-")
  } else {
    outname <- paste(outname,index,sep="-")
  }
  
  writeStatisticsExcel(subFrame,outname,paste("correctSAntibioticsVsMode",index,sep="-"),subfolder=subfolder)
}

meRateAntibioticsVsMode <- function(index  = 6,prefix=NA,subfolder = "")
{
  name <-  "SStatisticsAb"
  if(!is.na(prefix)){
    name <- paste(prefix,name,sep="-")
  }
  frame <- readStatisticsExcel(name)
  
  subFrame <- frame
  subFrame <- subFrame %>% 
    filterSignificanceRange() %>% 
    filterIndexRange(index) %>% 
    dplyr::select(-correct)
  
  subFrame <- subFrame %>% tidyr::pivot_wider(names_from = mode ,values_from = ME)
  outname <- "meRateAntibioticsVsMode"
  if(!is.na(prefix)){
    outname <- paste(prefix,outname,index,sep="-")
  } else {
    outname <- paste(outname,index,sep="-")
  }
  
  writeStatisticsExcel(subFrame,outname,paste("meRateAntibioticsVsMode",index,sep="-"),subfolder=subfolder)
}


correctRAntibioticsVsMode <- function(index  = 6,prefix=NA,subfolder = "")
{
  name <-  "RStatisticsAb"
  if(!is.na(prefix)){
    name <- paste(prefix,name,sep="-")
  }
  frame <- readStatisticsExcel(name)
  
  subFrame <- frame
  subFrame <- subFrame %>% 
    filterSignificanceRange() %>% 
    filterIndexRange(index) %>% 
    dplyr::select(-VME)
  
  subFrame <- subFrame %>% tidyr::pivot_wider(names_from = mode ,values_from = correct)
  outname <- "correctRAntibioticsVsMode"
  if(!is.na(prefix)){
    outname <- paste(prefix,outname,index,sep="-")
  } else {
    outname <- paste(outname,index,sep="-")
  }
  
  writeStatisticsExcel(subFrame,outname,paste("correctRAntibioticsVsMode",index,sep="-"),subfolder=subfolder)
}

vmeRateAntibioticsVsMode <- function(index  = 6,prefix=NA,subfolder = "")
{
  name <-  "RStatisticsAb"
  if(!is.na(prefix)){
    name <- paste(prefix,name,sep="-")
  }
  frame <- readStatisticsExcel(name)
  
  subFrame <- frame
  subFrame <- subFrame %>% 
    filterSignificanceRange() %>% 
    filterIndexRange(index) %>% 
    dplyr::select(-correct)
  
  subFrame <- subFrame %>% tidyr::pivot_wider(names_from = mode ,values_from = VME)
  outname <- "vmeRateAntibioticsVsMode"
  if(!is.na(prefix)){
    outname <- paste(prefix,outname,index,sep="-")
  } else {
    outname <- paste(outname,index,sep="-")
  }
  
  writeStatisticsExcel(subFrame,outname,paste("vmeRateAntibioticsVsMode",index,sep="-"),subfolder=subfolder)
}

vmeRateSampleVsAntibiotic <- function(index  = 6,prefix=NA,subfolder = "")
{
  name <-  "RStatisticsAbSample"
  if(!is.na(prefix)){
    name <- paste(prefix,name,sep="-")
  }
  frame <- readStatisticsExcel(name)
  
  subFrame <- frame
  subFrame <- subFrame %>% 
    filterSignificanceRange() %>%
    filterMode(MODE_A) %>% 
    filterIndexRange(index) %>% 
    dplyr::select(-correct)
  
  subFrame <- subFrame %>% tidyr::pivot_wider(names_from = antibiotic ,values_from = VME)
  outname <- "vmeRateSampleVsAntibiotic"
  if(!is.na(prefix)){
    outname <- paste(prefix,outname,index,sep="-")
  } else {
    outname <- paste(outname,index,sep="-")
  }
  
  writeStatisticsExcel(subFrame,outname,paste("vmeRateSampleVsAntibiotic",index,sep="-"),subfolder=subfolder)
}

meRateSampleVsAntibiotic <- function(index  = 6,prefix=NA,subfolder = "")
{
  name <-  "SStatisticsAbSample"
  if(!is.na(prefix)){
    name <- paste(prefix,name,sep="-")
  }
  frame <- readStatisticsExcel(name)
  
  subFrame <- frame
  subFrame <- subFrame %>% 
    filterSignificanceRange() %>%
    filterMode(MODE_A) %>% 
    filterIndexRange(index) %>% 
    dplyr::select(-correct)
  
  subFrame <- subFrame %>% tidyr::pivot_wider(names_from = antibiotic ,values_from = ME)
  outname <- "meRateSampleVsAntibiotic"
  if(!is.na(prefix)){
    outname <- paste(prefix,outname,index,sep="-")
  } else {
    outname <- paste(outname,index,sep="-")
  }
  
  writeStatisticsExcel(subFrame,outname,paste("meRateSampleVsAntibiotic",index,sep="-"),subfolder=subfolder)
}


sanitycheck <- function()
{
   fullFrame <- readStatisticsFrameCSV(name="abErrorStatisticsFull",mode="Mode-A")
   subFrame <- fullFrame %>% dplyr::filter(noinputab==1 & significanceLevel==0.1) %>% dplyr::select(correct, correctS, correctR, false, falseS, falseR, total,    S  , R, ambiguous, ambiguousS ,ambiguousR ,notpredicted ,notpredictedS, notpredictedR ,predicted)
   aSummary <- apply(X=subFrame,MARGIN=c(2),FUN=sum)
   VMErate <- (aSummary["falseS"] + aSummary["notpredictedR"])/aSummary["R"]
   MErate <- (aSummary["falseR"] + aSummary["notpredictedS"])/aSummary["S"]
   meRateModeVsSignificanceLevel(index=1)
   vmeRateModeVsSignificanceLevel(index=1)
   #(falseS+notpredictedR)/R
}

ALL <- function ()
{
  #merge data from different modes
  for(aPrefix in c(NA,SUBGROUPS)){
    mergeModesErrorStatistics(prefix = aPrefix)
    mergeRiskStratifiedStatistics(prefix = aPrefix)
  }
    
  #Figure 5A
  correctSignificanceLevelVsIndex(MODE_A)
  #Figure 5B
  metricVsSignificanceLevel(MODE_A)

  #Fig 6A
  correctModeVsIndex()
  #Fig 6B
  correctModeVsSignificanceLevel()

  # Fig 7A
  correctAntibioticsVsIndex(MODE_A)
  # Fig7B
  antibioticsVsMetric(MODE_A)

  #Fig 8A-D
  correctRAntibioticsVsSignificanceLevel()
  correctSAntibioticsVsSignificanceLevel()
  meRateAntibioticsVsSignificanceLevel()
  vmeRateAntibioticsVsSignificanceLevel()
  
  #Fig 7prim
#  meRateAntibioticsVsSignificanceLevel(prefix = "oneperabgroup")
#  vmeRateAntibioticsVsSignificanceLevel(prefix = "oneperabgroup")
  
  
# meRateAntibioticsVsSignificanceLevel(index=13)
# vmeRateAntibioticsVsSignificanceLevel(index=13)

  for(aPrefix in c(NA,ONEPERABGROUP)){
    
  meRateAntibioticsVsIndex(prefix = aPrefix)
  vmeRateAntibioticsVsIndex(prefix = aPrefix)
  }
  

  meRateAntibioticsVsIndex <- readStatisticsExcel(name ="meRateAntibioticsVsIndex-4_6_8-Mode-A")
  meRateAntibioticsVsIndexPerGroup <- readStatisticsExcel(name ="oneperabgroup-meRateAntibioticsVsIndex-4_6_8-Mode-A")
  meMergedRateAntibioticsVsIndex <- meRateAntibioticsVsIndex %>% dplyr::full_join(meRateAntibioticsVsIndexPerGroup,by="antibiotic",suffix = c("", "g"),) %>% dplyr::select(antibiotic,"4","4g","6","6g","8","8g")
  writeStatisticsExcel(meMergedRateAntibioticsVsIndex,name="meMergedRateAntibioticsVsIndex",shortName="meMergedRateAntibioticsVsIndex")

  vmeRateAntibioticsVsIndex <- readStatisticsExcel(name ="vmeRateAntibioticsVsIndex-4_6_8-Mode-A")
  vmeRateAntibioticsVsIndexPerGroup <- readStatisticsExcel(name ="oneperabgroup-vmeRateAntibioticsVsIndex-4_6_8-Mode-A")
  vmeMergedRateAntibioticsVsIndex <- vmeRateAntibioticsVsIndex %>% dplyr::full_join(vmeRateAntibioticsVsIndexPerGroup,by="antibiotic",suffix = c("", "g"),) %>% dplyr::select(antibiotic,"4","4g","6","6g","8","8g")
  writeStatisticsExcel(vmeMergedRateAntibioticsVsIndex,name="vmeMergedRateAntibioticsVsIndex",shortName="vmeMergedRateAntibioticsVsIndex")
  
  for(mode in MODES){
    indexRange = 4:13
    meRateAntibioticsVsIndex(indexRange=indexRange,mode=mode)
    vmeRateAntibioticsVsIndex(indexRange=indexRange,mode=mode)
    meRateAntibioticsVsIndex(prefix = "oneperabgroup",indexRange=indexRange,mode=mode)
    vmeRateAntibioticsVsIndex(prefix = "oneperabgroup",indexRange=indexRange,mode=mode)
  }
  
  
  
  #Fig S1
  meRateModeVsSignificanceLevel()
  vmeRateModeVsSignificanceLevel()
  
  # Fig S2
  correctSAntibioticsVsMode()
  meRateAntibioticsVsMode()
  correctRAntibioticsVsMode()
  vmeRateAntibioticsVsMode()
  
  # Heatmaps
  meRateSampleVsAntibiotic()
  vmeRateSampleVsAntibiotic()
  #meRateSampleVsAntibiotic(index = 13)
  #vmeRateSampleVsAntibiotic(index = 13)
  
  
  #Poster
  #meRateAntibioticsVsMode(prefix = "oneperabgroup")
  #vmeRateAntibioticsVsMode(prefix = "oneperabgroup")
}
