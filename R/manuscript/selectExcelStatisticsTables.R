source(file = 'manuscript/manuscriptcommon.R')
source(file = 'model/modelcommon.R')




mergeRiskStratifiedStatistics <-function()
{
  names <- c("riskForMEStatistics",
             "riskForVMEStatistics",
             "riskForMEStatisticsAb",
             "riskForVMEStatisticsAb")
  for(name in names){
    orgName <- name
    framesAllModes <- lapply(MODES,function(mode){readStatisticsFrameCSV(name,"",mode)})
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
  names <- c("abErrorStatistics",
             "abErrorStatisticsNormalized",
             "errorStatistics",
             "errorStatisticsNormalized")
  
  for(name in names){
    orgName <- name
    if(!is.na(prefix)){
      name <- paste(prefix,name,sep="-")
    }
    #name <- "abErrorStatistics"
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
  #significanceRange <- c(NA,"0.100")
  if(!myIsNa(significanceRange)){
    subFrame <- subFrame %>% dplyr::filter(significanslevel %in% significanceRange)
  } else {
    subFrame <- subFrame %>% dplyr::filter(is.na(significanslevel))
  }
  if(length(unique(subFrame$significanslevel))==1){
    subFrame <- subFrame %>% dplyr::select(!significanslevel) 
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

filterMetricRange <- function(frame,metricRange){
  #frame <- readStatisticsExcel(name)
  longSubframe <- frame %>% tidyr::pivot_longer(cols = METRICS_COLS, names_to = "metric")
  #metricRange <- c("ME","VME") 
  #metricRange <- "correct"
  #metricRange <- "ME"
  
  longSubframe <- longSubframe %>% dplyr::filter(metric %in% metricRange)
  longSubframe %>% tidyr::pivot_wider(names_from = metric ,values_from = value)
}

filterMetricRangeNorm <- function(frame,metricRange){
  #frame <- readStatisticsExcel(name)
  longSubframe <- frame %>% tidyr::pivot_longer(cols = c("ME","VME"), names_to = "metric")
  #metricRange <- c("ME","VME") 
  #metricRange <- "correct"
  #metricRange <- "ME"
  
  longSubframe <- longSubframe %>% dplyr::filter(metric %in% metricRange)
  longSubframe %>% tidyr::pivot_wider(names_from = metric ,values_from = value)
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
  frame <- readStatisticsExcel("abErrorStatistics")
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
  frame <- readStatisticsExcel("abErrorStatistics")
  subFrame <- frame
  subFrame <- subFrame %>% filterMode(mode) %>% filterSignificanceRange(NA) %>% reorderMetrics()
  
  indexRange <- c(4,6,8)
  aList <- lapply(indexRange, function(x) { subFrame %>% filterIndexRange(x)})
  names(aList) <-  as.character(indexRange)
  writeStatisticsExcel(aList,paste("antibioticsVsMetric",mode,sep="-"),paste("abVsMetric",substring(mode,nchar(mode),nchar(mode)),sep="-"))
}

correctAntibioticsVsIndex <- function(mode  = "Mode-A")
{
  frame <- readStatisticsExcel("abErrorStatistics")
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
metricsVsIndex <- function(mode  = "Mode-A")
{
  frame <- readStatisticsExcel("errorStatistics")
  subFrame <- frame
  subFrame <- subFrame %>% filterMode(mode) %>% filterSignificanceRange(NA) %>% filterIndexRange(c(4:8))
  longSubframe <- subFrame %>% tidyr::pivot_longer(cols = METRICS_COLS, names_to = "metric")
  subFrame <- longSubframe %>% tidyr::pivot_wider(names_from = noinputab ,values_from = value)

  
  writeStatisticsExcel(subFrame,paste("metricsVsIndex",mode,sep="-"),paste("abVsMetric",substring(mode,nchar(mode),nchar(mode)),sep="-"))
}

metricsVsModeForSixAb <- function(significanceLevel = NA)
{
  frame <- readStatisticsExcel("errorStatistics")
  subFrame <- frame
  subFrame <- subFrame %>% filterSignificanceRange(significanceLevel) %>% filterIndexRange(6)
  longSubframe <- subFrame %>% tidyr::pivot_longer(cols = METRICS_COLS, names_to = "metric")
  subFrame <- longSubframe %>% tidyr::pivot_wider(names_from = mode ,values_from = value)
  
  
  writeStatisticsExcel(subFrame,paste("metricsVsMode",6,significanceLevel,sep="-"),paste("metricsVsMode",6,significanceLevel,sep="-"))
                       
}
metricsVsModeForSignificaneLevelNA <- function(index = 6)
{
  frame <- readStatisticsExcel("errorStatistics")
  subFrame <- frame
  subFrame <- subFrame %>% filterSignificanceRange(NA) %>% filterIndexRange(index)
  longSubframe <- subFrame %>% tidyr::pivot_longer(cols = METRICS_COLS, names_to = "metric")
  subFrame <- longSubframe %>% tidyr::pivot_wider(names_from = mode ,values_from = value)
  
  
  writeStatisticsExcel(subFrame,paste("metricsVsMode",index,sep="-"),paste("metricsVsMode",index,sep="-"))
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
  frame <- readStatisticsExcel("errorStatistics")
  subFrame <- frame
  subFrame <- subFrame %>% 
    filterMode(mode) %>% 
    filterIndexRange(c(4:8)) %>% 
    filterMetricRange("correct") %>% 
    filterSignificanceRange(c(NA,"0.100","0.050","0.025"))
  subFrame <- subFrame %>% tidyr::pivot_wider(names_from = noinputab ,values_from = correct)
  
  
  
  subFrame <- subFrame %>% dplyr::slice(4:1)
  subFrame$significanslevel <- unlist(lapply(subFrame$significanslevel,function(x) {signiToPercent(x)}))
  writeStatisticsExcel(subFrame,paste("correctSignificanceLevelVsIndex",mode,sep="-"),paste("correctSignilVsIndex",substring(mode,nchar(mode),nchar(mode)),sep="-"))
}

correctModeVsSignificanceLevel <- function(index  = 6)
{
  frame <- readStatisticsExcel("errorStatistics")
  subFrame <- frame
  subFrame <- subFrame %>% 
    filterIndexRange(index) %>% 
    filterMetricRange("correct") %>% 
    filterSignificanceRange(c(NA,"0.100","0.050","0.025"))
  subFrame$significanslevel <- unlist(lapply(subFrame$significanslevel,function(x) {signiToPercent(x)}))
  
  subFrame <- subFrame %>% tidyr::pivot_wider(names_from = significanslevel ,values_from = correct)

  subFrame <- subFrame %>% dplyr::select(mode ,`10%`,`5%`,`2.5%`)  

  writeStatisticsExcel(subFrame,paste("correctModeVsSignificanceLevel",index,sep="-"),paste("correctModeVsSigni",index,sep="-"))
}


meNormalizedAntibioticsVsSignificanveLevel <- function(index  = 6,prefix=NA,subfolder = "")
{
  name <-  "abErrorStatisticsNormalized"
  if(!is.na(prefix)){
    name <- paste(prefix,name,sep="-")
  }
  frame <- readStatisticsExcel(name,subfolder=subfolder)
  #frame <- readStatisticsFrameCSV(name="nonATUstats-abErrorStatisticsNormalized",mode="Mode-A")
  subFrame <- frame
  subFrame <- subFrame %>% 
    filterMode(MODE_A) %>% 
    filterIndexRange(index) %>% 
    filterMetricRangeNorm("ME") %>% 
    filterSignificanceRange(c(NA,"0.100","0.050","0.025"))
  subFrame$significanslevel <- unlist(lapply(subFrame$significanslevel,function(x) {signiToPercent(x)}))
  
  subFrame <- subFrame %>% tidyr::pivot_wider(names_from = significanslevel ,values_from = ME)
  subFrame <- subFrame %>% dplyr::select(antibiotic,`non-conformal` ,`10%`,`5%`,`2.5%`)  
  
  #writeStatisticsExcel(subFrame,paste("nonATUstats-meNormalizedAntibioticsVsSignificanveLevel",index,sep="-"),paste("meNormAntibioticsVsSigni",index,sep="-"))
  outname <- "meNormalizedAntibioticsVsSignificanveLevel"
  if(!is.na(prefix)){
    outname <- paste(prefix,outname,index,sep="-")
  } else {
    outname <- paste(outname,index,sep="-")
  }
  
  writeStatisticsExcel(subFrame,outname,paste("meNormAntibioticsVsSigni",index,sep="-"),subfolder=subfolder)
}

vmeNormalizedAntibioticsVsSignificanveLevel <- function(index  = 6,prefix=NA,subfolder = "")
{
  name <-  "abErrorStatisticsNormalized"
  if(!is.na(prefix)){
    name <- paste(prefix,name,sep="-")
  }
  frame <- readStatisticsExcel(name,subfolder = subfolder)
  subFrame <- frame
  subFrame <- subFrame %>% 
    filterMode(MODE_A) %>% 
    filterIndexRange(index) %>% 
    filterMetricRangeNorm("VME") %>% 
    filterSignificanceRange(c(NA,"0.100","0.050","0.025"))
  subFrame$significanslevel <- unlist(lapply(subFrame$significanslevel,function(x) {signiToPercent(x)}))
  
  subFrame <- subFrame %>% tidyr::pivot_wider(names_from = significanslevel ,values_from = VME)
  subFrame <- subFrame %>% dplyr::select(antibiotic,`non-conformal` ,`10%`,`5%`,`2.5%`)  
  
  outname <- "vmeNormalizedAntibioticsVsSignificanveLevel"
  if(!is.na(prefix)){
    outname <- paste(prefix,outname,index,sep="-")
  } else {
    outname <- paste(outname,index,sep="-")
  }
  
  writeStatisticsExcel(subFrame,outname,paste("vmeNormAntibioticsVsSigni",index,sep="-"),subfolder = subfolder)
}


meMetricVsAntibiotic <- function(index = 6,significanceLevel = NA)
{
  name <-  "riskForMEStatisticsAb"
  
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

vmeMetricVsAntibiotic <- function(index = 6,significanceLevel = NA)
{
  name <-  "riskForVMEStatisticsAb"
  
  frame <- readStatisticsExcel(name)
  subFrame <- frame
  subFrame <- subFrame %>% 
    filterIndexRange(index) %>% 
    filterSignificanceRange(significanceLevel) %>%
    filterMode(modeRange = MODE_A)
  #Swap
  subFrame <- subFrame %>% tidyr::pivot_longer(cols = c("correct","VME"),names_to = "metric")  %>%
    tidyr::pivot_wider(names_from = antibiotic, values_from = value)
  name <- paste("vmeMetricVsAntibiotic",index,significanceLevel,sep="-")
  writeStatisticsExcel(subFrame,name,name)
}



meNormalizedModeVsSignificanceLevel <- function(index  = 6,prefix=NA,subfolder = "")
{
  name <-  "errorStatisticsNormalized"
  if(!is.na(prefix)){
    name <- paste(prefix,name,sep="-")
  }
  frame <- readStatisticsExcel(name,subfolder = subfolder)
  subFrame <- frame
  subFrame <- subFrame %>% 
    filterIndexRange(index) %>% 
    filterMetricRangeNorm("ME") %>% 
    filterSignificanceRange(c(NA,"0.100","0.050","0.025"))
  subFrame$significanslevel <- unlist(lapply(subFrame$significanslevel,function(x) {signiToPercent(x)}))
  
  subFrame <- subFrame %>% tidyr::pivot_wider(names_from = significanslevel ,values_from = ME)
  
  subFrame <- subFrame %>% dplyr::select(mode ,`10%`,`5%`,`2.5%`)  
  
  outname <- "meNormalizedModeVsSignificanceLevel"
  if(!is.na(prefix)){
    outname <- paste(prefix,outname,index,sep="-")
  } else {
    outname <- paste(outname,index,sep="-")
  }
  writeStatisticsExcel(subFrame,outname,paste("meNormalizedModeVsSigni",index,sep="-"),subfolder = subfolder)
}

normalizedErrorVsSignificanceLevelNonATU <- function(index  = 6)
{
  subfolder <- "subset"
  frame <- readStatisticsFrameCSV(name="nonATUstats-errorStatisticsNormalized",mode="Mode-A",subfolder = subfolder)
  subFrame <- frame
  subFrame <- subFrame %>% 
    filterIndexRange(index) %>% 
#    filterMetricRangeNorm("VME") %>% 
    #dplyr::filter(metric=="VME") %>% dplyr::select(-VME) %>% 
    filterSignificanceRange(c(NA,"0.100","0.050","0.025"))
  subFrame$significanslevel <- unlist(lapply(subFrame$significanslevel,function(x) {signiToPercent(x)}))
  
  subFrame <- subFrame %>% tidyr::pivot_longer(cols = c(ME,VME),names_to = "metric")
  subFrame <- subFrame %>% tidyr::pivot_wider(names_from = significanslevel ,values_from = value)
  
  subFrame <- subFrame %>% dplyr::select(metric ,`non-conformal`,`10%`,`5%`,`2.5%`)  
  
  writeStatisticsExcel(subFrame,paste("normalizedErrorVsSignificanceLevelNonATU",index,sep="-"),paste("normalizedErrorVsSigniNonATU",index,sep="-"),subfolder=subfolder)
}


vmeNormalizedModeVsSignificanceLevel <- function(index  = 6,prefix=NA,subfolder="")
{
  name <-  "errorStatisticsNormalized"
  if(!is.na(prefix)){
    name <- paste(prefix,name,sep="-")
  }
  frame <- readStatisticsExcel(name,subfolder=subfolder)
  subFrame <- frame
  subFrame <- subFrame %>% 
    filterIndexRange(index) %>% 
    filterMetricRangeNorm("VME") %>% 
    filterSignificanceRange(c(NA,"0.100","0.050","0.025"))
  subFrame$significanslevel <- unlist(lapply(subFrame$significanslevel,function(x) {signiToPercent(x)}))
  
  subFrame <- subFrame %>% tidyr::pivot_wider(names_from = significanslevel ,values_from = VME)
  
  subFrame <- subFrame %>% dplyr::select(mode ,`10%`,`5%`,`2.5%`)  
  
  outname <- "vmeNormalizedModeVsSignificanceLevel"
  if(!is.na(prefix)){
    outname <- paste(prefix,outname,sep="-")
  } else {
    outname <- paste(outname,index,sep="-")
  }
  writeStatisticsExcel(subFrame,outname,paste("vmeNormalizedModeVsSigni",index,sep="-"),subfolder=subfolder)
}



correctModeVsIndex <- function(significanceLevel  = "0.100")
{
  frame <- readStatisticsExcel("errorStatistics")
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
  frame <- readStatisticsExcel("errorStatistics")
  subFrame <- frame
  subFrame <- subFrame %>% 
    filterMode(mode) %>% 
    filterIndexRange(indexRange) %>% 
    filterSignificanceRange(c(NA,"0.100","0.050","0.025"))
  subFrame <- subFrame %>% tidyr::pivot_longer(cols = METRICS_COLS,names_to = "metric" )
  subFrame$significanslevel <- unlist(lapply(subFrame$significanslevel,function(x) {signiToPercent(x)}))
  subFrame <- subFrame %>% tidyr::pivot_wider(names_from = significanslevel ,values_from = value)
  subFrame <- subFrame %>% dplyr::select(noinputab, metric ,`non-conformal`,`10%`,`5%`,`2.5%`)  
  aList <- lapply(indexRange, function(x) { subFrame %>% filterIndexRange(x)})
  names(aList) <-  as.character(indexRange)
  
  writeStatisticsExcel(aList,paste("metricVsSignificanceLevel",mode,sep="-"),paste("metricVsSigni",substring(mode,nchar(mode),nchar(mode)),sep="-"))
}



ALL <- function ()
{
  #merge data from different modes
  mergeModesErrorStatistics()
  lapply(SELECTED,function(x){mergeModesErrorStatistics(prefix = paste(x,collapse="_"),subfolder = "subset")})
  mergeModesErrorStatistics()
  
  mergeModesErrorStatistics(prefix = "nonATUstats",subfolder = "subset")
  
  #general level
  correctSignificanceLevelVsIndex(MODE_A)
  metricsVsIndex(MODE_A)
  metricVsSignificanceLevel(MODE_A)
  correctModeVsSignificanceLevel(index = 6)
  correctModeVsIndex(significanceLevel = "0.100")

  metricsVsModeForSixAb("0.100")
  metricsVsModeForSixAb("0.050")
  metricsVsModeForSixAb("0.025")
  metricsVsModeForSignificaneLevelNA(4)
  metricsVsModeForSignificaneLevelNA(6)
  metricsVsModeForSignificaneLevelNA(8)
  
  meNormalizedModeVsSignificanceLevel()
  vmeNormalizedModeVsSignificanceLevel()
  lapply(SELECTED,function(x){meNormalizedModeVsSignificanceLevel(length(x),paste(x,collapse="_"),subfolder="subset")})
  lapply(SELECTED,function(x){vmeNormalizedModeVsSignificanceLevel(length(x),paste(x,collapse="_"),subfolder="subset")})
  
  #antibiotics level  
  antibioticsVsMetric(MODE_A)
  antibioticsVsMode(6)
  correctAntibioticsVsIndex(MODE_A)
  
  meNormalizedAntibioticsVsSignificanveLevel()
  vmeNormalizedAntibioticsVsSignificanveLevel()
  lapply(SELECTED,function(x){meNormalizedAntibioticsVsSignificanveLevel(length(x),paste(x,collapse="_"),subfolder="subset")})
  lapply(SELECTED,function(x){vmeNormalizedAntibioticsVsSignificanveLevel(length(x),paste(x,collapse="_"),subfolder="subset")})
  
  meNormalizedAntibioticsVsSignificanveLevel(prefix = "nonATUstats",subfolder="subset")
  vmeNormalizedAntibioticsVsSignificanveLevel(prefix = "nonATUstats",subfolder="subset")
  normalizedErrorVsSignificanceLevelNonATU()
  
  
  
  meMetricVsAntibiotic()
  vmeMetricVsAntibiotic()
  
}


