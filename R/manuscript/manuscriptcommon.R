###########

#source(file = 'common.R')
source(file = 'model/modelcommon.R')
source(file = 'manuscript/metric_helpers_common.R')

library(scales)

# FIXME migrate to readxl
library(openxlsx)

library(writexl)

MODE_A <- "Mode-A"
MODE_B <- "Mode-B"
MODE_C <- "Mode-C"

MODES <- c(MODE_A,MODE_B,MODE_C)

# call something else
#METRICS <- c("correct","ME","VME")
#METRICS_LOWER_CASE <- c("correct","me","vme")

DEFAULT_MODE <- "Mode-A"
DEFAULT_CPMODE <- "normal"

COLNAMES_LONG_ORDER <- c("cpmode","mode","word","sample","noinputab","antibiotic","significanceLevel", "metric","count")
COLNAMES_WIDE_ORDER <- c("cpmode","mode","word","sample","noinputab","antibiotic","significanceLevel", METRICS)

PENICILLINS <- "Penicillins"
CEPHALOSPORINS <- "Cephalosporins"
FLOUROQUINOLONS <- "Fluoroquinolones"
AMINOGLYCOSIDES <- "Aminoglycosides"

AB_GROUPS <- c(PENICILLINS,CEPHALOSPORINS,FLOUROQUINOLONS,AMINOGLYCOSIDES)


sortAb <- function(df)
{
  df %>% arrange(factor(antibiotic, levels = ANTIBIOTICS))
}

sortMetric <- function(df)
{
  df %>% arrange(factor(metric, levels = METRICS))
}

sortLong <-function(df)
{
  df %>% sortAb() %>% sortMetric()
}

sortWide <-function(df)
{
  df %>% sortAb()
}


normalizeLong <- function(df)
{
  df %>% select(any_of(COLNAMES_LONG_ORDER)) %>% sortLong()
}

normalizeWide <- function(df)
{
  df %>% select(any_of(COLNAMES_WIDE_ORDER)) %>% sortWide()
}


readStatisticsCommon <- function(name)
{
  readRDS(file.path(getCommonModelFolder(),sprintf("%s.rds",name))) %>% as_tibble()
}



readCountFrameLongGroups <- function()
{
  ab_groups <- tribble(
    ~antibiotic, ~ab_group,
    "AMP", PENICILLINS,
    "PIP", PENICILLINS,
    "AMC", PENICILLINS,
    "TZP", PENICILLINS,
    
    "CAZ", CEPHALOSPORINS,
    "CRO", CEPHALOSPORINS,
    "CTX", CEPHALOSPORINS,
    "FEP", CEPHALOSPORINS,
    
    "CIP", FLOUROQUINOLONS,
    "OFX", FLOUROQUINOLONS,
    "LVX", FLOUROQUINOLONS,
    "MFX", FLOUROQUINOLONS,
    
    "GEN", AMINOGLYCOSIDES,
    "TOB", AMINOGLYCOSIDES
  )
  
  readCountFrameLong() %>%
    left_join(ab_groups, by = "antibiotic")
  
}

readCountFrameWideSumByGroup <- function() {
  readCountFrameLongSumByGroup() %>% pivot_wider(names_from = metric, values_from = count)
}

readCountFrameLongSumByGroup <- function() {
  readCountFrameLong() %>%
    dplyr::mutate(
      ab_group = dplyr::case_when(
        antibiotic %in% c("AMP", "PIP", "AMC", "TZP") ~ "Penicillins",
        antibiotic %in% c("CAZ", "CRO", "CTX", "FEP") ~ "Cephalosporins",
        antibiotic %in% c("CIP", "OFX", "LVX", "MFX") ~ "Fluoroquinolones",
        antibiotic %in% c("GEN", "TOB")               ~ "Aminoglycosides",
        TRUE                                          ~ NA_character_
      )
    ) %>%
    dplyr::group_by(
      cpmode, mode, noinputab, ab_group, metric, significanceLevel
    ) %>%
    dplyr::summarise(
      count = sum(count, na.rm = TRUE),
      .groups = "drop"
    )
}



readCountFrameLongSum <- function(prefix=NA)
{
  readCountFrameLong(prefix) %>%
    group_by(cpmode, mode, noinputab, metric, significanceLevel) %>%
    summarise(
      count = sum(count, na.rm = TRUE),
      .groups = "drop"
    )
}


readCountSampleFrameLongSum_AggregationOnSample <-function()
  {
  readCountSampleFrameLong() %>%
    group_by(cpmode, mode, noinputab, antibiotic, metric, significanceLevel) %>%
    summarise(
      count = sum(count, na.rm = TRUE),
      .groups = "drop"
    )
}


readCountSampleFrameLongSum_AggregatedOnAntibiotic <-function(){
  readCountSampleFrameLong() %>%
    group_by(cpmode, mode, sample,noinputab, metric, significanceLevel) %>%
    summarise(
      count = sum(count, na.rm = TRUE),
      .groups = "drop"
    )
}



readCountFrameWideSum <- function(prefix=NA)
{
  readCountFrameLongSum(prefix) %>% pivot_wider(names_from = metric, values_from = count)
}


readCountFrameWide <- function(prefix=NA)
{
  dfCountWide <- readStatisticsCommon(paste_non_na(prefix,"abErrorStatisticsCount")) %>% normalizeWide()
  dfCountWide
}




readCountFrameLong <- function(prefix=NA)
{
  dfCountLong <- readStatisticsCommon(paste_non_na(prefix,"allStats_sum")) %>% normalizeLong()
  dfCountLong
}

readCountSampleFrameLong <- function(prefix=NA)
{
  dfCountLong <- readStatisticsCommon(paste_non_na(prefix,"allStats_sample")) %>% normalizeLong()
  dfCountLong
}

readCountWordFrameLong <- function(prefix=NA)
{
  dfCountLong <- readStatisticsCommon(paste_non_na(prefix,"allStats_word")) %>% normalizeLong()
  dfCountLong
}

readCountWordFrameLong_AggregateOnAntibiotic <-function()
{
  readCountWordFrameLong() %>%
    group_by(cpmode, mode, word, noinputab, metric, significanceLevel) %>%
    summarise(
      count = sum(count, na.rm = TRUE),
      .groups = "drop"
    )
}



filterNormal <- function(df){
  filter_and_drop(df,cpmode,"normal")
}

filterStrict <- function(df){
  filter_and_drop(df,cpmode,"strict")
}

defaultCpMode <- function(df)
{
  filter_and_drop(df,cpmode,DEFAULT_CPMODE)
}


defaultMode <- function(df){
  filter_and_drop(df,mode,DEFAULT_MODE)
}

defaults <-function(df)
{
  df %>% defaultCpMode() %>% defaultMode()
}

makeMetricLong <- function(df)
{
  df %>% pivot_longer(cols = all_of(METRICS),names_to = "metric",values_to = "count") %>% normalizeLong()
}

makeMetricWide <- function(df)
{
  df %>% pivot_wider(names_from = "metric",values_from = "count") %>% normalizeWide()
}




format_percent <- function(x) {
  x <- as.double(x)
  formatted <- format(round(x * 100, 10), trim = TRUE, scientific = FALSE)
  paste0(formatted, "%")
}


readModelMillimeterTable <- function()
{
  fileName <- "modelzonemillimeters.csv"
  millimeterTable <- read.csv2(file.path(getCommonModelFolder(),fileName))
  head(millimeterTable)
}

readCommonStatisticsCSV  <-
  function(name)
  {
    file <- paste(getCommonModelFolder(),paste(name,"csv",sep="."),  sep = "/")
    read.csv2(file = file,check.names = FALSE)
  }



readStatisticsFrameCSV  <-
  function(name,
           subfolder="",mode)
  {
    file <- paste(getStatisticsFolder(mode),subfolder,paste(name,"csv",sep="."),  sep = "/")
    
    read.csv2(file = file,check.names = FALSE)
  } 

readStatisticsExcel <- function(name,subfolder="")
{
  file <- paste(manuscriptDirectory,subfolder,paste(name,"xlsx",sep="."), sep = "/")
  readxl::read_xlsx(path=file)
}


writeStatisticsExcel <- function(frames,name,shortName,mode="",subfolder="")
{
  if(mode==""){
    shortFile <- paste(name,".xlsx",sep="") 
  } else {
    shortFile <- paste(name,"-",mode,".xlsx",sep="") 
  }
  
  file <- paste(manuscriptDirectory,subfolder,shortFile,  sep = "/")
  
  if(is.data.frame(frames)){
    aList <- list()
    if(mode==""){
      sheetName <- paste(shortName,sep="")
    } else {
      sheetName <- paste(shortName,"-",shortMode(mode),sep="")
    }
    aList[[sheetName]] <- frames
    frames <- aList
  } else { # list of frames
    if(mode==""){
      names(frames) <- paste(shortName,names(frames),sep="-")
    }
    else {
      names(frames) <- paste(shortName,names(frames),shortMode(mode),sep="-")
    }
  }
  #  print(file)
  #  print(frame)
  write_xlsx(x = frames, path=file)
  #write.xlsx( x= frames, file=file) 
  #  print("ok")
}

