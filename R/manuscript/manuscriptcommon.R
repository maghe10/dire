###########

#source(file = 'common.R')
source(file = 'model/modelcommon.R')
library(scales)

# FIXME migrate to readxl
library(openxlsx)

library(writexl)

MODE_A <- "Mode-A"
MODE_B <- "Mode-B"
MODE_C <- "Mode-C"

MODES <- c(MODE_A,MODE_B,MODE_C)

METRICS <- c("correct","ME","VME")
METRICS_LOWER_CASE <- c("correct","me","vme")

format_percent <- function(x) {
  x <- as.double(x)
  formatted <- format(round(x * 100, 10), trim = TRUE, scientific = FALSE)
  paste0(formatted, "%")
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

