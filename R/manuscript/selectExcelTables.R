###########

source(file = 'common.R')

library(openxlsx)
#library(writexl)

#install.packages('xlsx')     
#library(xlsx)   


MODE_A <- "ModeA"
MODE_B <- "ModeB"
MODE_C <- "ModeC"

MODES <- c(MODE_A,MODE_B,MODE_C)

 

getModeRootStatisticsFolder <- function(mode)
{
  paste(modelDirectory,
        "output",
        mode,
        "statistics",
        sep = "/")
}


readStatisticsFrameCSV  <-
  function(name,
           subfolder="",mode)
  {
    file <- paste(getModeRootStatisticsFolder(mode),subfolder,paste(name,"csv",sep="."),  sep = "/")
    
    read.csv2(file = file,check.names = FALSE)
  } 


writeStatisticsExcel <- function(frame,name,mode="")
  {
  file <- paste(manuscriptDirectory,paste(mode,name,".xlsx",sep=""),  sep = "/")
#  print(file)
#  print(frame)
  write.xlsx( x= frame, file=file) 
#  print("ok")
  }



writeAllModes <-function(name,subfolder="")
{
  aList <- list()
  for(mode in MODES){
    frame <- readStatisticsFrameCSV(name,subfolder,mode)
    aList[[mode]] <- frame
  }

  writeStatisticsExcel(aList,name)
}

writeSeveral <- function(names,commonname,subfolder,mode)
{
  aList <- list()
  for(name in names){
    frame <- readStatisticsFrameCSV(name,subfolder,mode)
    aList[[name]] <- frame
  }
  
  writeStatisticsExcel(aList,commonname,mode)
}


ALL <- function()
{
  writeAllModes("allStats")
  writeSeveral(c("abStatistics-4","abStatistics-6","abStatistics-8"),"VME-ME-CORRECT","",mode=MODE_A)
 
  
  
  clinicalSigniStatistics()
}



clinicalSigniStatistics <- function()
{
  aList <- list()
  name <- "clinicalSigniStatistics"
  for(mode in MODES){
    frame <- readStatisticsFrameCSV(name,"",mode)
    frame <- frame[frame$significanceLevel %in% c("crude","0.025","0.050","0.100"),colnames(frame) %in% c("significanceLevel",as.character(4:8))]                                                                                           
    frame[frame=="crude"] <- "nonconformal"
    aList[[mode]] <- frame
  }
  writeStatisticsExcel(aList,name)
}

technicalSigniStatistics <- function()
{
  aList <- list()
  name <- "technicalSigniStatistics"
  for(mode in MODES){
    frame <- readStatisticsFrameCSV(name,"",mode)
    frame <- frame[frame$significanceLevel %in% c("0.025","0.050","0.100"),colnames(frame) %in% c("significanceLevel",as.character(4:13))]                                                                                           
    aList[[mode]] <- frame
  }
  writeStatisticsExcel(aList,name)
}


