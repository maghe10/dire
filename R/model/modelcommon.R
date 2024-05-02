source(file = 'common.R')
library(stringr)


MODE <- "ModeA"
minusAMC <- FALSE

if(!minusAMC)
  {
  PREDICTION_BASE <-
  paste(modelDirectory,
        "output",
        MODE,
        sep = "/")
  NUMBER_OF_ANTIBIOTICS <- 14
  ANTIBIOTICS = c("AMP",	"AMC"	,"PIP"	,"TZP",	"CAZ",	"CRO",	"CTX"	,"FEP"	,"CIP"	,"OFX"	,"LVX"	,"MFX"	,"GEN"	,"TOB")
  PASCAL = c(14,91,364,1001,2002,3003,3432,3003,2002,1001,364,91,14)
} else {
  PREDICTION_BASE <-
    paste(modelDirectory,
          "output",
          MODE,
          "minusAMC",
          sep = "/")
  NUMBER_OF_ANTIBIOTICS <- 13
  ANTIBIOTICS = c("AMP",	"PIP"	,"TZP",	"CAZ",	"CRO",	"CTX"	,"FEP"	,"CIP"	,"OFX"	,"LVX"	,"MFX"	,"GEN"	,"TOB")
  PASCAL = c(13,78,286,715,1287,1716,1716,1287,715,286,78,13)
}

#PREDICTION_BASE <-
#  paste(modelDirectory,
#        "finished",
#        "I_to_R_noDate",
#        "SR_as_word",
#        sep = "/")



PREDICTION_FOLDER <-
  paste(PREDICTION_BASE,
        "preds",
        sep = "/")
COMPARE_FOLDER <- paste(PREDICTION_BASE,"compare", sep = "\\")

STATISTICS_FOLDER <- paste(PREDICTION_BASE,"statistics", sep = "\\")

STATISTICS_BEST_FOLDER <- paste(STATISTICS_FOLDER,"best", sep = "\\")
STATISTICS_SAMPLE_WORDS_FOLDER <- paste(STATISTICS_FOLDER,"words", sep = "\\")

STATISTICS_METRICS_FOLDER <- paste(STATISTICS_FOLDER,"metrics", sep = "\\")
STATISTICS_SAMPLE_METRICS_FOLDER <- paste(STATISTICS_METRICS_FOLDER,"samples", sep = "\\")
STATISTICS_ANTIBIOTICS_METRICS_FOLDER <- paste(STATISTICS_METRICS_FOLDER,"antibiotics", sep = "\\")

TEMP_FOLDER <- paste(modelDirectory, "temp", sep = "\\")

signis = c("001", "0025", "005", "01")
range <- 1:(NUMBER_OF_ANTIBIOTICS - 1)



maxPredictions <- function(frame, index)
{
  nrow(frame) * ncol(frame) * (NUMBER_OF_ANTIBIOTICS - index)
}


applyElementWise <- function(x,y,aFunction)
{
  resultMatrix <- mapply(function(u,v){mapply(aFunction,u,v)},x,y)
  resultFrame <- as.data.frame(resultMatrix)
  rownames(resultFrame) <- rownames(x)
  resultFrame
}  


countWords <- function(frame,word)
{
  count <- 0
  aVector <- lapply(frame, function(x) { sum(ifelse(str_detect(x, word), 1, 0))})
  sum(unlist(aVector))
}


fixSigni <- function(signi)
{
  if(is.null(signi) || is.na(signi)){
    rv <- "nonconformal"
  } else {
    sTmp <- paste("0.", signi, sep = "")
    dTmp <- as.double(sTmp) * 10
    rv <- format(round(dTmp, 3), nsmall = 3)
  }
  rv
}


readInputframe <- function(k)
{
  file <-
    paste(modelDirectory,
          paste("modelInput", k, ".csv", sep = ""),
          sep = "/")
  # sample name is row name
  frame <- read.csv2(file = file, row.names = "sample")
  frame
}

readOutputFrameShort <-
  function(shortfile,
           folder)
  {
    file <- paste(folder, shortfile,  sep = "/")
    frame <- read.csv2(file = file, row.names = "sample")
    cn <- colnames(frame)
    index  <- grep("^X$", colnames(frame))
    returnFrame <- frame
    #There might be a first index column
    if (length(index) > 0) {
      returnFrame <- frame[, -index]
    }
    returnFrame
  }

makeShortFileName <-function(type,index,signi=NULL)
{
  if(is.null(signi)){
    paste(paste("modelOutput",type,as.character(index),sep="_"),"csv",sep=".")
  }
  else {
    paste(paste("modelOutput",type,as.character(index),signi,sep="_"),"csv",sep=".")
    
  }
}


readOutputFrame <-
  function(k,
           signi = NULL ,
           type,
           folder = PREDICTION_FOLDER)
  {
    shortfile <- makeShortFileName(type,k,signi)
    readOutputFrameShort(shortfile = shortfile,folder = folder)
  }

# readOutputFrame <-
#   function(k,
#            signi = NULL ,
#            type,
#            folder = PREDICTION_FOLDER)
#   {
#     if (is.null(signi)) {
#       shortfile <-   paste("modelOutput_", type, "_", k, ".csv", sep = "")
#     }
#     else {
#       shortfile <-
#         paste("modelOutput_", type, "_", k, "_", signi, ".csv", sep = "")
#     }
#     file <- paste(folder, shortfile,  sep = "/")
#     frame <- read.csv2(file = file, row.names = "sample")
#     cn <- colnames(frame)
#     index  <- grep("^X$", colnames(frame))
#     returnFrame <- frame
#     #There might be a first index column
#     if (length(index) > 0) {
#       returnFrame <- frame[, -index]
#     }
#     returnFrame
#   }

writeOutputFrameShort  <-
  function(frame,
           shortfile,
           folder)
  {
    file <- paste(folder, shortfile,  sep = "/")
    
    write.csv2(sampleAsColumns(frame),
               file = file,
               row.names = FALSE)
  }    



writeOutputFrame  <-
  function(frame,
           k,
           signi = NULL ,
           type ,
           folder)
  {
    shortfile <- makeShortFileName(type,k,signi)
    writeOutputFrameShort(frame=frame,shortfile = shortfile,folder = folder)
}

# writeOutputFrame  <-
#   function(frame,
#            k,
#            signi = NULL ,
#            type ,
#            folder = paste(modelDirectory, "temp", sep = "/"))
#   {
#     if (is.null(signi)) {
#       shortfile <-   paste("modelOutput_", type, "_", k, ".csv", sep = "")
#     }
#     else {
#       shortfile <-
#         paste("modelOutput_", type, "_", k, "_", signi, ".csv", sep = "")
#     }
#     file <- paste(folder, shortfile,  sep = "/")
#     
#     write.csv2(sampleAsColumns(frame),
#                file = file,
#                row.names = FALSE)
#   }

writeStatisticsFrame  <-
  function(frame,
           name,
           folder,row.names = FALSE)
  {
    file <- paste(folder,paste(name,".csv",sep=""),  sep = "/")
    if(row.names){
      frame <- sampleAsColumns(frame)
    }
    write.csv2(frame,
               file = file,row.names = FALSE)
  }    

readStatisticsFrame  <-
  function(name,
           folder)
  {
    file <- paste(folder, paste(name,".csv",sep=""),  sep = "/")
    
    read.csv2(file = file)
  }    

