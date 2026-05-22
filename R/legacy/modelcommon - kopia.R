source(file = 'common - kopia.R')
library(stringr)


signis = c("0025", "005", "01")

# Uncomment to skip conformal predictions
#signis <- c()

METRICS_COLS <- c("correct","ME","VME")

BEST_ABS <- list (  
  c ("AMP","CRO","CIP","TOB") , 
  c("AMP","AMC","CRO","CIP","OFX","TOB"),  
  c("AMC","PIP","TZP","CAZ","CRO","CIP","OFX","TOB"))
BEST_ALT <- list(
  c("AMC","CTX","OFX","GEN"),
  c("PIP","TZP","CTX","LVX","MFX","GEN"),
  c("AMC" ,"TZP" ,"CAZ", "CTX", "FEP", "OFX","LVX" ,"GEN"),
  c("AMP", "PIP" ,"CRO","CTX","FEP", "CIP", "MFX" ,"TOB"))

SELECTED <- append(BEST_ABS,BEST_ALT)

BEST_SELECTION <- c(fourbest = "TZP_CAZ_LVX_TOB",fivebest = "AMC_TZP_CAZ_LVX_TOB",sixbest = "AMC_TZP_CAZ_LVX_MFX_TOB")
ONEPERABGROUP <- "oneperabgroup"
SUBGROUPS <- c(names(BEST_SELECTION),ONEPERABGROUP)

PENICILLINS <- c("AMP","AMC","PIP","TZP")
CEPHALOSPORINS <- c("CAZ","CRO","CTX","FEP")
FLOUROQUINOLONS <- c("CIP" ,"OFX" ,"LVX" ,"MFX")
AMINOGLYCOSIDES <- c("GEN" ,"TOB")
AB_GROUPS <- list(PENICILLINS,CEPHALOSPORINS,FLOUROQUINOLONS,AMINOGLYCOSIDES)



COMPARE_CORRECT <- "correct"
COMPARE_VME <- "VME"
COMPARE_ME <- "ME"
COMPARE_PREDICTED <- "predicted"
COMPARE_AMBIGUOUS <- "ambiguous"
COMPARE_PREDICTED_OR_AMBIGUOUS <- "predorambiguous"
COMPARE_NOT_PREDICTED <- "notpredicted"

#Names in the model, and basic order in which they should
ALL_ANTIBIOTICS_IN_MODEL <- c("AMP","AMC","PIP","TZP","CAZ","CRO","CTX","FEP","CIP","OFX","LVX","MFX","GEN","TOB")

paste_non_na <- function(..., sep = "_") {
  x <- c(...)
  x <- x[!is.na(x)]
  paste(x, collapse = sep)
}

METRICS <- c(
  "S", #answer S, pred any
  "R", #answer R, pred any
  "correctS", #answer S, pred S
  "correctR", #answer R, pred R
  "falseS", #answer R, pred S
  "falseR", #answer S, pred R
  "zerolabelS", #answer S, pred empty
  "zerolabelR", #answer R, pred empty
  "twolabelS",  #answer S, pred S/R
  "twolabelR", #answer R, pred S/R
  "total"      #answer any, pred any
)



NUMBER_OF_ANTIBIOTICS <- 14
ANTIBIOTICS = c("AMP",	"AMC"	,"PIP"	,"TZP",	"CAZ",	"CRO",	"CTX"	,"FEP"	,"CIP"	,"OFX"	,"LVX"	,"MFX"	,"GEN"	,"TOB")
PASCAL = c(14,91,364,1001,2002,3003,3432,3003,2002,1001,364,91,14)

range <- 4:(NUMBER_OF_ANTIBIOTICS - 1)


getModelInputFolder <- function()
{
  paste(modelDirectory,
        "input",
        sep = "/")
}


getCommonModelFolder <- function()
{
  paste(modelDirectory,
        "output",
        modelVersion,
        "common",
        sep = "/")
}


getPredictionBase <- function(mode){
	  PREDICTION_BASE <-
	  paste(modelDirectory,
	        "output",
	        modelVersion,
	        mode,
	        sep = "/")
	PREDICTION_BASE
}


getPredictionAltFolder <- function(mode)
{
  paste(getPredictionBase(mode),
        "predsalt",
        sep = "/")
  
}

getConformalPredictionFolder <- function(mode)
{
  paste(getPredictionBase(mode),
        "predscp",
        sep = "/")
}

getCommonModelFolder <- function()
{
  paste(modelDirectory,
        "output",
        modelVersion,
        "common",
        sep = "/")
}


getPredictionFolder <- function(mode)
{
	paste(getPredictionBase(mode),
        "preds",
        sep = "/")
		
}

getCompareFolder <-function(mode)
{
	paste(getPredictionBase(mode),"compare", sep = "\\")	
}


getStatisticsFolder <-function(mode)
{
	paste(getPredictionBase(mode),"statistics", sep = "\\")	
}

getStatisticsTmpFolder <-function(mode)
{
  paste(getStatisticsFolder(mode),"tmp", sep = "\\")
  
}

getStatisticsSubsetFolder <-function(mode)
{
  paste(getStatisticsFolder(mode),"subset", sep = "\\")
}



getStatisticsBestFolder <-function(mode)
{
	paste(getStatisticsFolder(mode),"best", sep = "\\")

}

getStatisticsSampleWordFolder <-function(mode)
{
	paste(getStatisticsFolder(mode), "words", sep = "\\")

}

getStatisticsMetricsFolder <- function(mode)
{
	paste(getStatisticsFolder(mode), "metrics", sep = "\\")
}

getStatisticsSampleMetricsFolder <- function(mode)
{
	paste(getStatisticsMetricsFolder(mode),"samples", sep = "\\")
}

getStatisticsAntibioticsMetricsFolder <- function(mode)
{
	paste(getStatisticsMetricsFolder(mode),"antibiotics", sep = "\\")
}


TEMP_FOLDER <- paste(modelDirectory, "temp", sep = "\\")




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
           folder)
  {
    shortfile <- makeShortFileName(type,k,signi)
    readOutputFrameShort(shortfile = shortfile,folder = folder)
  }


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

