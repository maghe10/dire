source(file='common.R')
library(stringr)

MODE <- "ModeA"

IN_FOLDER <- paste(modelDirectory,
      "output",
      MODE,
      "preds",
      sep = "/")
OUT_FOLDER  <- paste(modelDirectory,
           "output",
            MODE,
           "tmp",
           "preds",
           sep = "/")


asShortFileName <- function(k,signi = NULL ,type)
{
  if(is.null(signi)){
    shortfile <-   paste("modelOutput_",type,"_",k,".csv",sep="")
  }
  else {
    shortfile <-   paste("modelOutput_",type,"_",k,"_",signi,".csv",sep="")
  }
  shortfile
}

readOutputFrame <-function(shortfile,folder = IN_FOLDER )
{
  file <- paste(folder, shortfile,  sep="/")
  frame <- read.csv2(file=file,row.names = "sample")
  cn <- colnames(frame)
  index  <- grep("^X$", colnames(frame))
  returnFrame <- frame 
  #There might be a first index column
  if(length(index)>0){
    returnFrame <- frame[,-index]    
  }
  returnFrame
}


writeOutputFrame  <- function(frame,shortfile, folder = OUT_FOLDER)
{
  file <- paste(folder, shortfile,  sep="/")
  write.csv2(sampleAsColumns(frame),file=file,row.names = FALSE)
}


sieveColumns <- function(frame,antibiotic)
{
    returnFrame <- frame
    cols <- grep(antibiotic,colnames(returnFrame),invert=TRUE)
    if(length(cols)<2){
      returnFrame <- NULL
    } else {
      returnFrame <- returnFrame[,cols]
    }
    returnFrame
}

sieveValues <- function(frame,antibiotic)
{
  returnFrame <- frame
  if(!is.null(returnFrame))  {
    returnFrame <- as.data.frame(apply(returnFrame,c(1,2),FUN = function(value)
    {
      words <- unlist(strsplit(value," "))
      sievedWords <- words[grep(antibiotic,words,invert=TRUE)]
      value = paste(sievedWords,collapse=" ")
      value
    }))
  }
  returnFrame
}



ALL<-function(){
  stopifnot(dir.exists(IN_FOLDER))
  
  
  if(!dir.exists(OUT_FOLDER)){
    dir.create(OUT_FOLDER,recursive = TRUE )
  }
  SIEVE(antibiotic = "AMC")
}

SIEVE <- function(antibiotic="AMC")
{
  files <- list.files(IN_FOLDER,pattern = "modelOutput_")
  for (file in files){
    print(file)
    frame <- readOutputFrame(file)
    print(dim.data.frame(frame))
    frame <- sieveValues(sieveColumns(frame,antibiotic),antibiotic)
    if(!is.null(frame)){
      writeOutputFrame(frame,file)
    } else {
      print("skip")
    }
  }
}




