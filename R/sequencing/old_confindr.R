source(file='common.R')



dir <- confindrtrimmedDirectory
outdir <- paste(processedRootR, "confindrtrimmed", sep="/")


samplesReport <- list.files(dir,"confindr_report.csv",recursive = TRUE)
samplesAlleles <- list.files(dir,"*_alleles.csv",recursive = TRUE)
samplesContamination <- list.files(dir,"*_contaminaion.csv",recursive = TRUE)


#sample <- samplesReport[1]

allSamplesReport <- data.frame()

for (sample in samplesReport) {
  aSplit<- strsplit(sample, split = "/")[[1]]  
  aSampleName <- aSplit[[1]]
  print(aSampleName)
  aTable <- read.csv(comment.char = "#",file = paste(sep="",dir, "\\" , sample),sep=",")
  if(nrow(allSamplesReport) == 0) {
    allSamplesReport <- aTable
  }
  else {
    allSamplesReport <- rbind(allSamplesReport,aTable)
  }
}




allSamplesTableAlleles <- data.frame()

for (sample in samplesAlleles) {
  aSplit<- strsplit(sample, split = "/")[[1]]  
  aSampleName <- aSplit[[1]]
  aTable <- read.csv(comment.char = "#",file = paste(sep="",dir, "\\" , sample),sep=",")
  if(nrow(aTable)!=0) {
  listOfSampleName <- rep(aSampleName,nrow(aTable))
  aTableWithSampleName <- cbind(listOfSampleName,aTable)
  colnames(aTableWithSampleName)[1] <- "sample"
  
  if(nrow(allSamplesTableAlleles) == 0) {
    allSamplesTableAlleles <- aTableWithSampleName
  }
  else {
    allSamplesTableAlleles <- rbind(allSamplesTableAlleles,aTableWithSampleName)
  }
  }
}

allSamplesTableContamination <- data.frame()
for (sample in samplesAlleles) {
  aSplit<- strsplit(sample, split = "/")[[1]]  
  aSampleName <- aSplit[[1]]
  aTable <- read.csv(comment.char = "#",file = paste(sep="",dir, "\\" , sample),sep=",")
  if(nrow(aTable)!=0) {
    listOfSampleName <- rep(aSampleName,nrow(aTable))
    aTableWithSampleName <- cbind(listOfSampleName,aTable)
    colnames(aTableWithSampleName)[1] <- "sample"
    
    if(nrow(allSamplesTableContamination) == 0) {
      allSamplesTableContamination <- aTableWithSampleName
    }
    else {
      allSamplesTableContamination <- rbind(allSamplesTableContamination,aTableWithSampleName)
    }
  }
}


write.csv2(x=allSamplesTableContamination,row.names = FALSE,file= paste(sep="",outdir, "\\" ,"allSamplesContaminationConfindr.csv"))
write.csv2(x=allSamplesTableAlleles,row.names = FALSE,file= paste(sep="",outdir, "\\" ,"allSamplesTableAllelesConfindr.csv"))
write.csv2(x=allSamplesReport,row.names = FALSE,file= paste(sep="",outdir, "\\" ,"allSamplesReportConfindr.csv"))


