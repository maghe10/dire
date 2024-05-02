source(file='common.R')



dir <- jspecieswsDirectory
outdir <- paste(processedRootR, "jspeciesws", sep="/")

#fixme Csv
#samples <- list.files(dir,"^A[\\w]+\\.csv$")



samples <- list.files(dir,"^ANIbMA.")

allSamplesVsRef <- data.frame()
allRefVsSamples <- data.frame()


for (sample in samples) {
  aTable <- read.csv(comment.char = "#",file = paste(sep="",dir, "\\" , sample),sep="\t",skip=1)

  samplesVsRef <- aTable[grep("sample",aTable$X), ]
  rownames(samplesVsRef) <- samplesVsRef$X
  samplesVsRef <- data.frame(samplesVsRef[grep("esche",colnames(samplesVsRef))])
  if(nrow(allSamplesVsRef) == 0) {
    allSamplesVsRef <- samplesVsRef
  } else {
    allSamplesVsRef <- rbind(allSamplesVsRef,samplesVsRef)
  }
  refVsSamples <- aTable[grep("esche",aTable$X), ]
  rownames(refVsSamples) <- refVsSamples$X
  refVsSamples <- data.frame(refVsSamples[grep("sample",colnames(refVsSamples))])
  if(nrow(allRefVsSamples) == 0) {
    allRefVsSamples <- refVsSamples
  }
  else {
    allRefVsSamples <- cbind(allRefVsSamples,refVsSamples)
  }
}

value <- allRefVsSamples[1,]
#substring(allRefVsSamples,1,5)
#substring(allRefVsSamples,8,11)

allSamplesVsRef <- cbind(allSamplesVsRef, as.numeric(substring(allSamplesVsRef[,1],1,5)))
allSamplesVsRef <- cbind(allSamplesVsRef, as.numeric(substring(allSamplesVsRef[,1],8,12)))
colnames(allSamplesVsRef)[2] <- "ANIb"
colnames(allSamplesVsRef)[3] <- "Aligned"

allRefVsSamples <- rbind(allRefVsSamples, as.numeric(substring(allRefVsSamples[1,],1,5)))
allRefVsSamples <- rbind(allRefVsSamples, as.numeric(substring(allRefVsSamples[1,],8,12)))
rownames(allRefVsSamples)[2] <- "ANIb"
rownames(allRefVsSamples)[3] <- "Aligned"

#Excludes

allSamplesVsRef <-allSamplesVsRef[rownames(allSamplesVsRef)!="sample38.fasta",]
allSamplesVsRef <-allSamplesVsRef[rownames(allSamplesVsRef)!="sample14.fasta",]
allRefVsSamples <- allRefVsSamples[,colnames(allRefVsSamples)!="sample38.fasta"]
allRefVsSamples <- allRefVsSamples[,colnames(allRefVsSamples)!="sample14.fasta"]

#Fix dot to comma
allRefVsSamples <- format(allRefVsSamples, decimal.mark = ',')
allSamplesVsRef <- format(allSamplesVsRef, decimal.mark = ',')



write.csv2(x=allSamplesVsRef,row.names = TRUE,file= paste(sep="",outdir, "\\" ,"allSamplesVsRef.csv"))
write.csv2(x=allRefVsSamples,row.names = TRUE,file= paste(sep="",outdir, "\\" ,"allRefVsSamples.csv"))

print(min(as.numeric(allSamplesVsRef$ANIb)))
print(min( as.numeric(allRefVsSamples[2,])))

print(min(as.numeric(allSamplesVsRef$Aligned)))
print(min( as.numeric(allRefVsSamples[3,])))

sampleNamesJWS <- sort(rownames(allSamplesVsRef))
sampleNamesJWS <- substring(sampleNamesJWS,1,nchar(sampleNamesJWS)-6)
print(sampleNamesJWS)

unique(rownames(allSamplesVsRef))


