source(file='common.R')

dir <- tygsDirectory
outdir <- paste(processedRootR, "TYGS", sep="/")




samples <- list.files(dir,"Type Strain Genome Server.csv",recursive = TRUE)
sample <- samples[1]

mergedTable <- data.frame()
for (sample in samples) {
  aTable <- read.csv(comment.char = "#",file = paste(sep="",dir, "\\" , sample),sep=",")
  if(nrow(mergedTable) == 0) {
    mergedTable <- aTable
  }
  else {
    mergedTable <- rbind(mergedTable,aTable)
  }
}

## FIXME sort mergedTable2 <- order(x=mergedTable)

write.csv2(x=mergedTable,row.names = TRUE,file= paste(sep="",outdir, "\\" ,"fullMergedTable.csv"))

orderedMergedTable <- mergedTable[order(mergedTable[,1]),]
write.csv2(x=orderedMergedTable,row.names = TRUE,file= paste(sep="",outdir, "\\" ,"fullMergedTable.csv"))

filteredMergedTable <- orderedMergedTable[grep("sample",orderedMergedTable[,2],invert=TRUE),]
write.csv2(x=filteredMergedTable,row.names = TRUE,file= paste(sep="",outdir, "\\" ,"referenceMergedTable.csv"))

#rownames(filteredMergedTable) <- filteredMergedTable[,1]   



samplesTYGS <- sort(unique(mergedTable$Query.strain[grep("sample",mergedTable$Query.strain)]))
#exclude 14
samplesStart <- samplesTYGS[samplesTYGS != "'sample14'"]

eColi <- orderedMergedTable[grep("Escherichia coli DSM 30083",orderedMergedTable[,2],invert=FALSE),]
eColi <- eColi[eColi[,1] %in% samplesStart ,]
eColi <- eColi[order(eColi[,5],decreasing = TRUE),]
sSonneri <- orderedMergedTable[grep("Shigella sonnei ATCC 29930",orderedMergedTable[,2],invert=FALSE),]
sSonneri <- sSonneri[sSonneri[,1] %in% samplesStart ,]
sSonneri <- sSonneri[order(sSonneri[,5],decreasing = TRUE),]


# Put 50 in different buckets, depending on d4-closeness t
cladeA <- {}
cladeB <- {}
samplesLeft <- samplesStart
#Order by d4
eColiLeft <- eColi
sSonneriLeft <- sSonneri
length(samplesLeft)
odd <- FALSE

while(length(samplesLeft)>0){
  sample <- NULL
  if(!odd){
    sample <- eColiLeft[1,1]
    cladeA <- c(cladeA, sample)
  } else {
    sample <- sSonneriLeft[1,1]
    cladeB <- c(cladeB, sample)
  }
  eColiLeft <- eColiLeft[eColiLeft[,1] != sample ,]
  sSonneriLeft <- sSonneriLeft[sSonneriLeft[,1] != sample ,]
  samplesLeft <- samplesLeft[samplesLeft != sample]
  odd <- !odd
}

cladeATable <- eColi[eColi[,1] %in% cladeA,]
write.csv2(x=cladeATable,row.names = TRUE,file= paste(sep="",outdir, "\\" ,"cladeATable.csv"))
cladeBTable <- sSonneri[sSonneri[,1] %in% cladeB,]
write.csv2(x=cladeBTable,row.names = TRUE,file= paste(sep="",outdir, "\\" ,"cladeBTable.csv"))



sampleNamesTYGS <- substring(samplesTYGS,2,nchar(samplesTYGS)-1)
print(sampleNamesTYGS)


