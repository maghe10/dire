library(readr)
library(dplyr)
library(openxlsx) 
source(file='common.R')
dir <- amrfinderDirectory
outdir <- paste(processedRootRassembly, "amrfinder", sep="/")


samples <- list.files(dir)

# remove excluded
samples <- samples[samples != "sample14.tsv" ]
samples <- samples[samples != "sample38.tsv" ]

#sample <- samples[1]
#anotherSample <- samples[2]
#sample <- anotherSample

summaryColumns <- c('Gene symbol','Sequence name','Scope','Element type','Element subtype','Class','Subclass','Method','Accession of closest sequence')
allSummary <- data.frame()

AMRs <- {}
for (sample in samples) {
  aSample <- read_tsv( paste(sep="",dir, "\\" ,sample) ,col_names=TRUE)
  aSampleAMR <- aSample[aSample$`Element type` == 'AMR',]
  aSampleAMRSummary <- aSampleAMR[,summaryColumns]

  AMRs <- union(aSampleAMR$`Gene symbol`,AMRs)

  toAppend <-aSampleAMRSummary[!aSampleAMRSummary$`Gene symbol` %in% allSummary$`Gene symbol` ,]
  allSummary <- rbind (allSummary,toAppend) 
}

allSummary <- allSummary %>% 
  arrange(across(everything()))


geneTable <- data.frame(AMRs)
colnames(geneTable)[ncol(geneTable)] <- 'Gene symbol'

for (sample in samples) {
  aSample <- read_tsv( paste(sep="",dir, "\\" ,sample) ,col_names=TRUE)
  aSampleAMR <- aSample[aSample$`Element type` == 'AMR',]
  
  geneTable <- cbind(geneTable,AMRs %in% aSampleAMR$`Gene symbol`)
  colnames(geneTable)[ncol(geneTable)] <- sample
}



statistics = data.frame(geneTable$`Gene symbol`,rowSums(geneTable[-1]))
colnames(statistics)[1] <- 'Gene symbol' 
colnames(statistics)[2] <- 'Count'



sampleNamesAMR <- sort(colnames(geneTable[,c(2:100)]))
sampleNamesAMR <- substring(sampleNamesAMR,1,nchar(sampleNamesAMR)-4)

AMRClasses <- unique(allSummary$Class)
AMRClass <- AMRClasses[4]
ClassTable <- data.frame()
ClassCoreTable <- data.frame()
for (AMRClass in AMRClasses){
  genesFrameMathcingAMRClass <- allSummary[allSummary$`Class`== AMRClass,]
  genesCoreFrameMathcingAMRClass <- genesFrameMathcingAMRClass[genesFrameMathcingAMRClass$`Scope`== 'core',]
  
  aGeneSubtable <- geneTable
  aGeneCoreSubtable <- geneTable
  aGeneSubtable <-aGeneSubtable[which(aGeneSubtable$`Gene symbol` %in% genesFrameMathcingAMRClass$`Gene symbol`),]
  aGeneCoreSubtable <-aGeneSubtable[which(aGeneSubtable$`Gene symbol` %in% genesCoreFrameMathcingAMRClass$`Gene symbol`),]
  
  row <- data.frame(AMRClass)
  rowCore <- data.frame(AMRClass)
  for (n in 2:ncol(aGeneSubtable)){
    bVal <- any(aGeneSubtable[,n])
    bValCore <- any(aGeneCoreSubtable[,n])
#    print(bVal)
    row = cbind(row,bVal)
    rowCore = cbind(rowCore,bValCore)
  }
  colnames(row) <- colnames(geneTable)
  colnames(rowCore) <- colnames(geneTable)
  if(nrow(ClassTable) == 0) {
    ClassTable <- data.frame(row)
    colnames(ClassTable) <- colnames(geneTable)
  }
  else {
    ClassTable <- rbind(ClassTable,row)
  }
  if(nrow(ClassCoreTable) == 0) {
    ClassCoreTable <- data.frame(rowCore)
    colnames(ClassCoreTable) <- colnames(geneTable)
  }
  else {
    ClassCoreTable <- rbind(ClassCoreTable,rowCore)
  }
  
}
colnames(ClassTable)[1] <- 'Class'
colnames(ClassCoreTable)[1] <- 'Class'



statistics <- statistics[order(statistics$`Gene symbol`),]
             

write.csv2(x= geneTable,row.names = FALSE,file= paste(sep="",outdir, "\\" ,assemblymethod,"_genes_",amrfinderDatabase,".csv"))
write.csv2(x= statistics,row.names = FALSE,file= paste(sep="",outdir, "\\" ,assemblymethod,"_statistics_",amrfinderDatabase,".csv"))
write.csv2(x= allSummary,row.names = FALSE,file= paste(sep="",outdir, "\\" ,assemblymethod,"_summary_",amrfinderDatabase,".csv"))
write.csv2(x= ClassTable,row.names = FALSE,file= paste(sep="",outdir, "\\" ,assemblymethod,"_class_",amrfinderDatabase,".csv"))
write.csv2(x= ClassCoreTable,row.names = FALSE,file= paste(sep="",outdir, "\\" ,assemblymethod,"_classcore_",amrfinderDatabase,".csv"))

write.xlsx(x= geneTable,rowNames = FALSE,file= paste(sep="",outdir, "\\" ,assemblymethod,"_genes_",amrfinderDatabase,".xlsx"))
write.xlsx(x= statistics,rowNames = FALSE,file= paste(sep="",outdir, "\\" ,assemblymethod,"_statistics_",amrfinderDatabase,".xlsx"))
write.xlsx(x= allSummary,rowNames = FALSE,file= paste(sep="",outdir, "\\" ,assemblymethod,"_summary_",amrfinderDatabase,".xlsx"))
write.xlsx(x= ClassTable,rowNames = FALSE,file= paste(sep="",outdir, "\\" ,assemblymethod,"_class_",amrfinderDatabase,".xlsx"))
write.xlsx(x= ClassCoreTable,rowNames = FALSE,file= paste(sep="",outdir, "\\" ,assemblymethod,"_classcore_",amrfinderDatabase,".xlsx"))


