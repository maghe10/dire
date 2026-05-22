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

sampleNameToNumber <- function(sampleName)
{
  as.integer(stringr::str_extract(sampleName, "(\\d+)"))
}


# Older version summaryColumns <- c('Gene symbol','Sequence name','Scope','Element type','Element subtype','Class','Subclass','Method','Accession of closest sequence')
oldSummaryColumns <- c('Gene symbol','Sequence name','Scope','Element type','Element subtype','Class','Subclass','Method','Accession of closest sequence')
summaryColumns <- c('Element symbol','Element name','Scope','Type','Subtype','Class','Subclass','Method','Closest reference accession')
allSummary <- data.frame()

#geneList <- c("pmrB_E123D","glpT_E448K")
#gene <- geneList[2]
#sample <- samples[2]
#aSample <- read_tsv( paste(sep="",dir, "\\" ,sample) ,col_names=TRUE)
#sampleGeneRow <- cbind(sample=sample,aSample[aSample$`Gene symbol`==gene,])

allSamplesAndAllGenes <- data.frame()
AMRs <- {}
for (sample in samples) {
  # sample <- "sample100.tsv"
  print(sample)
  aSample <- read_tsv( paste(sep="",dir, "\\" ,sample) ,col_names=TRUE)
  #older version aSampleAMR <- aSample[aSample$`Element type` == 'AMR',]
  aSampleAMR <- aSample[aSample$`Type` == 'AMR',]
  aSampleAMRSummary <- aSampleAMR[,summaryColumns]
  colnames(aSampleAMRSummary) <- oldSummaryColumns
  
  aSampleWithSampleId <- cbind(sampleid=sampleNameToNumber(sample),aSample)
  allSamplesAndAllGenes <- rbind(allSamplesAndAllGenes,aSampleWithSampleId)

  AMRs <- union(aSampleAMRSummary$`Gene symbol`,AMRs)

  toAppend <-aSampleAMRSummary[!aSampleAMRSummary$`Gene symbol` %in% allSummary$`Gene symbol` ,]
  allSummary <- rbind (allSummary,toAppend) 
}

allSamplesAndAllGenes <- allSamplesAndAllGenes[order(allSamplesAndAllGenes$sampleid),]

allSummary <- allSummary %>% 
  arrange(across(everything()))



names(samples) <- unlist(sampleNameToNumber(samples))
sortedSamples <- samples[order(as.integer(names(samples)))]

geneTable <- data.frame(AMRs)
colnames(geneTable)[ncol(geneTable)] <- 'Gene symbol'
for (sample in sortedSamples) {
  aSample <- read_tsv( paste(sep="",dir, "\\" ,sample) ,col_names=TRUE)
  aSampleAMR <- aSample[aSample$`Type` == 'AMR',]
  
  geneTable <- cbind(geneTable,AMRs %in% aSampleAMR$`Element symbol`)
  colnames(geneTable)[ncol(geneTable)] <- sampleNameToNumber(sample)
}


statistics = data.frame(geneTable$`Gene symbol`,rowSums(geneTable[-1]))
colnames(statistics)[1] <- 'Gene symbol' 
colnames(statistics)[2] <- 'Count'

sampleNamesAMR <- colnames(geneTable[,c(2:ncol(geneTable))])

AMRClasses <- unique(allSummary$Class)
#AMRClass <- AMRClasses[4]
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
geneTable <- geneTable[order(geneTable$`Gene symbol`),]

betalactamGenes <- allSummary[allSummary$Class=="BETA-LACTAM",]$`Gene symbol`
betalactamCoreGenes <- allSummary[allSummary$Class=="BETA-LACTAM" & allSummary$Scope=="core",]$`Gene symbol`
betalactamPlusGenes <- allSummary[allSummary$Class=="BETA-LACTAM" & allSummary$Scope=="plus",]$`Gene symbol`
otherCoreGenes <- allSummary[allSummary$Class!="BETA-LACTAM" & allSummary$Scope=="core",]$`Gene symbol`
otherPlusGenes <- allSummary[allSummary$Class!="BETA-LACTAM" & allSummary$Scope=="plus",]$`Gene symbol`

betalactamGeneTable <- geneTable[geneTable$`Gene symbol` %in% betalactamGenes,]
betalactamCoreGeneTable <- geneTable[geneTable$`Gene symbol` %in% betalactamCoreGenes,]
betalactamPlusGeneTable <- geneTable[geneTable$`Gene symbol` %in% betalactamPlusGenes,]
otherCoreGeneTable <- geneTable[geneTable$`Gene symbol` %in% otherCoreGenes,]
otherPlusGeneTable <- geneTable[geneTable$`Gene symbol` %in% otherPlusGenes,]

betalactamCoreSummary <- allSummary[allSummary$Class=="BETA-LACTAM" & allSummary$Scope=="core",]


WRITE_TABLES <- function()
{
  write.xlsx(x= geneTable,rowNames = FALSE,file= paste(sep="",outdir, "\\" ,assemblymethod,"_genes-in-samples_",amrfinderDatabase,".xlsx"))
  write.xlsx(x= betalactamGeneTable,rowNames = FALSE,file= paste(sep="",outdir, "\\" ,assemblymethod,"_betalactam-genes-in-samples_",amrfinderDatabase,".xlsx"))
  write.xlsx(x= betalactamCoreGeneTable,rowNames = FALSE,file= paste(sep="",outdir, "\\" ,assemblymethod,"_betalactam-core-genes-in-samples_",amrfinderDatabase,".xlsx"))
  write.xlsx(x= betalactamPlusGeneTable,rowNames = FALSE,file= paste(sep="",outdir, "\\" ,assemblymethod,"_betalactam-plus-genes-in-samples_",amrfinderDatabase,".xlsx"))
  
  write.xlsx(x= statistics,rowNames = FALSE,file= paste(sep="",outdir, "\\" ,assemblymethod,"_gene-statistics_",amrfinderDatabase,".xlsx"))
  write.xlsx(x= allSummary,rowNames = FALSE,file= paste(sep="",outdir, "\\" ,assemblymethod,"_gene-summary_",amrfinderDatabase,".xlsx"))
  write.xlsx(x= betalactamCoreSummary,rowNames = FALSE,file= paste(sep="",outdir, "\\" ,assemblymethod,"_betalactam-core-summary_",amrfinderDatabase,".xlsx"))
  write.xlsx(x= ClassTable,rowNames = FALSE,file= paste(sep="",outdir, "\\" ,assemblymethod,"_class-in-samples_",amrfinderDatabase,".xlsx"))
  write.xlsx(x= ClassCoreTable,rowNames = FALSE,file= paste(sep="",outdir, "\\" ,assemblymethod,"_classcore-in-samples_",amrfinderDatabase,".xlsx"))
  
  write.xlsx(x= allSamplesAndAllGenes,rowNames = FALSE,file= paste(sep="",outdir, "\\" ,assemblymethod,"_allSamplesAndAllGenes_",amrfinderDatabase,".xlsx"))
}



ALL <- function()
{
  WRITE_TABLES()
}
