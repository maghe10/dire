source(file = 'manuscript/manuscriptcommon.R')


amrStatistics <-function()
{
  indir <- paste(processedRootRassembly, "amrfinder", sep="/")
  aList <- list()
  
  frame <- read.xlsx(xlsxFile = paste(sep="",indir, "\\" ,assemblymethod,"_allSamplesAndAllGenes_",amrfinderDatabase,".xlsx"))
  colnames(frame) <- gsub("\\.", " ", colnames(frame))
  aList[["allSamplesAndAllGenes"]] <- frame
  
  frame <- read.xlsx(xlsxFile = paste(sep="",indir, "\\" ,assemblymethod,"_classcore-in-samples_",amrfinderDatabase,".xlsx"))
  colnames(frame) <- gsub("\\.", " ", colnames(frame))
  aList[["classCoreInSamples"]] <- frame
  
  frame <- read.xlsx(xlsxFile = paste(sep="",indir, "\\" ,assemblymethod,"_genes-in-samples_",amrfinderDatabase,".xlsx"))
  colnames(frame) <- gsub("\\.", " ", colnames(frame))
  aList[["allGenesInSamples"]] <- frame
  
  frame <- read.xlsx(xlsxFile = paste(sep="",indir, "\\" ,assemblymethod,"_betalactam-core-genes-in-samples_",amrfinderDatabase,".xlsx"))
  colnames(frame) <- gsub("\\.", " ", colnames(frame))
  aList[["betalactamCoreInSamples"]] <- frame
  
  
  writeStatisticsExcel(aList,"amrStatistics","amr") 
}

compareGenotypePhenotype <- function()
{
  indir <- paste(processedRootRassembly, "genotype", sep="/")
  names <- c("betalactamases","genotype_phenotype","enzymeFamily_phenotype", "samples_and_genotypes")
  for(name in names){
    frame <- read.csv2(check.names=FALSE,file= paste(sep="",indir, "\\" ,assemblymethod,"_",name,"_",amrfinderDatabase,".csv"))
    writeStatisticsExcel(frame=frame,name,name)
  }
  
}


ALL <-function()
{
  compareGenotypePhenotype()
  amrStatistics()
}

