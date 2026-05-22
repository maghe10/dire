source(file='common.R')
library(openxlsx) 

dir <- resfinderDirectory
outdir <- paste(processedRootRassembly, "resfinder", sep="/")

subdirs <- list.files(dir,pattern = "^sample[0-9]+")

#excluded
subdirs <- subdirs[subdirs != "sample14" ]
subdirs <- subdirs[subdirs != "sample38" ]


sampleNameToNumber <- function(sampleName)
{
  as.integer(stringr::str_extract(sampleName, "(\\d+)"))
}


#subdirs <- c("sample2","sample3")
#subdir <- "sample2"

readPhenoTable <- function(subdir)
{
    f <- read_delim(file = paste(dir,subdir,"pheno_table_escherichia_coli.txt",sep="/"),  delim ="\t",skip = 16)
    f
}

readGeneTable <- function(subdir)
{
  f <- read_delim(file = paste(dir,subdir,"Resfinder_results_tab.txt",sep="/"),  delim ="\t",skip = 0)
  f
}


allGenes <- function()
{
  allGenes <- c()
  for(subdir in subdirs){
    f <- readGeneTable(subdir)
    genes <- f$`Resistance gene`
    allGenes <- unique(append(allGenes,genes))
  }
  sort(allGenes)
}
  

fetchAggregatedPhenoTable <- function()
{
  result <- NULL
  
  for(subdir in subdirs){
    f <- readPhenoTable(subdir)
    sampleID <- sampleNameToNumber(subdir)
    
    AMP <-f[f$`# Antimicrobial`=='ampicillin',]$`WGS-predicted phenotype`
    AMC <-f[f$`# Antimicrobial`=='ampicillin+clavulanic acid',]$`WGS-predicted phenotype`
    TZP  <-  f[f$`# Antimicrobial`=='piperacillin+tazobactam',]$`WGS-predicted phenotype`
    CTX <-  f[f$`# Antimicrobial`=='cefotaxime',]$`WGS-predicted phenotype`
    CAZ <-  f[f$`# Antimicrobial`=='cefotaxime',]$`WGS-predicted phenotype`
    FEP <-  f[f$`# Antimicrobial`=='cefepime',]$`WGS-predicted phenotype`
    AMPGene <-f[f$`# Antimicrobial`=='ampicillin',]$`Genetic background`
    AMCGene <-f[f$`# Antimicrobial`=='ampicillin+clavulanic acid',]$`Genetic background`
    row <- list(sampleID,AMP,AMC,TZP,CTX,CAZ,FEP,AMPGene,AMCGene)
    if(is.null(result)){
      result <- data.frame(row)
      colnames(result) <- c("sampleID","AMP","AMC","TZP","CTX","CAZ","FEP","AMPGene","AMCGene")
    } else {
      result <- rbind(result,row)
    }
  }
  sampleSort(result)    
}

fetchAggregatedGenoTable <- function(allGenes=allGenes())
{
  result <- data.frame(matrix(ncol = length(allGenes)+1, nrow = 0))
  
  
  for(subdir in subdirs){
    f <- readGeneTable(subdir)
    sampleID <- sampleNameToNumber(subdir)
    row <- as.character(allGenes %in% f$`Resistance gene`)
    row <- append(row,sampleID,after=0)
    result <- rbind(result,row)
  }
  colnames(result) <- append(allGenes,"sampleID",after = 0)
  sampleSort(result)    
}



sampleSort <-function(aTable)
{
  aTable[order(as.integer(aTable$sampleID)),]
}




readAmrFinderGeneTable <-function()
{
  amrFinderGenes <- read.xlsx(xlsxFile= paste(sep="",paste(processedRootRassembly, "amrfinder", sep="/"), "\\" ,assemblymethod,"_betalactam-core-genes-in-samples_",amrfinderDatabase,".xlsx"))
  sampleID <- colnames(amrFinderGenes[-1])
  fixedAmrFinderGenes <- t(amrFinderGenes)
  colnames(fixedAmrFinderGenes) <- fixedAmrFinderGenes[1,]
  fixedAmrFinderGenes <- fixedAmrFinderGenes[-1,]
  fixedAmrFinderGenes <- cbind(sampleID,fixedAmrFinderGenes)
  as.data.frame(fixedAmrFinderGenes)
}

isBetalactamaseGene <-function(genename)
{
  #genename = "blaTEM"
  substr(genename,1,3)=="bla"
}

expandWithFalse <- function(fromTable, expandGenes)
{
  for(gene in expandGenes) {
    fromTable <- cbind(fromTable,rep(FALSE,nrow(fromTable)))
    colnames(fromTable)[ncol(fromTable)] <- gene
  }
  fromTable[,order(colnames(fromTable))]
}

diffAmrFinderGenes <- function()
{
  resFinderGeneTable <- read.xlsx(xlsxFile= paste(sep="",outdir, "\\" ,assemblymethod,"_resfinder_genotype_",resfinderDatabase,".xlsx"))
  amrFinderGeneTable <- readAmrFinderGeneTable()
  resFinderGenes <- colnames(resFinderGeneTable[-1])
#  print(length(resFinderGenes))
  amrFinderGenes <- colnames(amrFinderGeneTable[-1])
#  print(length(amrFinderGenes))

 betalactamaseAmrFinderGenes <- amrFinderGenes
  
 #betalactamaseAmrFinderGenes <- amrFinderGenes[unlist(lapply(amrFinderGenes,isBetalactamaseGene))]
  betalactamaseResFinderGenes <- resFinderGenes[unlist(lapply(resFinderGenes,isBetalactamaseGene))]

  onlyAmrFinderGenes = setdiff(betalactamaseAmrFinderGenes,betalactamaseResFinderGenes)
  onlyResFinderGenes = setdiff(betalactamaseResFinderGenes,betalactamaseAmrFinderGenes)
  
  
  resFinderBetalactamaseTable <- expandWithFalse(resFinderGeneTable[,betalactamaseResFinderGenes],onlyAmrFinderGenes)
  amrFinderBetalactamaseTable <- expandWithFalse(amrFinderGeneTable[,betalactamaseAmrFinderGenes],onlyResFinderGenes)
  
 # colnames(resFinderBetalactamaseTable)
#  colnames(amrFinderBetalactamaseTable)
#  combinedGenes
#  dim(resFinderBetalactamaseTable)
#  dim(amrFinderBetalactamaseTable)
  onlyResFinder = resFinderBetalactamaseTable==TRUE & amrFinderBetalactamaseTable==FALSE
  onlyAmrFinder = amrFinderBetalactamaseTable==TRUE & resFinderBetalactamaseTable==FALSE
#  dim(onlyAmrFinder)

  amrFinder <- apply(
    amrFinderBetalactamaseTable=="TRUE", 1, 
    function(u) paste( names(which(u)), collapse="," ) 
  )
  resFinder <- apply(
    resFinderBetalactamaseTable, 1, 
    function(u) paste( names(which(u)), collapse="," ) 
  )
  
  amrFinderOnly <- apply(
    onlyAmrFinder, 1, 
    function(u) paste( names(which(u)), collapse="," ) 
  )
  resFinderOnly <- apply(
    onlyResFinder, 1, 
    function(u) paste( names(which(u)), collapse="," ) 
  )
  
  
  listDiff = list()
  for(i in 1:length(amrFinderOnly)){
    aDiff = ""
    if(amrFinderOnly[i]!=""){
      aDiff = paste("+",amrFinderOnly[i],sep="")
    }
    if(resFinderOnly[i]!=""){
      aDiff = paste(aDiff,"-",resFinderOnly[i],sep="")
    }
    listDiff = append(listDiff,aDiff)
  } 
  
  result <- data.frame(matrix(ncol = 0, nrow = length(amrFinderOnly)))
  result <- cbind(result,resFinderGeneTable$sampleID)
  colnames(result)[ncol(result)] <- "sampleID"
  result <- cbind(result,unlist(listDiff))
  colnames(result)[ncol(result)] <- "diff"
  result <- cbind(result,amrFinderOnly)
  result <- cbind(result,resFinderOnly)
  result <- cbind(result,amrFinder)
  result <- cbind(result,resFinder)
  
  result
}


WRITE_TABLES <- function()
{
  
  write.xlsx(x= fetchAggregatedPhenoTable(),rowNames = FALSE,file= paste(sep="",outdir, "\\" ,assemblymethod,"_resfinder_phenotype_",resfinderDatabase,".xlsx"))
  allGenes <- allGenes()
  write.xlsx(x= fetchAggregatedGenoTable(allGenes),rowNames = FALSE,file= paste(sep="",outdir, "\\" ,assemblymethod,"_resfinder_genotype_",resfinderDatabase,".xlsx"))
  write.xlsx(x = readAmrFinderGeneTable(),,rowNames = FALSE,file= paste(sep="",outdir, "\\" ,assemblymethod,"_amrfinder_genotype_",resfinderDatabase,".xlsx"))  

  write.xlsx(x = diffAmrFinderGenes(),rowNames = FALSE,file= paste(sep="",outdir, "\\" ,assemblymethod,"_amrfinder_resfinder_diff_",resfinderDatabase,".xlsx"))  
  
}


ALL <- function()
{
  if(!dir.exists(outdir)){
    dir.create(outdir)
  }
  WRITE_TABLES()
}


