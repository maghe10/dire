#source(file='genotype/amrfinder.R')
source(file='model/modelCommon.R')

outdir <- paste(processedRootRassembly, "genotype", sep="/")


GENE_SYMBOL <- "Gene_symbol"
CLAVULANIC_ACID_EXPECTED_EFECTIVE <- "Clavulanic_acid_expected_effective"
MOLECULAR_CLASS <- "Molecular_class"
FUNCTIONAL_GROUP <- "Functional_group"

################### curate a bit #############################
getCuratedGeneTable <- function(orgGeneTable)
{
  tmp <- orgGeneTable
  #  dim(tmp)
  
  tmp[tmp$`Gene symbol`=="blaCTX-M-15",-1] <- tmp[tmp$`Gene symbol`=="blaCTX-M",-1] | tmp[tmp$`Gene symbol`=="blaCTX-M-15",-1]
  tmp[tmp$`Gene symbol`=="blaTEM-1",-1] <- tmp[tmp$`Gene symbol`=="blaTEM",-1] | tmp[tmp$`Gene symbol`=="blaTEM-1",-1]
  tmp <- tmp[tmp$`Gene symbol`!="blaCTX-M",]
  tmp <- tmp[tmp$`Gene symbol`!="blaTEM",]
  #  dim(tmp)
  tmp
}

readGeneTable <- function()
{
  outdir <- paste(processedRootRassembly, "amrfinder", sep="/")
  betalactamCoreGeneTable <- read.xlsx(check.names = FALSE,sep.names = " ",xlsxFile= paste(sep="",outdir, "\\" ,assemblymethod,"_betalactam-core-genes-in-samples_",amrfinderDatabase,".xlsx"))
  betalactamCoreGeneTable
}

readBetalactamasesTable <- function()
{
  read.xlsx(paste(processedRootExcel,"Betalactamases.xlsx",sep="\\"),check.names=FALSE,sep.names=" ")
}


orgGeneTable <- readGeneTable()
#orgGeneTable[orgGeneTable$`Gene symbol`=="blaTEM",colnames(orgGeneTable)=="50"]
#orgGeneTable[orgGeneTable$`Gene symbol`=="blaLAP-2",colnames(orgGeneTable)=="50"]

geneTable <- getCuratedGeneTable(orgGeneTable)
geneTable[geneTable$`Gene symbol`=="blaTEM-1",colnames(geneTable)=="50"]
geneTable[geneTable$`Gene symbol`=="blaLAP-2",colnames(geneTable)=="50"]

betalactamCoreGenes <- geneTable$`Gene symbol`
betalactamCoreGeneTable <- geneTable[geneTable$`Gene symbol` %in% betalactamCoreGenes,]

###############################################################

MODE_A <- "Mode-A"
MODE_C <- "Mode-C"

extractFenotype <- function()
{
  mode <- MODE_A
  frameA <- readSirTable(mode)
  mode <- MODE_C
  frameC <- readSirTable(mode)
  
  countTrue <- function(x) {sum(x, na.rm = TRUE)}
  
  SS_NOTATU = frameA$AMP=="S" & frameC$AMC=="S"
  print(countTrue(SS_NOTATU))
  SS_ATU = frameA$AMP=="S" & frameA$AMC=="S" &frameC$AMC=="R"
  print(countTrue(SS_ATU))
  RS_NOTATU = frameA$AMP=="R" & frameC$AMC=="S"
  print(countTrue(RS_NOTATU))
  RS = frameA$AMP=="R" & frameA$AMC=="S" & frameC$AMC=="R"
  print(countTrue(RS))
  RR = frameA$AMP=="R" & frameA$AMC=="R"
  print(countTrue(RR))
  
  aFrame <- data.frame(SS_NOTATU,SS_ATU,RS_NOTATU,RS,RR)
  names(aFrame) <- c ("AMP=S, AMC=S (Not ATU)",	"AMP=S, AMC=S (ATU only)",	"AMP=R, AMC=S (Not ATU)",	"AMP=R, AMC=S (ATU only)",	"AMP=R, AMC=R")
  )
rownames(aFrame)<-frameA$sampleID
aFrame
}





analyzeClavulanicAcidNew <- function()
{
  # preanalysed core betalactamase genes on their effectiveness on clavulanic acid
  betalactamasesTable <-readBetalactamasesTable()
  
  # curate  
  betalactamasesTable <- betalactamasesTable[betalactamasesTable$Gene_symbol!="blaCTX-M",]
  betalactamasesTable <- betalactamasesTable[betalactamasesTable$Gene_symbol!="blaTEM",]
  
  # genotypes
  betalactamaseGenesEffective <- betalactamasesTable[betalactamasesTable[,colnames(betalactamasesTable)==CLAVULANIC_ACID_EXPECTED_EFECTIVE]==TRUE,colnames(betalactamasesTable)==GENE_SYMBOL]
  betalactamaseGenesNonEffective <- betalactamasesTable[betalactamasesTable[,colnames(betalactamasesTable)==CLAVULANIC_ACID_EXPECTED_EFECTIVE]==FALSE,colnames(betalactamasesTable)==GENE_SYMBOL]
  
  # bushjacobi
  bushJacobi <- getBetalactamaseBushJacobiGroupForSample(betalactamasesTable)
  
  # genes
  transposedGenes <- t(geneTable)
  transposedGenes <- data.frame(transposedGenes)
  colnames(transposedGenes) <- transposedGenes[1,]
  transposedGenes <- transposedGenes[-1,]  
  transposedGenes <- apply (transposedGenes,c(1,2),as.logical)
  
  sampleID <- colnames(geneTable[-1])
  resultFrame <- data.frame(sampleID)
  colnames(resultFrame)[length(resultFrame)] <- "sampleID"
  
  # appendGene <- function(aGeneTable,resultFrame) {
  #   for(sample in resultFrame$sampleID){
  #     genesForSample <- aGeneTable$`Gene symbol`[aGeneTable[,colnames(aGeneTable)==sample]]
  #     value <- ""
  #     if(length(genesForSample)>0){
  #       value <- paste(genesForSample,collapse = ",")
  #     }
  #     resultFrame[resultFrame$sampleID==sample,length(resultFrame)] <- value
  #   }
  #   resultFrame
  # }
  
  genotype <- apply(
    transposedGenes, 1, 
    function(u) paste( names(which(u)), collapse="," ) 
  )
  resultFrame <- cbind(resultFrame,genotype)
  colnames(resultFrame)[length(resultFrame)] <- "genotype"
  
  bushJacobiCollapsed <- apply(
    bushJacobi, 1, 
    function(u) paste( names(which(u)), collapse="," ) 
  )
  resultFrame <- cbind(resultFrame,bushJacobiCollapsed)
  colnames(resultFrame)[length(resultFrame)] <- "bush-jacobi"
  
  samplesWithCoreBetalactamasesNonSuceptibleToClavulanicAcidGenetable <- transposedGenes[,geneTable$`Gene symbol` %in% betalactamaseGenesNonEffective]
  samplesWithCoreBetalactamasesNonSuceptibleToClavulanicAcid <- apply(samplesWithCoreBetalactamasesNonSuceptibleToClavulanicAcidGenetable,1,any)
  
  samplesWithoutCoreBetalactamases <- !apply(transposedGenes,1,any)
  samplesWithOnlyCoreBetalactamasesSuceptibleToClavulanicAcid <- !samplesWithCoreBetalactamasesNonSuceptibleToClavulanicAcid & !samplesWithoutCoreBetalactamases
  
  # sum(samplesWithoutCoreBetalactamases)
  # sum(samplesWithCoreBetalactamasesNonSuceptibleToClavulanicAcid)
  # sum(samplesWithOnlyCoreBetalactamasesSuceptibleToClavulanicAcid)
  
  resultFrame <- cbind(resultFrame,samplesWithoutCoreBetalactamases)
  resultFrame <- cbind(resultFrame,samplesWithOnlyCoreBetalactamasesSuceptibleToClavulanicAcid)
  resultFrame <- cbind(resultFrame,samplesWithCoreBetalactamasesNonSuceptibleToClavulanicAcid)
  
  resultFrame <- cbind(resultFrame,bushJacobi)
  resultFrame <- cbind(resultFrame,transposedGenes)
  
  write.csv2(x= resultFrame,row.names = FALSE,file= paste(sep="", outdir, "\\" ,assemblymethod,"_samples_and_genotypes_",amrfinderDatabase,".csv"))
  
  ### now genotypecount
  samplesAndGenotypesFrame <- resultFrame 
  genotypes <- sort(unique(samplesAndGenotypesFrame$genotype))
  
  
  sortedGenotypes <-  genotypes[c(2:length(genotypes),1)]
  
  resultFrame <- data.frame(sortedGenotypes)
  colnames(resultFrame)[length(resultFrame)] <- "genotype"
  
  
  effective <- c()
  for(genotype in resultFrame$genotype){
    # genotype <- "blaTEM-1"
    effectiveGenotype <- all(samplesAndGenotypesFrame[samplesAndGenotypesFrame$genotype==genotype,]$samplesWithOnlyCoreBetalactamasesSuceptibleToClavulanicAcid)
    effective <- append(effective,effectiveGenotype)
  }
  
  resultFrame <- cbind(resultFrame,effective)
  colnames(resultFrame)[length(resultFrame)] <- "Betalactamase genotype"
  
  
  phenotype <- extractFenotype()
  for(name in colnames(phenotype)){
    
    subFrame <- samplesAndGenotypesFrame[phenotype[,name],]
    count <- c()
    for(genotype in resultFrame$genotype){
      # genotype <- "blaTEM-1"
      count <- append(count,sum(subFrame$genotype==genotype))
    }
    resultFrame <- cbind(resultFrame,count)
    colnames(resultFrame)[length(resultFrame)] <- name
  }
  
  write.csv2(x= resultFrame,row.names = FALSE,file= paste(sep="", outdir, "\\" ,assemblymethod,"_genotypes_",amrfinderDatabase,".csv"))
} 


readSirTable <- function(mode = "Mode-A")
{                         
  print(mode)
  millimeterTable <- read.csv2(paste(modelDirectory,"input","millimeterTable.csv",sep="/"),check.names = FALSE)
  sirTable <- read.csv2(paste(modelDirectory,"input",paste("sirAntibioticsModel_", mode ,".csv",sep=""),sep="/"))
  
  allSamples <- colnames(geneTable[-1])
  phenotypicallyAMC_S_AMP_R <-  as.integer(sirTable[sirTable$AMP=="R" & sirTable$AMC=="S",]$sample)
  phenotypicallyAMC_AMP_S <-  as.integer(sirTable[sirTable$AMP=="S" & sirTable$AMC=="S",]$sample)
  phenotypicallyAMC_AMP_R <-   as.integer(sirTable[sirTable$AMP=="R" & sirTable$AMC=="R",]$sample)
  
  
  AMP_HEADER <- "Studie-1-AMP"
  AMC_HEADER <- "KISS I-AMC"
  
  resultFrame <- data.frame(allSamples)
  colnames(resultFrame)[length(resultFrame)] <- "sampleID"
  
  resultFrame <- cbind(resultFrame,rep(NA,99))
  colnames(resultFrame)[length(resultFrame)] <- "AMP"
  resultFrame[,length(resultFrame)] <-  sirTable[,"AMP"]
  
  resultFrame <- cbind(resultFrame,rep(NA,99))
  colnames(resultFrame)[length(resultFrame)] <- "AMP-mm"
  resultFrame[,length(resultFrame)] <- millimeterTable[,AMP_HEADER]
  
  
  resultFrame <- cbind(resultFrame,rep(NA,99))
  colnames(resultFrame)[length(resultFrame)] <- "AMC"
  resultFrame[,length(resultFrame)] <-  sirTable[,"AMC"]
  
  resultFrame <- cbind(resultFrame,rep(NA,99))
  colnames(resultFrame)[length(resultFrame)] <- "AMC-mm"
  resultFrame[,length(resultFrame)] <- millimeterTable[,AMC_HEADER]
  
  resultFrame
}

analyzeClavulanicAcid <- function(mode = "Mode-A")
{
  print(mode)
  # preanalysed core betalactamase genes on their effectiveness on clavulanic acid
  caTable <-read.xlsx(paste(processedRootExcel,"Betalactamases.xlsx",sep="\\"))
  
  # curate  
  caTable <- caTable[caTable$Gene_symbol!="blaCTX-M",]
  caTable <- caTable[caTable$Gene_symbol!="blaTEM",]
  
  # genotypes
  betalactamaseGenesEffective <- caTable[caTable[,colnames(caTable)==CLAVULANIC_ACID_EXPECTED_EFECTIVE]==TRUE,colnames(caTable)==GENE_SYMBOL]
  betalactamaseGenesNonEffective <- caTable[caTable[,colnames(caTable)==CLAVULANIC_ACID_EXPECTED_EFECTIVE]==FALSE,colnames(caTable)==GENE_SYMBOL]
  allBetalactamases <- union(betalactamaseGenesEffective,betalactamaseGenesNonEffective)
  
  potentialEffectiveGeneTable <- geneTable[geneTable$`Gene symbol` %in% betalactamaseGenesEffective,]
  nonEffectiveGeneTable <- geneTable[geneTable$`Gene symbol` %in% betalactamaseGenesNonEffective,]
  
  allBetalactamasesGentable <- geneTable[geneTable$`Gene symbol` %in% allBetalactamases,]
  xxx <- apply(potentialEffectiveGeneTable,2,any)
  xxx <- xxx[-1]
  samplesPotentiallyEffective <- names(xxx[xxx==TRUE])
  
  yyy <- apply(nonEffectiveGeneTable,2,any)
  yyy <- yyy[-1]
  samplesNonEffective <-names(yyy[yyy==TRUE])
  
  zzz <- apply(allBetalactamasesGentable,2,any)
  zzz <- zzz[-1]
  samplesWithBetalactamases <- names(zzz[zzz==TRUE])
  
  
  intersect <- samplesNonEffective[samplesNonEffective %in% samplesPotentiallyEffective]
  # clavulanic acid is assessed to be effective if there is a betalactame that clavulavulanic acid
  # has effect on, but there is no betalactamase that it is not effective on
  samplesEffective <- samplesPotentiallyEffective[!samplesPotentiallyEffective %in% intersect]
  
  genotypicallyPossiblyAMC_S_AMP_R = unlist(lapply(samplesEffective, sampleNameToNumber))
  
  
  
  # phenotypes
  
  millimeterTable <- read.csv2(paste(modelDirectory,"input","millimeterTable.csv",sep="/"))
  sirTable <- read.csv2(paste(modelDirectory,"input",paste("sirAntibioticsModel_", mode ,".csv",sep=""),sep="/"))
  
  allSamples <- colnames(geneTable[-1])
  phenotypicallyAMC_S_AMP_R <-  as.integer(sirTable[sirTable$AMP=="R" & sirTable$AMC=="S",]$sample)
  phenotypicallyAMC_AMP_S <-  as.integer(sirTable[sirTable$AMP=="S" & sirTable$AMC=="S",]$sample)
  phenotypicallyAMC_AMP_R <-   as.integer(sirTable[sirTable$AMP=="R" & sirTable$AMC=="R",]$sample)
  
  #compare phenotype with genotype
  potentialWrongAMC_S_AMP_RGenotype <- phenotypicallyAMC_S_AMP_R[!phenotypicallyAMC_S_AMP_R %in% genotypicallyPossiblyAMC_S_AMP_R]       
  print("Samples with AMC_S_AMP_R phenotype but not genotype")
  print(potentialWrongAMC_S_AMP_RGenotype)
  
  potentialWrongAMC_S_AMP_SGenotype <- phenotypicallyAMC_AMP_S[phenotypicallyAMC_AMP_S %in% samplesWithBetalactamases]       
  print("Samples with AMC_S_AMP_S phenotype but not genotype")
  print(potentialWrongAMC_S_AMP_SGenotype)
  
  
  # equalAMC_S_AMP_R_genotype_phenotype <- intersect(phenotypicallyAMC_S_AMP_R,genotypicallyPossiblyAMC_S_AMP_R)
  # genotypicallyEffectiveWithClavulanicAcidButNotAMC_S_AMP_RPhenotype <- genotypicallyPossiblyAMC_S_AMP_R[!genotypicallyPossiblyAMC_S_AMP_R %in% phenotypicallyAMC_S_AMP_R]       
  
  # geneTableResistantSamples <- geneTable[,colnames(geneTable) %in% genotypicallyEffectiveWithClavulanicAcidButNotAMC_S_AMP_RPhenotype]
  # rownames(geneTableResistantSamples) <- geneTable$`Gene symbol`
  # geneTableResistantSamples[betalactamaseGenesEffective,]
  
  
  allSamplesWithCoreBetalactamases <- samplesWithBetalactamases
  allSamplesWithOnlyCoreBetalactamasesSuceptibleToClavulanicAcid <- samplesEffective
  
  AMP_HEADER <- "Studie.1.AMP"
  AMC_HEADER <- "KISS.I.AMC"
  
  # sirTable[,c("sample","AMP")]
  # sirTable[,"AMC"]
  # sirTable[,c("sample","AMC")]
  
  resultFrame <- data.frame(allSamples)
  colnames(resultFrame)[length(resultFrame)] <- "sampleID"
  
  resultFrame <- cbind(resultFrame,rep(NA,99))
  colnames(resultFrame)[length(resultFrame)] <- "AMP"
  resultFrame[,length(resultFrame)] <-  sirTable[,"AMP"]
  
  resultFrame <- cbind(resultFrame,rep(NA,99))
  colnames(resultFrame)[length(resultFrame)] <- "AMP-mm"
  resultFrame[,length(resultFrame)] <- millimeterTable[,AMP_HEADER]
  
  
  resultFrame <- cbind(resultFrame,rep(NA,99))
  colnames(resultFrame)[length(resultFrame)] <- "AMC"
  resultFrame[,length(resultFrame)] <-  sirTable[,"AMC"]
  
  resultFrame <- cbind(resultFrame,rep(NA,99))
  colnames(resultFrame)[length(resultFrame)] <- "AMC-mm"
  resultFrame[,length(resultFrame)] <- millimeterTable[,AMC_HEADER]
  
  aFrame <- readStatisticsFrame("samplemetrics-AMC-AMP-given-nonconformal-6",folder = getStatisticsSampleMetricsFolder(mode))
  if(!all(aFrame[order(aFrame$sample),'sample']==resultFrame$sampleID)){
    errorCondition("Samples not matching",!all(aFrame[order(aFrame$sample),'sample']==resultFrame$sampleID))
  }
  
  resultFrame <- cbind(resultFrame,rep(NA,99))
  colnames(resultFrame)[length(resultFrame)] <- "correct-AMC-AMP-given"
  resultFrame[,length(resultFrame)] = aFrame[order(aFrame$sample),]$correct
  resultFrame <- cbind(resultFrame,rep(NA,99))
  colnames(resultFrame)[length(resultFrame)] <- "ME-AMC-AMP-given"
  resultFrame[,length(resultFrame)] = aFrame[order(aFrame$sample),]$ME
  resultFrame <- cbind(resultFrame,rep(NA,99))
  colnames(resultFrame)[length(resultFrame)] <- "VME-AMC-AMP-given"
  resultFrame[,length(resultFrame)] = aFrame[order(aFrame$sample),]$VME
  
  aFrame <- readStatisticsFrame("samplemetrics-AMC-AMP-notgiven-nonconformal-6",folder = getStatisticsSampleMetricsFolder(mode))
  if(!all(aFrame[order(aFrame$sample),'sample']==resultFrame$sampleID)){
    errorCondition("Samples not matching",!all(aFrame[order(aFrame$sample),'sample']==resultFrame$sampleID))
  }
  
  resultFrame <- cbind(resultFrame,rep(NA,99))
  colnames(resultFrame)[length(resultFrame)] <- "correct-AMC-AMP-notgiven"
  resultFrame[,length(resultFrame)] = aFrame[order(aFrame$sample),]$correct
  resultFrame <- cbind(resultFrame,rep(NA,99))
  colnames(resultFrame)[length(resultFrame)] <- "ME-AMC-AMP-notgiven"
  resultFrame[,length(resultFrame)] = aFrame[order(aFrame$sample),]$ME
  resultFrame <- cbind(resultFrame,rep(NA,99))
  colnames(resultFrame)[length(resultFrame)] <- "VME-AMC-AMP-notgiven"
  resultFrame[,length(resultFrame)] = aFrame[order(aFrame$sample),]$VME
  
  
  resultFrame <- cbind(resultFrame,allSamples %in% allSamplesWithCoreBetalactamases)
  colnames(resultFrame)[length(resultFrame)] <- "allSamplesWithCoreBetalactamases"
  resultFrame <- cbind(resultFrame,allSamples %in% allSamplesWithOnlyCoreBetalactamasesSuceptibleToClavulanicAcid)
  colnames(resultFrame)[length(resultFrame)] <- "allSamplesWithOnlyCoreBetalactamasesSuceptibleToClavulanicAcid"
  
  
  #  sample <- resultFrame$sampleID[1]
  #  aGeneTable <- betalactamCoreGeneTable
  appendGene <- function(aGeneTable,resultFrame) {
    for(sample in resultFrame$sampleID){
      genesForSample <- aGeneTable$`Gene symbol`[aGeneTable[,colnames(aGeneTable)==sample]]
      value <- ""
      if(length(genesForSample)>0){
        value <- paste(genesForSample,collapse = ",")
      }
      resultFrame[resultFrame$sampleID==sample,length(resultFrame)] <- value
    }
    resultFrame
  }
  
  resultFrame <- cbind(resultFrame,rep(NA,99))
  colnames(resultFrame)[length(resultFrame)] <- "betalactamases-core"
  resultFrame <- appendGene(betalactamCoreGeneTable,resultFrame)
  
  # resultFrame <- cbind(resultFrame,rep(NA,99))
  # colnames(resultFrame)[length(resultFrame)] <- "betalactamases-plus"
  # resultFrame <- appendGene(betalactamPlusGeneTable,resultFrame)
  
  #  resultFrame <- cbind(resultFrame,rep(NA,99))
  #  colnames(resultFrame)[length(resultFrame)] <- "other-core"
  #  resultFrame <- appendGene(otherCoreGeneTable,resultFrame)
  
  #  resultFrame <- cbind(resultFrame,rep(NA,99))
  #  colnames(resultFrame)[length(resultFrame)] <- "other-plus"
  #  resultFrame <- appendGene(otherPlusGeneTable,resultFrame)
  
  
  betalactamGenesCoreFirst <- c(betalactamCoreGenes,betalactamPlusGenes)  
  #  gene <- betalactamGenesCoreFirst[1]
  for(gene in betalactamGenesCoreFirst)
  {
    samplesWithGene <- as.logical(geneTable[geneTable$`Gene symbol`==gene,][,-1])
    #    print(length(samplesWithGene))
    resultFrame <- cbind(resultFrame,samplesWithGene)
    colnames(resultFrame)[length(resultFrame)] <- gene
  }
  
  
  write.csv2(x= resultFrame,row.names = FALSE,file= paste(sep="",outdir, "\\" ,assemblymethod,"_clavulanic-acid_",mode,"_",amrfinderDatabase,".csv"))
  
  
  fetchGenes <- function(frame,listOfGenes)
  {
    result <- ""
    for(gene in listOfGenes) {
      print(gene)
      if(any(frame[,gene])){
        if(result==""){
          result = gene
        }
        else {
          result = paste(result,gene,sep=",")
        }
      }
    }
    result
  }
  
  phenotypeFrame <- NULL
  
  for ( a in c("S","R")){
    for ( b in c("S","R")){
      frame <- resultFrame[resultFrame$AMP==a & resultFrame$AMC==b,]
      genes <- fetchGenes(frame,betalactamCoreGenes)
      count <- dim(frame)[[1]]
      samples <- paste(frame$sampleID,sep=",",collapse=",")
      
      rbind(phenotypeFrame,c(
        a,
        b,
        count,
        samples,
        genes,
        mean(frame$`correct-AMC-AMP-given`),
        mean(frame$`ME-AMC-AMP-given`),
        mean(frame$`VME-AMC-AMP-given`),
        mean(frame$`correct-AMC-AMP-notgiven`),
        mean(frame$`ME-AMC-AMP-notgiven`),
        mean(frame$`VME-AMC-AMP-notgiven`)
      ))->phenotypeFrame      
    }
  }
  colnames(phenotypeFrame) <- c(
    "AMP",
    "AMC",
    "count",
    "samples",
    "genes",
    "correct-AMC-AMP-given",
    "ME-AMC-AMP-given",
    "VME-AMC-AMP-given",
    "correct-AMC-AMP-notgiven",
    "ME-AMC-AMP-notgiven",
    "VME-AMC-AMP-notgiven"
  )
  
  write.csv2(x= phenotypeFrame,row.names = FALSE,file= paste(sep="",outdir, "\\" ,assemblymethod,"_clavulanic-acid-summary_",mode,"_",amrfinderDatabase,".csv"))
  
  betalactamasesOutputTable <- caTable[,colnames(caTable) %in% c(GENE_SYMBOL,
                                                                 CLAVULANIC_ACID_EXPECTED_EFECTIVE,
                                                                 MOLECULAR_CLASS,
                                                                 FUNCTIONAL_GROUP)]
  colnames(betalactamasesOutputTable)
  aRow <- rep(NA,length(colnames(betalactamasesOutputTable)))
  names(aRow) <- colnames(betalactamasesOutputTable)
  #  AMP<-'R'
  #  AMC<-'S'
  for(AMP in c('S','R')) {
    for(AMC in c('S','R')) {
      columnName <- paste("AMP",AMP,"AMC",AMC,sep = "_")
      col <- rep("",nrow(betalactamasesOutputTable))
      names(col) <- betalactamasesOutputTable$Gene_symbol
      betalactamasesOutputTable$Gene_symbol
      #      gene <- "blaTEM"
      for(gene in betalactamasesOutputTable$Gene_symbol) {
        sampleAndPhenotypeTable <- resultFrame[resultFrame[,colnames(resultFrame)==gene]==TRUE,c("sampleID","AMP","AMC")]
        samples <- sampleAndPhenotypeTable[sampleAndPhenotypeTable$AMP==AMP & sampleAndPhenotypeTable$AMC==AMC,colnames(sampleAndPhenotypeTable)=="sampleID"]
        #length(samples)
        samplesString <- paste(samples,collapse=",")
        col[gene] <- paste(samples,collapse=",")
      }
      betalactamasesOutputTable <- cbind(betalactamasesOutputTable,col)
      colnames(betalactamasesOutputTable)[ncol(betalactamasesOutputTable)] <- columnName
    }
  }
  
  
  
  
  nobetalactamasesResultFrame <- resultFrame[resultFrame$`betalactamases-core`=="",]
  row <- aRow
  for(AMP in c('S','R')) {
    for(AMC in c('S','R')) {
      columnName <- paste("AMP",AMP,"AMC",AMC,sep = "_")
      sampleAndPhenotypeTable <- nobetalactamasesResultFrame[,c("sampleID","AMP","AMC")]
      samples <- sampleAndPhenotypeTable[sampleAndPhenotypeTable$AMP==AMP & sampleAndPhenotypeTable$AMC==AMC,colnames(sampleAndPhenotypeTable)=="sampleID"]
      #length(samples)
      samplesString <- paste(samples,collapse=",")
      row[columnName] <- samplesString
    }
  }
  
  betalactamasesOutputTable <- rbind(betalactamasesOutputTable,row)
  
  write.csv2(x= betalactamasesOutputTable,row.names = FALSE,file= paste(sep="",outdir, "\\" ,assemblymethod,"_betalactamases_",mode,"_",amrfinderDatabase,".csv"))
  
}

getBetalactamaseGeneForSample <-function(betalactamasesTable = readBetalactamasesTable())
{
  genes <- unique(betalactamasesTable$Gene_symbol)
  rv <- NULL
  for(geneSymbol in genes){
    #geneSymbol <- "blaDHA-1"
    subSet <- betalactamCoreGeneTable[betalactamCoreGeneTable$`Gene symbol` %in% geneSymbol,]
    hasGeneInGroup <- apply(subSet[-1], 2, any)
    if(is.null(rv)){
      rv <- data.frame(hasGeneInGroup)
      row.names(rv) <- colnames(subSet[-1])
    }else {
      rv <- cbind(rv,hasGeneInGroup) 
    }
    colnames(rv)[ncol(rv)] <- geneSymbol
  }
  rv[,order(colnames(rv))]
}


getBetalactamaseBushJacobiGroupForSample <-function(betalactamasesTable = readBetalactamasesTable())
{
  bushJacobi <- unique(betalactamasesTable$Functional_group)
  rv <- NULL
  for(group in bushJacobi){
    # group <- "2be"
    geneSymbol <- betalactamasesTable[betalactamasesTable$Functional_group==group,]$Gene_symbol
    subSet <- betalactamCoreGeneTable[betalactamCoreGeneTable$`Gene symbol` %in% geneSymbol,]
    hasGeneInGroup <- apply(subSet[-1], 2, any)
    if(is.null(rv)){
      rv <- data.frame(hasGeneInGroup)
      row.names(rv) <- colnames(subSet[-1])
    }else {
      rv <- cbind(rv,hasGeneInGroup) 
    }
    colnames(rv)[ncol(rv)] <- group
  }
  rv[,order(colnames(rv))]
}

getBetalactamaseAmblerGroupForSample <-function(betalactamasesTable = readBetalactamasesTable())
{
  caTable <-read.xlsx(paste(processedRootExcel,"Betalactamases.xlsx",sep="\\"))
  amblerGroup <- unique(betalactamasesTable$Molecular_class)
  rv <- NULL
  for(group in amblerGroup){
    # group <- "C"
    geneSymbol <- betalactamasesTable[betalactamasesTable$Molecular_class==group,]$Gene_symbol
    subSet <- betalactamCoreGeneTable[betalactamCoreGeneTable$`Gene symbol` %in% geneSymbol,]
    hasGeneInGroup <- apply(subSet[-1], 2, any)
    if(is.null(rv)){
      rv <- data.frame(hasGeneInGroup)
      row.names(rv) <- colnames(subSet[-1])
    }else {
      rv <- cbind(rv,hasGeneInGroup) 
    }
    colnames(rv)[ncol(rv)] <- group
  }
  rv[,order(colnames(rv))]
}

dumpClassperSample<- function()
{
  x <- getBetalactamaseBushJacobiGroupForSample()
  x <- cbind(rownames(x),x)
  colnames(x)[1] <- "sampleID"
  write.csv2(x,row.names = FALSE,file= paste(sep="",outdir, "\\" ,assemblymethod,"_bushJacobi_",amrfinderDatabase, ".csv"))
  
  x <- getBetalactamaseAmblerGroupForSample()
  x <- cbind(rownames(x),x)
  colnames(x)[1] <- "sampleID"
  write.csv2(x,row.names = FALSE,file= paste(sep="",outdir, "\\" ,assemblymethod,"_ambler_",amrfinderDatabase, ".csv"))
}

getSampleMetricsAMC <-function(mode)
{
  readStatisticsFrame("samplemetrics-AMC-nonconformal-6",folder = getStatisticsSampleMetricsFolder(mode))
}

ALL <- function()
{
  stopifnot(dir.exists(processedRootRassembly))
  
  
  if(!dir.exists(outdir)){
    dir.create(outdir)
  }
  analyzeClavulanicAcidNew()
  
  dumpClassperSample()
  outdir 
  
  mode = "Mode-A"
  analyzeClavulanicAcid(mode)
  
  mode = "Mode-B"
  analyzeClavulanicAcid(mode)
  
  mode = "Mode-C"
  analyzeClavulanicAcid(mode)
}