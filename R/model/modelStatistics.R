source(file = 'model/modelcommon.R')
library(stringr)




errorRatesNormalized <-function()
{
  phenotypeExcel =  paste(processedRootExcel, "ModelInput.xlsx", sep="/")
  modelLimitTable <- sheetToDataFrame(read_xlsx(phenotypeExcel,sheet = SHEET_MODEL_LIMIT))
  
  modes <- colnames(modelLimitTable[2:ncol(modelLimitTable)])
  mode <- "Mode A"
  for(mode in modes){
    modelFile <- paste("sirAntibioticsModel","_",sub(" ", "-", mode),".csv",sep="")
    sirAntibioticsModel <- read.csv2(row.names=1,paste(modelDirectory,"input",modelFile,sep="/"))
    demographicsModel <- read.csv2(row.names=1,paste(modelDirectory,"input","demographicsModel.csv",sep="/"))
    antibioticsNames <- colnames(sirAntibioticsModel)
    inputAntibioticsNames <- antibioticsNames
    numberOfInputAntibiotics <- length(inputAntibioticsNames)
    sirDataFrame <- sirAntibioticsModel
    demographicsDataframe <- demographicsModel
    
    sirDataFrameWords <- makeWordsDataFrame(sirDataFrame,demographicsModel)
    csvTable <- sampleAsColumns(sirDataFrameWords)
    file = paste(modelDirectory,"input",paste("sirAntibioticsModelWords_",sub(" ", "-", mode),".csv",sep=""),  sep="/")
#    write.csv2(x=csvTable,file=file,row.names = FALSE)
    
    demographicsModelNoDate <- demographicsModel
    demographicsModelNoDate$date <- "<unk>"
    sirDataFrameWords <- makeWordsDataFrame(sirDataFrame,demographicsModelNoDate)
    csvTable <- sampleAsColumns(sirDataFrameWords)
    file = paste(modelDirectory,"input",paste("sirAntibioticsModelWordsNoDate_",sub(" ", "-", mode),".csv",sep=""),  sep="/")
#    write.csv2(x=csvTable,file=file,row.names = FALSE)
  }
  
  
}

GENERATE_STATISTICS_ALT <- function()
{
  #range <- 1:4
  range <- 1:13
  allStats = NULL
  for(index in range){
    allstatsForIndex <- generateAllStatsAlternative(index)
    writeStatisticsFrame(allstatsForIndex, name = paste("allstats_alt",index,sep="_"),getStatisticsTmpFolder(MODE))
    if(is.null(allStats)){
      allStats <- allstatsForIndex
    } else {
      allStats <- rbind(allStats,allstatsForIndex)
    }
    writeStatisticsFrame(allStats, name = paste("allstats_alt_part",index,sep="_"),getStatisticsTmpFolder(MODE))
  }
  writeStatisticsFrame(allStats, name = "allstats_alt",getStatisticsFolder(MODE))
}

GENERATE_STATISTICS_ALT_FROM_WIDER <- function()
{
  altWider <-readStatisticsFrame("allstats_alt_wider",getStatisticsFolder(MODE))
  sum_wider <-  altWider %>% 
    dplyr::group_by(noinputab,metric,antibiotic,significanslevel) %>%
    dplyr::summarize(count=sum(count),.groups = "drop")
  
  sum_wider <- sum_wider %>% dplyr::select(noinputab, metric, count, antibiotic, significanslevel)
  sum_wider <- sum_wider  %>% dplyr::arrange_all()
  
  # alt <-readStatisticsFrame("allstats_alt",getStatisticsFolder(MODE))
  # alt <- alt %>% dplyr::arrange_all()
  # 
  # all(alt==sum_wider,na.rm=TRUE)
  # any(alt!=sum_wider,na.rm=TRUE)
  
  writeStatisticsFrame(sum_wider, name = "allstats_alt_from_wider",getStatisticsFolder(MODE))
}


GENERATE_STATISTICS_ALT_WIDER <- function()
{
  #range <- 1:4
  range <- 1:13
  allStats = NULL
  for(index in range){
    allstatsForIndex <- generateAllStatsAlternativeWider(index)
    writeStatisticsFrame(allstatsForIndex, name = paste("allstats_alt_wider",index,sep="_"),getStatisticsTmpFolder(MODE))
    if(is.null(allStats)){
      allStats <- allstatsForIndex
    } else {
      allStats <- rbind(allStats,allstatsForIndex)
    }
    writeStatisticsFrame(allStats, name = paste("allstats_alt_wider_part",index,sep="_"),getStatisticsTmpFolder(MODE))
  }
  writeStatisticsFrame(allStats, name = "allstats_alt_wider",getStatisticsFolder(MODE))
}

GENERATE_STATISTICS_NON_ATU <- function()
{
  isolatesNonATU <- as.character(c( 2, 3, 4, 9, 11, 12, 18, 19, 20, 24, 25, 26, 28, 29, 32, 34, 35, 36, 37, 39, 41, 44, 45, 46, 48, 49, 51, 52, 53, 56, 57, 58, 60, 61, 65, 67, 69, 72, 75, 76, 78, 79, 82, 84, 86, 87, 88, 95, 96, 98, 101, 102, 105, 110, 112, 113, 116, 120, 121, 122, 123, 124, 125 ))
  range <- 1:13
  #range <- 1:2
  allStats = NULL
  for(index in range){
    allstatsForIndex <- generateAllStatsAlternative(index,isolatesNonATU)
    writeStatisticsFrame(allstatsForIndex, name = paste("nonATUstats",index,sep="_"),getStatisticsTmpFolder(MODE))
    if(is.null(allStats)){
      allStats <- allstatsForIndex
    } else {
      allStats <- rbind(allStats,allstatsForIndex)
    }
    writeStatisticsFrame(allStats, name = paste("nonATUstats_part",index,sep="_"),getStatisticsTmpFolder(MODE))
  }
  writeStatisticsFrame(allStats, name = paste("nonATUstats"),getStatisticsSubsetFolder(MODE))
}


GENERATE_METRICS_ALT <- function(name = "allStats_alt",prefix=NA,folder = getStatisticsFolder(MODE))
{

  errorStatisticsFull <- errorStatistics(name,folder)
  errorStatistics <- errorStatisticsFull %>% dplyr::select(noinputab,significanslevel, correct=correct_frac,VME=VME_frac,ME=ME_frac)
  errorStatisticsNormalized <- errorStatisticsFull %>% dplyr::select(noinputab,significanslevel,VME=VME_frac_norm,ME=ME_frac_norm)
  namebase <- "errorStatistics"
  if(!is.na(prefix)){
    namebase <- paste (prefix,namebase,sep="-")
  }
  writeStatisticsFrame(errorStatisticsFull, name = paste(namebase, "Full",sep=""),folder)
  writeStatisticsFrame(errorStatistics, name = namebase,folder)
  writeStatisticsFrame(errorStatisticsNormalized, name = paste(namebase,"Normalized",sep=""),folder)

  abErrorStatisticsFull <- errorStatisticsPerAntibiotic(name,folder)
  abErrorStatistics <- abErrorStatisticsFull %>% dplyr::select(noinputab,significanslevel,antibiotic,correct=correct_frac,VME=VME_frac,ME=ME_frac)
  abErrorStatisticsNormalized <- abErrorStatisticsFull %>% dplyr::select(noinputab,significanslevel,antibiotic,VME=VME_frac_norm,ME=ME_frac_norm)
  namebase <- "abErrorStatistics"
  if(!is.na(prefix)){
    namebase <- paste (prefix,namebase,sep="-")
  }
  writeStatisticsFrame(abErrorStatisticsFull,  name = paste(namebase, "Full",sep=""),folder)
  writeStatisticsFrame(abErrorStatistics, name = namebase,folder)
  writeStatisticsFrame(abErrorStatisticsNormalized, name = paste(namebase,"Normalized",sep=""),folder)
}




generateAllStatsAlternativeWider <- function(index,isolateSubset=NA)
{
  #isolateSubset <- c("2","125")
  #index <- 2
  #inputAntibiotics <- c("AMP","AMC","CRO","CIP","OFX","TOB")
  
  allAbs <- c("AMP","AMC","PIP","TZP","CAZ","CRO","CTX","FEP","CIP","OFX","LVX","MFX","GEN","TOB")
  correctFrame <- readOutputFrameShort(shortfile="modelOutput_answer_13.csv",folder=getPredictionFolder(MODE))
  if(length(isolateSubset)==1 && is.na(isolateSubset)){
    isolateSubset <- rownames(correctFrame)
  }
  correctFrame <- correctFrame[isolateSubset,]
  
  predsFiles <- makeShortFileName("preds",index)
  for(signi in signis){
    predsFiles <- c(predsFiles,makeShortFileName("confpreds",index,signi))
  }
  names(predsFiles) <- c("crude",signis)
  
  frames <- lapply(predsFiles,function(x){ 
    aFrame <- readOutputFrameShort(shortfile=x,folder=getPredictionFolder(MODE))[isolateSubset,]
  })
  
  returnCols <- c("noinputab","metric","count", "antibiotic","significanslevel","word")
  returnFrame <- data.frame(matrix(nrow=0,ncol=length(returnCols)))
  colnames(returnFrame) <- returnCols

  for(ab in allAbs){
    #ab <- "TZP"
    abR <- paste(ab,"R",sep="_")
    abS <- paste(ab,"S",sep="_")
    abSR <- paste(ab,"SR",sep="_")
    columnWithoutAb <- paste(allAbs[allAbs!=ab],collapse = "_")
    correctList=correctFrame[,columnWithoutAb]
    
    for(name in names(frames)){
      #name = "crude"
      frame <- frames[[name]]
      
      subframeShouldBeS <- frame[correctList==abS,]
      subframeShouldBeR <- frame[correctList==abR,]
      
      sumMatch <-function(subframe,pattern){
        enLista <- apply(
          apply(subframe,c(1,2),
                function(x) {
                  grepl(pattern,x)})
          ,2,sum)
        
      }
      
      correctS <- sumMatch(subframeShouldBeS,abS) - sumMatch(subframeShouldBeS,abSR)
      ME <- sumMatch(subframeShouldBeS,abR)
      riskForME <-unlist(lapply(colnames(subframeShouldBeS),function(x){if(grepl(ab,x)==FALSE){nrow(subframeShouldBeS)}else{0}}))
      names(riskForME) <- colnames(subframeShouldBeS)
      #(correctS + ME)==riskForME
      #ME/riskForME
      
      correctR <-sumMatch(subframeShouldBeR,abR)
      VME <-  sumMatch(subframeShouldBeR,abS) -sumMatch(subframeShouldBeR,abSR)
      riskForVME <- unlist(lapply(colnames(subframeShouldBeR),function(x){if(grepl(ab,x)==FALSE){nrow(subframeShouldBeR)}else{0}}))
      names(riskForVME) <- colnames(subframeShouldBeR) 

      #min(VME/riskForVME,na.rm=TRUE)
      
      
      correct <- correctS + correctR
      
      predicted <- correct + VME + ME
      ambiguous <-  sumMatch(frame,abSR)
      predorambiguous <- predicted + ambiguous
      total <- riskForME + riskForVME
      notpredicted <- total - predorambiguous
      #total <- sum(grepl(ab,colnames(frame))==FALSE)*nrow(frame)
      
      #cols  "noinputab"        "metric"           "count"            "antibiotic"       "significanslevel"
      #metric "correct"         "VME"             "ME"              "predicted"       "ambiguous"       "predorambiguous" "notpredicted"       
      #significance level      NA      "0.010" "0.025" "0.050" "0.100" 
      significansLevel <- NA
      if(name!="crude"){
        significansLevel <- fixSigni(name)
      }
      
      #returnCols <- c("noinputab","metric","count", "antibiotic","significanslevel","word")
      tempFrame <- data.frame(matrix(nrow=0,ncol=length(returnCols)))

      
      #word <- "AMP_AMC"
      for(word in names(correct)){
        tempFrame <- rbind(tempFrame,
              c(index,"correct",correct[word],ab,significansLevel,word),
              c(index,"correctS",correctS[word],ab,significansLevel,word),
              c(index,"correctR",correctR[word],ab,significansLevel,word),
              c(index,"VME",VME[word],ab,significansLevel,word),
              c(index,"riskForVME",riskForVME[word],ab,significansLevel,word),
              c(index,"ME",ME[word],ab,significansLevel,word),
              c(index,"riskForME",riskForME[word],ab,significansLevel,word),
              c(index,"predicted",predicted[word],ab,significansLevel,word),
              c(index,"ambiguous",ambiguous[word],ab,significansLevel,word),
              c(index,"predorambiguous",predorambiguous[word],ab,significansLevel,word),
              c(index,"notpredicted",notpredicted[word],ab,significansLevel,word),
              c(index,"total",total[word],ab,significansLevel,word))
      }
      colnames(tempFrame) <- c("noinputab","metric","count", "antibiotic","significanslevel","word")
      colnames(returnFrame) <- returnCols

      returnFrame <- rbind(returnFrame, tempFrame)
    }
  }
  colnames(returnFrame) <- returnCols
  returnFrame$count <- as.integer(returnFrame$count)
  returnFrame  
}


generateAllStatsAlternative <- function(index,isolateSubset=NA)
{
  #isolateSubset <- c("2","125")
  #index <- 2
  #inputAntibiotics <- c("AMP","AMC","CRO","CIP","OFX","TOB")

  allAbs <- c("AMP","AMC","PIP","TZP","CAZ","CRO","CTX","FEP","CIP","OFX","LVX","MFX","GEN","TOB")
  correctFrame <- readOutputFrameShort(shortfile="modelOutput_answer_13.csv",folder=getPredictionFolder(MODE))
  if(length(isolateSubset)==1 && is.na(isolateSubset)){
    isolateSubset <- rownames(correctFrame)
  }
  correctFrame <- correctFrame[isolateSubset,]
  
  predsFiles <- makeShortFileName("preds",index)
  for(signi in signis){
    predsFiles <- c(predsFiles,makeShortFileName("confpreds",index,signi))
  }
  names(predsFiles) <- c("crude",signis)
  
  frames <- lapply(predsFiles,function(x){ 
    aFrame <- readOutputFrameShort(shortfile=x,folder=getPredictionFolder(MODE))[isolateSubset,]
    })
  
  returnCols <- c("noinputab","metric","count", "antibiotic","significanslevel")
  returnFrame <- data.frame(matrix(nrow=0,ncol=length(returnCols)))
  colnames(returnFrame) <- returnCols
  
  for(ab in allAbs){
    #ab <- "TZP"
    abR <- paste(ab,"R",sep="_")
    abS <- paste(ab,"S",sep="_")
    abSR <- paste(ab,"SR",sep="_")
    columnWithoutAb <- paste(allAbs[allAbs!=ab],collapse = "_")
    correctList=correctFrame[,columnWithoutAb]
    
    for(name in names(frames)){
      #name = "crude"
      frame <- frames[[name]]
      
      subframeShouldBeS <- frame[correctList==abS,]
      subframeShouldBeR <- frame[correctList==abR,]
      
      correctS <- sum(apply(subframeShouldBeS,c(1,2),function(x) {grepl(abS,x)})) - sum(apply(subframeShouldBeS,c(1,2),function(x) {grepl(abSR,x)}))
      ME <-  sum(apply(subframeShouldBeS,c(1,2),function(x) {grepl(abR,x)}))
      riskForME <- sum(grepl(ab,colnames(subframeShouldBeS))==FALSE)*nrow(subframeShouldBeS)
      #ME/riskForME
      
      correctR <- sum(apply(subframeShouldBeR,c(1,2),function(x) {grepl(abR,x)}))
      VME <-  sum(apply(subframeShouldBeR,c(1,2),function(x) {grepl(abS,x)})) - sum(apply(subframeShouldBeR,c(1,2),function(x) {grepl(abSR,x)}))
      riskForVME <- sum(grepl(ab,colnames(subframeShouldBeR))==FALSE)*nrow(subframeShouldBeR)

      #VME/riskForVME
      
            
      correct <- correctS + correctR

      predicted <- correct + VME + ME
      ambiguous <-  sum(apply(frame,c(1,2),function(x) {grepl(abSR,x)}))
      predorambiguous <- predicted + ambiguous
      total <- riskForME + riskForVME
      notpredicted <- total - predorambiguous
      #total <- sum(grepl(ab,colnames(frame))==FALSE)*nrow(frame)
      
      #cols  "noinputab"        "metric"           "count"            "antibiotic"       "significanslevel"
      #metric "correct"         "VME"             "ME"              "predicted"       "ambiguous"       "predorambiguous" "notpredicted"       
      #significance level      NA      "0.010" "0.025" "0.050" "0.100" 
      significansLevel <- NA
      if(name!="crude"){
        significansLevel <- fixSigni(name)
      }
      returnFrame <- rbind(returnFrame,
                           c(index,"correct",correct,ab,significansLevel),
                           c(index,"correctS",correctS,ab,significansLevel),
                           c(index,"correctR",correctR,ab,significansLevel),
                           c(index,"VME",VME,ab,significansLevel),
                           c(index,"riskForVME",riskForVME,ab,significansLevel),
                           c(index,"ME",ME,ab,significansLevel),
                           c(index,"riskForME",riskForME,ab,significansLevel),
                           c(index,"predicted",predicted,ab,significansLevel),
                           c(index,"ambiguous",ambiguous,ab,significansLevel),
                           c(index,"predorambiguous",predorambiguous,ab,significansLevel),
                           c(index,"notpredicted",notpredicted,ab,significansLevel),
                           c(index,"total",total,ab,significansLevel)
      )
    }
  }
  colnames(returnFrame) <- returnCols
  returnFrame$count <- as.integer(returnFrame$count)
  returnFrame  
}


compareAllstats <- function()
{
  # allStatsOrg <-readStatisticsFrame("allstats",getStatisticsFolder(MODE))
  # allStatsNew <-readStatisticsFrame("allstats_alt",getStatisticsFolder(MODE))
  # 
  # allStatsOrgArranged <- allStatsOrg %>% dplyr::arrange_all()
  # allStatsNewArranged <- allStatsNew %>% dplyr::arrange_all()
  # 
  # diff <- allStatsOrgArranged!=allStatsNewArranged
  # same <- allStatsOrgArranged==allStatsNewArranged
  # any(diff,na.rm = TRUE)
  # all(same,na.rm=TRUE)
  

  
  #  allStatsAltWider <-readStatisticsFrame("allstats_alt",getStatisticsFolder(MODE))
  # sum_wider <-  allStatsAltWider %>% 
  #   dplyr::group_by(noinputab,metric,antibiotic,significanslevel) %>%
  #   dplyr::summarize(count=sum(count),.groups = "drop")
  # 
  # sum_wider <- sum_wider %>% select(noinputab, metric, count, antibiotic, significanslevel)
  # sum_wider <- sum_wider  %>% dplyr::arrange_all()
  
    
  allStatsAlt <-readStatisticsFrame("allstats_alt",getStatisticsFolder(MODE))
  allStatsAltFromWider <-readStatisticsFrame("allstats_alt_from_wider",getStatisticsFolder(MODE))
  
  allStatsAltSorted <- allStatsAlt %>% dplyr::arrange_all()
  allStatsAltFromWiderSorted <- allStatsAltFromWider %>% dplyr::arrange_all()

  head(allStatsAltSorted)
  head(allStatsAltFromWiderSorted)
  tail(allStatsAltSorted)
  tail(allStatsAltFromWiderSorted)
  
  allSame <- all(allStatsAltSorted==allStatsAltFromWiderSorted,na.rm=TRUE)
  noneDifferent <- any(allStatsAltSorted!=allStatsAltFromWiderSorted,na.rm=TRUE)
  ok <- allSame && (!noneDifferent)
  if(ok)
    print("Ok")
  else
    print("Not ok")
  ok
}








abSusceptableFraction <- function()
{
  allAbs <- c("AMP","AMC","PIP","TZP","CAZ","CRO","CTX","FEP","CIP","OFX","LVX","MFX","GEN","TOB")
  correctFrame <- readOutputFrameShort(shortfile="modelOutput_answer_13.csv",folder=getPredictionFolder(MODE))
  rv <- unlist(lapply(allAbs,function(ab){
    #ab <- "CTX"
    columnWithoutAb <- paste(allAbs[allAbs!=ab],collapse = "_")
    correctList=correctFrame[,columnWithoutAb]
    abR <- paste(ab,"R",sep="_")
    abS <- paste(ab,"S",sep="_")
    shouldBeS <- correctList==abS
    shouldBeR <- correctList==abR
    fractionS <- sum(shouldBeS)/sum(shouldBeS | shouldBeR )
    fractionS
  }))
  names(rv) <- allAbs
  rv
}

susceptableFraction <- function()
{
  correctFrame <- readOutputFrameShort(shortfile="modelOutput_answer_13.csv",folder=getPredictionFolder(MODE))
  numberOfR <- sum(apply(correctFrame,c(1,2),function(x)grepl(x=x,pattern="_R")))
  numberOfS <- sum(apply(correctFrame,c(1,2),function(x)grepl(x=x,pattern="_S")))
  rv = numberOfS/(numberOfR+numberOfS)
  rv
}


errorStatisticsPerAntibiotic <- function(name="allstats_alt",folder = getStatisticsFolder(MODE))
{
  #name <- "statsAMC_CTX_OFX_GEN"
  allAbs <- c("AMP","AMC","PIP","TZP","CAZ","CRO","CTX","FEP","CIP","OFX","LVX","MFX","GEN","TOB")
  selectedMetrics <- c("correct","ME","VME")
  allMetrics <- c("correct","ME","VME","ambiguous","notpredicted")

  allStats <- readStatisticsFrame(name,folder)
  susceptibleFractionAntibiotics <- abSusceptableFraction()
  #unique(allStats$antibiotic)
  
  rv <- NULL
      
  for(ab in unique(allStats$antibiotic)){
    #ab <- "AMC"
    abStats <- allStats %>% dplyr::filter(antibiotic==ab)
    #abStats %>% dplyr::filter(noinputab==4)

    allMetricsCount <- abStats %>% dplyr::group_by(noinputab,significanslevel) %>% tidyr::pivot_wider(names_from=metric,values_from = count)
    #allMetricsCount %>% dplyr::filter(noinputab==4)
    summary_mutated <- allMetricsCount %>% 
      dplyr::mutate(
        correct_frac=correct/total,
        VME_frac=VME/total,
        ME_frac=ME/total,
        VME_frac_norm=VME/riskForVME,
        ME_frac_norm=ME/riskForME
        )
    #summary_mutated %>% dplyr::filter(noinputab==4) %>% dplyr::select(VME_frac,VME_frac_norm,VME_frac_norm_alt,ME_frac,ME_frac_norm,ME_frac_norm_alt)
    if(is.null(rv)){
      rv <- summary_mutated
    } else {
      rv <- rbind(rv,summary_mutated)
    }
  }
  #rv %>% dplyr::filter(noinputab==4,is.na(significanslevel))
  rv
}


errorStatistics <- function(name="allstats_alt",folder = getStatisticsFolder(MODE))
{
  selectedMetrics <- c("correct","ME","VME")
  allMetrics <- c("correct","ME","VME","ambiguous","notpredicted")
  
  allStats <- readStatisticsFrame(name,folder)


  allMetricsCount <- allStats %>% dplyr::group_by(noinputab,significanslevel,metric) %>% dplyr::summarize(count=sum(count), .groups = "drop") %>% tidyr::pivot_wider(names_from=metric,values_from = count)

  summary_mutated <- allMetricsCount %>% 
      dplyr::mutate(
        correct_frac=correct/total,
        VME_frac=VME/total,
        ME_frac=ME/total,
        VME_frac_norm=VME/riskForVME,
        ME_frac_norm=ME/riskForME
      )
  #summary_mutated %>%  dplyr::select(noinputab,significanslevel,correct,VME,ME,VME_frac_norm_alt,ME_frac_norm_alt)
    
    
  summary_mutated
}




withPlyr <-function()
{
#  dplyr

  pAllMetrics <- c("correct","ME","VME","ambiguous","notpredicted")
  pInputabs= 6

  significansLevels <- c(NA,"0.100","0.050","0,025")

  pAntibiotic = "CTX"
  VMES <- allStats %>% dplyr::filter(noinputab==pInputabs,metric=="ME")
  MES <- allStats %>% dplyr::filter(noinputab==pInputabs,metric=="ME")
  ALL <- allStats %>% dplyr::filter(noinputab==pInputabs,metric %in% pAllMetrics )

  MES <- MES %>% dplyr::filter(antibiotic==pAntibiotic)
  ALL <- ALL %>% dplyr::filter(antibiotic==pAntibiotic)

  ALLSummary <- ALL %>% dplyr::group_by(significanslevel) %>%
    dplyr::summarize(count_all_metrics = sum(count))
  MESummary <- MES %>% dplyr::group_by(significanslevel) %>%
    dplyr::summarize(count_all_metrics = sum(count))
  MESummary

  percentageByLevel <- function(aTibble,allTibble) {
    #aTibble <- MES
    nominator <- aTibble %>% dplyr::group_by(significanslevel) %>%
      dplyr::summarize(count_all_metrics = sum(count))
    denominator <- allTibble %>% dplyr::group_by(significanslevel) %>%
      dplyr::summarize(count_all_metrics = sum(count))
    rv <- cbind(nominator[,1],nominator[,2]/denominator[,2])
    rv
  }

  ME <- percentageByLevel(MES,ALL)
  VME <- percentageByLevel(VMES,ALL)

  MES %>%  dplyr::group_by(significanslevel) %>%
    dplyr::summarize(count_all_metrics = sum(count))

}

GENERATE_BEST_AND_ALT_METRICS <-function()
{
  altWider <-readStatisticsFrame("allstats_alt_wider",getStatisticsFolder(MODE))
  #tail(altWider)
  BEST_ALT <- list(
    c("AMC","CTX","OFX","GEN"),
    c("PIP","TZP","CTX","LVX","MFX","GEN"),
    c("AMC" ,"TZP" ,"CAZ", "CTX", "FEP", "OFX","LVX" ,"GEN"),
    c("AMP", "PIP" ,"CRO","CTX","FEP", "CIP", "MFX" ,"TOB"))
  
  for(bestWord in append(BEST_ABS,BEST_ALT)){
    #bestWord <- c("AMP", "AMC", "CRO", "CIP", "OFX", "TOB")
    bestWordCollapsed <- paste(bestWord,collapse="_")
    altWiderBestWord <- altWider %>% dplyr::filter(word==bestWordCollapsed) %>% dplyr::select(-word)
    altWiderBestWord <- altWiderBestWord %>% dplyr::filter(!antibiotic %in% bestWord)
    writeStatisticsFrame(altWiderBestWord,paste("stats",bestWordCollapsed,sep=""),getStatisticsSubsetFolder(MODE))
    GENERATE_METRICS_ALT(name=paste("stats",bestWordCollapsed,sep=""),prefix=bestWordCollapsed,getStatisticsSubsetFolder(MODE))
  }
  
}

checkDirs <-function()
{
  stopifnot(dir.exists(getPredictionBase(MODE)))
  
  
  if(!dir.exists(getStatisticsFolder(MODE))){
    dir.create(getStatisticsFolder(MODE))
  }
  if(!dir.exists(getStatisticsSubsetFolder(MODE))){
    dir.create(getStatisticsSubsetFolder(MODE))
  }
  
  if(!dir.exists(getStatisticsTmpFolder(MODE))){
    dir.create(getStatisticsTmpFolder(MODE))
  }
  
}

MOST <-  function()
{
  checkDirs()

  GENERATE_STATISTICS_ALT_FROM_WIDER()
  GENERATE_BEST_AND_ALT_METRICS()
  GENERATE_METRICS_ALT()
  GENERATE_METRICS_ALT(name="nonATUstats",prefix="nonATUstats",getStatisticsSubsetFolder(MODE))
  GENERATE_RISK_STRATIFIED_STATISTICS()
}

ALL <- function()
{
  checkDirs()
  GENERATE_STATISTICS_ALT()
  GENERATE_STATISTICS_ALT_WIDER()
  GENERATE_STATISTICS_NON_ATU()
  MOST()
}


GENERATE_RISK_STRATIFIED_STATISTICS<- function()
{
  riskForMEStatistics()
  riskForVMEStatistics()
  riskForMEStatisticsAb()
  riskForVMEStatisticsAb()
}

riskForMEStatistics <- function()
{
  allStats <-readStatisticsFrame("allstats_alt",getStatisticsFolder(MODE))  
  subset <- allStats %>% 
    dplyr::group_by(noinputab,significanslevel,metric) %>% dplyr::summarize(count=sum(count), .groups = "drop")  %>%
    dplyr::filter(metric %in% c("riskForME","correctS","ME")) %>%
    tidyr::pivot_wider(names_from=metric,values_from = count)
  fractions <- subset %>% dplyr::mutate(MEfrac=ME/riskForME,correctfrac=correctS/riskForME)  %>% dplyr::select (noinputab,significanslevel,correct=correctfrac,ME=MEfrac)
  writeStatisticsFrame(fractions,"riskForMEStatistics",folder = getStatisticsFolder(MODE))  
}

riskForVMEStatistics <- function()
{
  allStats <-readStatisticsFrame("allstats_alt",getStatisticsFolder(MODE))  
  subset <- allStats %>% 
    dplyr::group_by(noinputab,significanslevel,metric) %>% dplyr::summarize(count=sum(count), .groups = "drop")  %>%
    dplyr::filter(metric %in% c("riskForVME","correctR","VME")) %>%
    tidyr::pivot_wider(names_from=metric,values_from = count)
  fractions <- subset %>% dplyr::mutate(VMEfrac=VME/riskForVME,correctfrac=correctR/riskForVME)  %>% dplyr::select (noinputab,significanslevel,correct=correctfrac,VME=VMEfrac)
  writeStatisticsFrame(fractions,"riskForVMEStatistics",folder = getStatisticsFolder(MODE))  
}


riskForMEStatisticsAb <- function()
{
  allStats <-readStatisticsFrame("allstats_alt",getStatisticsFolder(MODE))  
  subset <- allStats %>% 
        dplyr::filter(metric %in% c("riskForME","correctS","ME")) %>%
        tidyr::pivot_wider(names_from=metric,values_from = count)
  fractions <- subset %>% dplyr::mutate(MEfrac=ME/riskForME,correctfrac=correctS/riskForME)  %>% dplyr::select (noinputab,antibiotic , significanslevel,correct=correctfrac,ME=MEfrac)
  writeStatisticsFrame(fractions,"riskForMEStatisticsAb",folder = getStatisticsFolder(MODE))  
}



riskForVMEStatisticsAb <- function()
{
  allStats <-readStatisticsFrame("allstats_alt",getStatisticsFolder(MODE))  
  subset <- allStats %>% 
    dplyr::filter(metric %in% c("riskForVME","correctR","VME")) %>%
    tidyr::pivot_wider(names_from=metric,values_from = count)
  fractions <- subset %>% dplyr::mutate(VMEfrac=VME/riskForVME,correctfrac=correctR/riskForVME)  %>% dplyr::select (noinputab,antibiotic , significanslevel,correct=correctfrac,VME=VMEfrac)
  writeStatisticsFrame(fractions,"riskForVMEStatisticsAb",folder = getStatisticsFolder(MODE))  
}


