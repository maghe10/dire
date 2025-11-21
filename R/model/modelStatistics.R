source(file = 'model/modelcommon.R')
library(stringr)

RANGE <- 1:13
#RANGE <- 1:2


reorderAntibioticsRows <-function(tibble)
{
  tibble %>% dplyr::arrange(match(antibiotic,ALL_ANTIBIOTICS_IN_MODEL))
}

createInputTables <-function()
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

#@deprecated
GENERATE_STATISTICS_ALT <- function()
{
  allStats = NULL
  for(index in RANGE){
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

#@deprecated
GENERATE_STATISTICS_ALT_FROM_WIDER <- function()
{
  altWider <-readStatisticsFrame("allstats_alt_wider",getStatisticsFolder(MODE))
  sum_wider <-  altWider %>% 
    dplyr::group_by(noinputab,metric,antibiotic,significanceLevel) %>%
    dplyr::summarize(count=sum(count),.groups = "drop")
  
  sum_wider <- sum_wider %>% dplyr::select(noinputab, metric, count, antibiotic, significanceLevel)
  sum_wider <- sum_wider  %>% dplyr::arrange_all()
  
  # alt <-readStatisticsFrame("allStats_sum",getStatisticsFolder(MODE))
  # alt <- alt %>% dplyr::arrange_all()
  # 
  # all(alt==sum_wider,na.rm=TRUE)
  # any(alt!=sum_wider,na.rm=TRUE)
  
  writeStatisticsFrame(sum_wider, name = "allstats_alt_from_wider",getStatisticsFolder(MODE))
}


#@deprecated
GENERATE_STATISTICS_ALT_WIDER <- function()
{
  allStats = NULL
  for(index in RANGE){
    allstatsForIndex <- generateAllStatsAlternativeWider(index)
#    writeStatisticsFrame(allstatsForIndex, name = paste("allstats_alt_wider",index,sep="_"),getStatisticsTmpFolder(MODE))
    if(is.null(allStats)){
      allStats <- allstatsForIndex
    } else {
      allStats <- rbind(allStats,allstatsForIndex)
    }
#    writeStatisticsFrame(allStats, name = paste("allstats_alt_wider_part",index,sep="_"),getStatisticsTmpFolder(MODE))
  }
  writeStatisticsFrame(allStats, name = "allstats_alt_wider",getStatisticsFolder(MODE))
}


#@deprecated
GENERATE_STATISTICS_ALT_WIDEST <- function()
{
  allStats = NULL
  for(index in RANGE){
    print(index)
    allstatsForIndex <- generateAllStatsAlternativeWidest(index)
    #    writeStatisticsFrame(allstatsForIndex, name = paste("allstats_alt_widest",index,sep="_"),getStatisticsTmpFolder(MODE))
    if(is.null(allStats)){
      allStats <- allstatsForIndex
    } else {
      allStats <- rbind(allStats,allstatsForIndex)
    }
    #    writeStatisticsFrame(allStats, name = paste("allstats_alt_widerst_part",index,sep="_"),getStatisticsTmpFolder(MODE))
  }
  writeStatisticsFrame(allStats, name = "allstats_alt_widest",getStatisticsFolder(MODE))
}




GENERATE_STATISTICS_NON_ATU <- function()
{
  isolatesNonATU <- as.character(c( 2, 3, 4, 9, 11, 12, 18, 19, 20, 24, 25, 26, 28, 29, 32, 34, 35, 36, 37, 39, 41, 44, 45, 46, 48, 49, 51, 52, 53, 56, 57, 58, 60, 61, 65, 67, 69, 72, 75, 76, 78, 79, 82, 84, 86, 87, 88, 95, 96, 98, 101, 102, 105, 110, 112, 113, 116, 120, 121, 122, 123, 124, 125 ))
  allStats = NULL
  for(index in RANGE){
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


GENERATE_SAMPLE_METRICS <- function(name = "allStats_sample",prefix=NA,folder = getStatisticsFolder(MODE))
{
  abSampleErrorStatisticsAll <- errorStatisticsPerAntibiotic(name,folder) %>% reorderAntibioticsRows()
  abSampleErrorStatisticsCount <- abSampleErrorStatisticsAll %>% dplyr::select(-c(correct_frac,VME_frac,ME_frac))
  abSampleErrorStatisticsFrac <- abSampleErrorStatisticsAll %>% dplyr::select(noinputab,significanceLevel,antibiotic,sample,correct=correct_frac,VME=VME_frac,ME=ME_frac)
  
  namebase <- "abSampleErrorStatistics"
  if(!is.na(prefix)){
    namebase <- paste (prefix,namebase,sep="-")
  }
  writeStatisticsFrame(abSampleErrorStatisticsCount,  name = paste(namebase, "Count",sep=""),folder)
  writeStatisticsFrame(abSampleErrorStatisticsFrac, name = paste(namebase, "Frac",sep=""),folder)
  
}
  
GENERATE_METRICS_ALT <- function(name = "allStats_sum",prefix=NA,folder = getStatisticsFolder(MODE))
{

  errorStatisticsAll <- errorStatistics(name,folder)
  errorStatisticsFrac <- errorStatisticsAll %>% dplyr::select(noinputab,significanceLevel, correct=correct_frac,VME=VME_frac,ME=ME_frac)
  errorStatisticsCount <- errorStatisticsAll %>% dplyr::select(-c(correct_frac,VME_frac,ME_frac))
  
  namebase <- "errorStatistics"
  if(!is.na(prefix)){
    namebase <- paste (prefix,namebase,sep="-")
  }
  writeStatisticsFrame(errorStatisticsCount, name = paste(namebase, "Count",sep=""),folder)
  writeStatisticsFrame(errorStatisticsFrac, name = paste(namebase, "Frac",sep=""),folder)

  abErrorStatisticsAll <- errorStatisticsPerAntibiotic(name,folder) %>% reorderAntibioticsRows()
  abErrorStatisticsFrac <- abErrorStatisticsAll %>% dplyr::select(noinputab,significanceLevel,antibiotic,correct=correct_frac,VME=VME_frac,ME=ME_frac)
  abErrorStatisticsCount <- abErrorStatisticsAll %>% dplyr::select(-c(correct_frac,VME_frac,ME_frac))
  
  namebase <- "abErrorStatistics"
  if(!is.na(prefix)){
    namebase <- paste (prefix,namebase,sep="-")
  }
  writeStatisticsFrame(abErrorStatisticsCount,  name = paste(namebase, "Count",sep=""),folder)
  writeStatisticsFrame(abErrorStatisticsFrac, name = paste(namebase, "Frac",sep=""),folder)
}




generateAllStatsAlternativeWider <- function(index,isolateSubset=NA,columns = NA)
{
  #isolateSubset <- c("2","125")
  #index <- 1
  #inputAntibiotics <- c("AMP","AMC","CRO","CIP","OFX","TOB")
  #
  #isolateSubset=NA
  #index <- 4
  #columns <- c("AMP_CAZ_CIP_GEN" ,"AMP_CAZ_CIP_TOB" ,"AMP_CAZ_OFX_GEN")

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
    if(length(columns)>1){
      aFrame <- aFrame[,columns]
    }
    aFrame
  })
  
  returnCols <- c("noinputab","metric","count", "antibiotic","significanceLevel","word")
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
      
      # among should be S
      correctS <- sumMatch(subframeShouldBeS,abS) - sumMatch(subframeShouldBeS,abSR)
      falseR <- sumMatch(subframeShouldBeS,abR)
      S <-unlist(lapply(colnames(subframeShouldBeS),function(x){if(grepl(ab,x)==FALSE){nrow(subframeShouldBeS)}else{0}}))
      names(S) <- colnames(subframeShouldBeS)
      ambiguousS <- sumMatch(subframeShouldBeS,abSR)
      notpredictedS <- S - (correctS + falseR + ambiguousS)
      
      # among should be S
      correctR <-sumMatch(subframeShouldBeR,abR)
      falseS <-  sumMatch(subframeShouldBeR,abS) -sumMatch(subframeShouldBeR,abSR)
      R <- unlist(lapply(colnames(subframeShouldBeR),function(x){if(grepl(ab,x)==FALSE){nrow(subframeShouldBeR)}else{0}}))
      names(R) <- colnames(subframeShouldBeR) 
      ambiguousR <- sumMatch(subframeShouldBeR,abSR)
      notpredictedR <- R - (correctR + falseS + ambiguousR)
      
      #among all
      correct <- correctS + correctR
      false <- falseR + falseS
      total <- S + R
      ambiguous <-  ambiguousS + ambiguousR
      notpredicted <- notpredictedS + notpredictedR
         
      predicted <- correct + false
      predorambiguous <- predicted + ambiguous

      significansLevel <- NA
      if(name!="crude"){
        significansLevel <- fixSigni(name)
      }
      
      tempFrame <- data.frame(matrix(nrow=0,ncol=length(returnCols)))

      
      #word <- "AMP_AMC"
      for(word in names(correct)){
        tempFrame <- rbind(tempFrame,
              c(index,"correct",correct[word],ab,significansLevel,word),
              c(index,"correctS",correctS[word],ab,significansLevel,word),
              c(index,"correctR",correctR[word],ab,significansLevel,word),
              c(index,"false",false[word],ab,significansLevel,word),
              c(index,"falseS",falseS[word],ab,significansLevel,word),
              c(index,"falseR",falseR[word],ab,significansLevel,word),
              c(index,"total",total[word],ab,significansLevel,word),
              c(index,"S",S[word],ab,significansLevel,word),
              c(index,"R",R[word],ab,significansLevel,word),
              c(index,"ambiguous",ambiguous[word],ab,significansLevel,word),
              c(index,"ambiguousS",ambiguousS[word],ab,significansLevel,word),
              c(index,"ambiguousR",ambiguousR[word],ab,significansLevel,word),
              c(index,"notpredicted",notpredicted[word],ab,significansLevel,word),
              c(index,"notpredictedS",notpredictedS[word],ab,significansLevel,word),
              c(index,"notpredictedR",notpredictedR[word],ab,significansLevel,word),
              c(index,"predicted",predicted[word],ab,significansLevel,word),
              c(index,"predorambiguous",predorambiguous[word],ab,significansLevel,word)
              )
      }
      colnames(tempFrame) <- c("noinputab","metric","count", "antibiotic","significanceLevel","word")
      colnames(returnFrame) <- returnCols

      returnFrame <- rbind(returnFrame, tempFrame)
    }
  }
  colnames(returnFrame) <- returnCols
  returnFrame$count <- as.integer(returnFrame$count)
  returnFrame  
}

generateAllStatsAlternativeWidest <- function(index,isolateSubset=NA,columns = NA)
{
  #isolateSubset <- c("2","125")
  #index <- 1
  #inputAntibiotics <- c("AMP","AMC","CRO","CIP","OFX","TOB")
  #
  #isolateSubset=NA
  #index <- 4
  #columns <- c("AMP_CAZ_CIP_GEN" ,"AMP_CAZ_CIP_TOB" ,"AMP_CAZ_OFX_GEN")
  
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
    if(length(columns)>1){
      aFrame <- aFrame[,columns]
    }
    aFrame
  })

  
  
  # sumMatch <-function(subframe,pattern){
  #   enLista <- apply(
  #     apply(subframe,c(1,2),
  #           function(x) {
  #             grepl(pattern,x)})
  #     ,2,sum)
  #   
  # }
  
  itemMatch <-function(sample,frame,pattern){
    row <- frame[sample,] 
    match <- as.integer(grepl(pattern,row))
    names(match) <- colnames(frame)
    match
  }
  returnCols <- c("noinputab","metric","count", "antibiotic","significanceLevel","word","sample")
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
      
    
      for(sample in rownames(frame))
      {
        # among should be S
        correctS <- itemMatch(sample,subframeShouldBeS,abS) - itemMatch(sample,subframeShouldBeS,abSR)
        falseR <- itemMatch(sample,subframeShouldBeS,abR)
        S <- itemMatch(sample,subframeShouldBeS,ab)
        ambiguousS <- itemMatch(sample,subframeShouldBeS,abSR)
        notpredictedS <- S - (correctS + falseR + ambiguousS)

        # among should be S
        correctR <-itemMatch(sample,subframeShouldBeR,abR)
        falseS <-  itemMatch(sample,subframeShouldBeR,abS) - itemMatch(sample,subframeShouldBeR,abSR)
        R <- itemMatch(sample,subframeShouldBeR,ab)
        
        ambiguousR <- itemMatch(sample,subframeShouldBeR,abSR)
        notpredictedR <- R - (correctR + falseS + ambiguousR)

        correct <- correctS + correctR
        false <- falseR + falseS
        total <- S + R
        ambiguous <-  ambiguousS + ambiguousR
        notpredicted <- notpredictedS + notpredictedR
        
        predicted <- correct + false
        predorambiguous <- predicted + ambiguous
        
        significansLevel <- NA
        if(name!="crude"){
          significansLevel <- fixSigni(name)
        }
        tempFrame <- data.frame(matrix(nrow=0,ncol=length(returnCols)))

        
        #word <- "AMP_AMC"
        for(word in names(frame)){
          tempFrame <- rbind(tempFrame,
                             c(index,"correct",correct[word],ab,significansLevel,word,sample),
                             c(index,"correctS",correctS[word],ab,significansLevel,word,sample),
                             c(index,"correctR",correctR[word],ab,significansLevel,word,sample),
                             c(index,"false",false[word],ab,significansLevel,word,sample),
                             c(index,"falseS",falseS[word],ab,significansLevel,word,sample),
                             c(index,"falseR",falseR[word],ab,significansLevel,word,sample),
                             c(index,"total",total[word],ab,significansLevel,word,sample),
                             c(index,"S",S[word],ab,significansLevel,word,sample),
                             c(index,"R",R[word],ab,significansLevel,word,sample),
                             c(index,"ambiguous",ambiguous[word],ab,significansLevel,word,sample),
                             c(index,"ambiguousS",ambiguousS[word],ab,significansLevel,word,sample),
                             c(index,"ambiguousR",ambiguousR[word],ab,significansLevel,word,sample),
                             c(index,"notpredicted",notpredicted[word],ab,significansLevel,word,sample),
                             c(index,"notpredictedS",notpredictedS[word],ab,significansLevel,word,sample),
                             c(index,"notpredictedR",notpredictedR[word],ab,significansLevel,word,sample),
                             c(index,"predicted",predicted[word],ab,significansLevel,word,sample),
                             c(index,"predorambiguous",predorambiguous[word],ab,significansLevel,word,sample)
          )
          
        }
        colnames(tempFrame) <- c("noinputab","metric","count", "antibiotic","significanceLevel","word","sample")
        colnames(returnFrame) <- returnCols
        
        returnFrame <- rbind(returnFrame, tempFrame)
      }
    }
  }
  colnames(returnFrame) <- returnCols
  returnFrame$count <- as.integer(returnFrame$count)
  returnFrame  
}



generateAllStatsAlternative <- function(index,isolateSubset=NA)
{
  #isolateSubset <- c("2","125")
  #index <- 1
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
  
  returnCols <- c("noinputab","metric","count", "antibiotic","significanceLevel")
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
          sum(apply(subframe,c(1,2),
                function(x) {
                  grepl(pattern,x)}))
      }
      # sumMatchList <-function(subframe,pattern){
      #   enLista <- apply(
      #     apply(subframe,c(1,2),
      #           function(x) {
      #             grepl(pattern,x)})
      #     ,2,sum)
      #   
      # }
      sum(grepl(ab,colnames(subframeShouldBeR))==FALSE)*nrow(subframeShouldBeR)
      
      # among should be S
      correctS <- sumMatch(subframeShouldBeS,abS) - sumMatch(subframeShouldBeS,abSR)
      falseR <- sumMatch(subframeShouldBeS,abR)
      S <-sum(grepl(ab,colnames(subframeShouldBeS))==FALSE)*nrow(subframeShouldBeS)
      
      ambiguousS <- sumMatch(subframeShouldBeS,abSR)
      notpredictedS <- S - (correctS + falseR + ambiguousS)
      
      # among should be R
      correctR <-sumMatch(subframeShouldBeR,abR)
      falseS <-  sumMatch(subframeShouldBeR,abS) -sumMatch(subframeShouldBeR,abSR)
      R <- sum(grepl(ab,colnames(subframeShouldBeR))==FALSE)*nrow(subframeShouldBeR)
      
      ambiguousR <- sumMatch(subframeShouldBeR,abSR)
      notpredictedR <- R - (correctR + falseS + ambiguousR)
      
      #among all
      correct <- correctS + correctR
      false <- falseR + falseS
      total <- S + R
      ambiguous <-  ambiguousS + ambiguousR
      notpredicted <- notpredictedS + notpredictedR
      
      predicted <- correct + false
      predorambiguous <- predicted + ambiguous
      
      significansLevel <- NA
      if(name!="crude"){
        significansLevel <- fixSigni(name)
      }

      returnFrame <- rbind(returnFrame,
                           c(index,"correct",correct,ab,significansLevel),
                           c(index,"correctS",correctS ,ab,significansLevel),
                           c(index,"correctR",correctR ,ab,significansLevel),
                           c(index,"false",false ,ab,significansLevel),
                           c(index,"falseS",falseS ,ab,significansLevel),
                           c(index,"falseR",falseR ,ab,significansLevel),
                           c(index,"total",total ,ab,significansLevel),
                           c(index,"S",S  ,ab,significansLevel),
                           c(index,"R",R,ab,significansLevel),
                           c(index,"ambiguous",ambiguous,ab,significansLevel),
                           c(index,"ambiguousS",ambiguousS,ab,significansLevel),
                           c(index,"ambiguousR",ambiguousR,ab,significansLevel),
                           c(index,"notpredicted",notpredicted,ab,significansLevel),
                           c(index,"notpredictedS",notpredictedS,ab,significansLevel),
                           c(index,"notpredictedR",notpredictedR,ab,significansLevel),
                           c(index,"predicted",predicted,ab,significansLevel),
                           c(index,"predorambiguous",predorambiguous,ab,significansLevel)      )
    }
  }
  colnames(returnFrame) <- returnCols
  returnFrame$count <- as.integer(returnFrame$count)
  returnFrame  
}


compareAllstats <- function()
{
  # allStatsOrg <-readStatisticsFrame("allstats",getStatisticsFolder(MODE))
  # allStatsNew <-readStatisticsFrame("allStats_sum",getStatisticsFolder(MODE))
  # 
  # allStatsOrgArranged <- allStatsOrg %>% dplyr::arrange_all()
  # allStatsNewArranged <- allStatsNew %>% dplyr::arrange_all()
  # 
  # diff <- allStatsOrgArranged!=allStatsNewArranged
  # same <- allStatsOrgArranged==allStatsNewArranged
  # any(diff,na.rm = TRUE)
  # all(same,na.rm=TRUE)
  

  
  #  allStatsAltWider <-readStatisticsFrame("allStats_sum",getStatisticsFolder(MODE))
  # sum_wider <-  allStatsAltWider %>% 
  #   dplyr::group_by(noinputab,metric,antibiotic,significanceLevel) %>%
  #   dplyr::summarize(count=sum(count),.groups = "drop")
  # 
  # sum_wider <- sum_wider %>% select(noinputab, metric, count, antibiotic, significanceLevel)
  # sum_wider <- sum_wider  %>% dplyr::arrange_all()
  
    
  allStatsAlt <-readStatisticsFrame("allStats_sum",getStatisticsFolder(MODE))
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


errorStatisticsPerAntibiotic <- function(name="allStats_sum",folder = getStatisticsFolder(MODE))
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

    allMetricsCount <- abStats %>% dplyr::group_by(noinputab,significanceLevel) %>% tidyr::pivot_wider(names_from=metric,values_from = count)
    #allMetricsCount %>% dplyr::filter(noinputab==4)
    summary_mutated <- allMetricsCount %>% 
      dplyr::mutate(
        correct_frac=correct/total,
        VME_frac=(falseS+notpredictedR)/total,
        ME_frac=(falseR+notpredictedS)/total
      )
    #summary_mutated %>% dplyr::filter(noinputab==4) %>% dplyr::select(VME_frac,VME_frac_norm,VME_frac_norm_alt,ME_frac,ME_frac_norm,ME_frac_norm_alt)
    if(is.null(rv)){
      rv <- summary_mutated
    } else {
      rv <- rbind(rv,summary_mutated)
    }
  }
  rv
}


errorStatistics <- function(name="allStats_sum",folder = getStatisticsFolder(MODE))
{
  selectedMetrics <- c("correct","ME","VME")
  allMetrics <- c("correct","ME","VME","ambiguous","notpredicted")
  
  allStats <- readStatisticsFrame(name,folder)

  collapse_antibiotics <- function(df) {
    df %>%
      dplyr::group_by(noinputab, metric, significanceLevel) %>%
      dplyr::summarise(count = sum(count, na.rm = TRUE), .groups = "drop")
  }
  
  collapsed <- collapse_antibiotics(allStats)
  
  allMetricsCount <- allStats %>% dplyr::group_by(noinputab,significanceLevel,metric) %>% dplyr::summarize(count=sum(count), .groups = "drop") %>% tidyr::pivot_wider(names_from=metric,values_from = count)
  
#  sanitycheck <- allMetricsCount %>% dplyr::summarize(noinputab, significanceLevel,total, ME_VME_NOTPREDICTED = falseS+falseR+notpredicted, correctAmbigous = correct + ambiguous)
#  sanitycheck2 <- sanitycheck  %>% dplyr::summarize(noinputab, significanceLevel,total,totalAgain = ME_VME_NOTPREDICTED + correctAmbigous)
  summary_mutated <- allMetricsCount %>% 
      dplyr::mutate(
        correct_frac=correct/total,
        VME_frac=(falseS+notpredictedR)/total,
        ME_frac=(falseR+notpredictedS)/total
      )
  #summary_mutated %>%  dplyr::select(noinputab,significanceLevel,correct,VME,ME,VME_frac_norm_alt,ME_frac_norm_alt)
  summary_mutated %>%
    dplyr::arrange(noinputab, significanceLevel)
  
  # F1_metrics <- calcPRMicro(collapsed)
  # rvWithF1 <- cbind(summary_mutated %>%
  #                     arrange(noinputab, significanceLevel),F1_metrics %>% select(precision, recall,f1=f1_micro))
  # #rv %>% dplyr::filter(noinputab==4,is.na(significanceLevel))
  # rvWithF1
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




ALL <-  function()
{
  checkDirs()
  
  GENERATE_METRICS_ALT()
  GENERATE_RISK_STRATIFIED_STATISTICS()
  GENERATE_SAMPLE_METRICS()
  SStatisticsAbSample()
  RStatisticsAbSample()
  
  #subgroup stats
  for(aPrefix in SUBGROUPS){
    GENERATE_METRICS_ALT(paste(aPrefix,"allStats_sum",sep="_"),prefix=aPrefix)
    GENERATE_RISK_STRATIFIED_STATISTICS(paste(aPrefix,"allStats_sum",sep="_"),prefix=aPrefix)
    GENERATE_SAMPLE_METRICS(paste(aPrefix,"allStats_sample",sep="_"),prefix=aPrefix)
    SStatisticsAbSample(paste(aPrefix,"allStats_sample",sep="_"),prefix=aPrefix)
    RStatisticsAbSample(paste(aPrefix,"allStats_sample",sep="_"),prefix=aPrefix)
  }

}


EXTENDED <- function()
{
  checkDirs()
  print("GENERATE_STATISTICS_ALT")
  GENERATE_STATISTICS_ALT()
  print("GENERATE_STATISTICS_ALT_WIDER")
  GENERATE_STATISTICS_ALT_WIDER()
  print("GENERATE_STATISTICS_NON_ATU")
  GENERATE_STATISTICS_NON_ATU()
  print("GENERATE_STATISTICS_ALT_WIDER_ATLEAST_ONE_PER_AB_GROUP")
  GENERATE_STATISTICS_ALT_FROM_WIDER()
  compareAllstats()
  GENERATE_BEST_AND_ALT_METRICS()
  print("nonATUstats")
  GENERATE_METRICS_ALT(name="nonATUstats",prefix="nonATUstats",getStatisticsSubsetFolder(MODE))
  ALL()
}


GENERATE_RISK_STRATIFIED_STATISTICS<- function(name = "allStats_sum",prefix=NA,folder = getStatisticsFolder(MODE))
{
  SStatistics(name,prefix,folder)
  RStatistics(name,prefix,folder)
  SStatisticsAb(name,prefix,folder)
  RStatisticsAb(name,prefix,folder)
}

SStatistics <- function(name = "allStats_sum",prefix=NA,folder = getStatisticsFolder(MODE))
{
  allStats <-readStatisticsFrame(name,folder)  
  subset <- allStats %>% 
    dplyr::group_by(noinputab,significanceLevel,metric) %>% dplyr::summarize(count=sum(count), .groups = "drop")  %>%
    dplyr::filter(metric %in% c("S","correctS","falseR","notpredictedS")) %>%
    tidyr::pivot_wider(names_from=metric,values_from = count)
  fractions <- subset %>% dplyr::mutate(MEfrac=(falseR+notpredictedS)/S,correctfrac=correctS/S)  %>% dplyr::select (noinputab,significanceLevel,correct=correctfrac,ME=MEfrac)
  namebase <- "SStatistics"
  if(!is.na(prefix)){
    namebase <- paste (prefix,namebase,sep="-")
  }
  writeStatisticsFrame(fractions,namebase,folder)  
}

RStatistics <- function(name = "allStats_sum",prefix=NA,folder = getStatisticsFolder(MODE))
{
  allStats <-readStatisticsFrame(name,folder)  
  subset <- allStats %>% 
    dplyr::group_by(noinputab,significanceLevel,metric) %>% dplyr::summarize(count=sum(count), .groups = "drop")  %>%
    dplyr::filter(metric %in% c("R","correctR","falseS","notpredictedR")) %>%
    tidyr::pivot_wider(names_from=metric,values_from = count)
  fractions <- subset %>% dplyr::mutate(VMEfrac=(falseS+notpredictedR)/R,correctfrac=correctR/R)  %>% dplyr::select (noinputab,significanceLevel,correct=correctfrac,VME=VMEfrac)
  namebase <- "RStatistics"
  if(!is.na(prefix)){
    namebase <- paste (prefix,namebase,sep="-")
  }
  writeStatisticsFrame(fractions,namebase,folder)  
}


SStatisticsAb <- function(name = "allStats_sum",prefix=NA,folder = getStatisticsFolder(MODE))
{
  allStats <-readStatisticsFrame(name,folder)  
  subset <- allStats %>% 
    dplyr::filter(metric %in% c("S","correctS","falseR","notpredictedS")) %>%
        tidyr::pivot_wider(names_from=metric,values_from = count) 
  fractions <- subset %>% 
    dplyr::mutate(MEfrac=(falseR+notpredictedS)/S,correctfrac=correctS/S)  %>% 
    dplyr::select (noinputab,significanceLevel,antibiotic,correct=correctfrac,ME=MEfrac) %>% 
    reorderAntibioticsRows()
  namebase <- "SStatisticsAb"
  if(!is.na(prefix)){
    namebase <- paste (prefix,namebase,sep="-")
  }
  writeStatisticsFrame(fractions,namebase,folder)  
}



RStatisticsAb <- function(name = "allStats_sum",prefix=NA,folder = getStatisticsFolder(MODE))
{
  allStats <-readStatisticsFrame(name,folder)  
  subset <- allStats %>% 
    dplyr::filter(metric %in% c("R","correctR","falseS","notpredictedR")) %>%
    tidyr::pivot_wider(names_from=metric,values_from = count)
  fractions <- subset %>% dplyr::mutate(VMEfrac=(falseS+notpredictedR)/R,correctfrac=correctR/R)  %>% dplyr::select (noinputab,significanceLevel,antibiotic,correct=correctfrac,VME=VMEfrac) %>% 
    reorderAntibioticsRows()
  namebase <- "RStatisticsAb"
  if(!is.na(prefix)){
    namebase <- paste (prefix,namebase,sep="-")
  }
  writeStatisticsFrame(fractions,namebase,folder)  
}


RStatisticsAbSample <- function(name = "allStats_sample",prefix=NA,folder = getStatisticsFolder(MODE))
{
  allStats <-readStatisticsFrame(name,folder)  
  subset <- allStats %>% 
    dplyr::filter(metric %in% c("R","correctR","falseS","notpredictedR")) %>%
    tidyr::pivot_wider(names_from=metric,values_from = count)
  fractions <- subset %>% 
    dplyr::mutate(VMEfrac=ifelse(test=R>0,yes=(falseS+notpredictedR)/R,no=0),correctfrac=ifelse(test=(R>0),yes=correctR/R,no=0))  %>% 
    dplyr::select (noinputab,sample,significanceLevel,antibiotic,correct=correctfrac,VME=VMEfrac) %>%
    reorderAntibioticsRows()
  namebase <- "RStatisticsAbSample"
  if(!is.na(prefix)){
    namebase <- paste (prefix,namebase,sep="-")
  }
  writeStatisticsFrame(fractions,namebase,folder)  
}

SStatisticsAbSample <- function(name = "allStats_sample",prefix=NA,folder = getStatisticsFolder(MODE))
{
  allStats <-readStatisticsFrame(name,folder)  
  subset <- allStats %>% 
    dplyr::filter(metric %in% c("S","correctS","falseR","notpredictedS")) %>%
    tidyr::pivot_wider(names_from=metric,values_from = count)
  fractions <- subset %>% 
    dplyr::mutate(MEfrac=ifelse(test=S>0,yes=(falseR+notpredictedS)/S,no=0),correctfrac=ifelse(test=(S>0),yes=correctS/S,no=0))  %>% 
    dplyr::select (noinputab,sample,significanceLevel,antibiotic,correct=correctfrac,ME=MEfrac) %>%
    reorderAntibioticsRows()
  namebase <- "SStatisticsAbSample"
  if(!is.na(prefix)){
    namebase <- paste (prefix,namebase,sep="-")
  }
  writeStatisticsFrame(fractions,namebase,folder)  
}

# Moved metric calc to manuscript.
#
# calcPRMicro <- function(
#     df = readStatisticsFrame("allStats_sum", getStatisticsFolder(MODE))
# ) {
#   library(dplyr)
#   library(tidyr)
#   
#   safe_div <- function(num, den) ifelse(den == 0, NA_real_, num / den)
#   
#   # Detect count column
#   cnt_col <- intersect(names(df), c("value", "count", "n"))[1]
#   if (is.na(cnt_col)) stop("No count column found (expected one of: value, count, n).")
#   
#   has_abx <- "antibiotic" %in% names(df)
#   
#   # Build grouping variables dynamically
#   group_vars <- c(if (has_abx) "antibiotic", "significanceLevel", "noinputab", "metric")
#   
#   # Keep only relevant metrics and aggregate
#   core <- df %>%
#     dplyr::filter(metric %in% c("correctR", "correctS", "falseR", "falseS")) %>%
#     dplyr::group_by(dplyr::across(all_of(group_vars))) %>%
#     dplyr::summarise(n = sum(.data[[cnt_col]], na.rm = TRUE), .groups = "drop")
#   
#   # Pivot to wide
#   wide <- tidyr::pivot_wider(core,
#                              names_from = metric,
#                              values_from = n,
#                              values_fill = 0)
#   
#   # Ensure columns exist
#   for (nm in c("correctR","correctS","falseR","falseS")) {
#     if (!nm %in% names(wide)) wide[[nm]] <- 0
#   }
#   
#   # Compute precision, recall, micro-F1
#   if (has_abx) {
#     out <- wide %>%
#       dplyr::transmute(
#         antibiotic,
#         significanceLevel,
#         noinputab,
#         TP = correctR,
#         TN = correctS,
#         FP = falseR,
#         FN = falseS,
#         precision = safe_div(TP, TP + FP),
#         recall    = safe_div(TP, TP + FN),
#         f1_micro  = safe_div(2 * precision * recall, precision + recall)
#       ) %>%
#       dplyr::arrange(noinputab, antibiotic, significanceLevel)
#   } else {
#     out <- wide %>%
#       dplyr::transmute(
#         significanceLevel,
#         noinputab,
#         TP = correctR,
#         TN = correctS,
#         FP = falseR,
#         FN = falseS,
#         precision = safe_div(TP, TP + FP),
#         recall    = safe_div(TP, TP + FN),
#         f1_micro  = safe_div(2 * precision * recall, precision + recall)
#       ) %>%
#       dplyr::arrange(noinputab, significanceLevel)
#   }
#   
#   out
# }




# F1_Score <- function()
# {
#   metrics <- calcF1Metrics()
#   NA_n6 <- metrics %>% filter(is.na(significanceLevel))  %>% filter(noinputab == 6) %>% select(antibiotic,f1=f1_micro)
#   NA_inputscore <- metrics %>% filter(is.na(significanceLevel)) %>% select(antibiotic,noinputab,f1=f1_micro)
#   
#   anothermetrics <- calcPRMicro()
#   writeStatisticsFrame(frame=NA_inputscore,name = "f1_noinputab",getStatisticsFolder(MODE))
# }

