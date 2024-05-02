source(file = 'model/modelcommon.R')
library(stringr)

ROWS <- 99

PATTERN_S <- "^([[:alpha:]]+)_S$"
PATTERN_R <- "^([[:alpha:]]+)_R$"
PATTERN_SR <- "^([[:alpha:]]+)_SR$"

unit_test <- function()
{
  a <- "CRO_S LVX_R TZP_R AMC_R MFX_S GEN_S CTX_S PIP_R FEP_R CAZ_S CIP_S OFX_S TOB_R"
  p <- "CRO_SR LVX_SR TZP_R AMC_S MFX_R GEN_S CTX_SR PIP_SR FEP_R CAZ_R CIP_S"

  ok = TRUE
  ok = ok & getAmbiguous(a,p) == "CRO LVX CTX PIP"
  ok = ok & getCorrect(a,p) == "TZP FEP GEN CIP"
  ok = ok & getME(a,p) == "MFX CAZ"
  ok = ok & getVME(a,p) == "AMC"
  ok = ok & getPredictedOrAmbiguous(a,p) == "AMC GEN CIP TZP MFX FEP CAZ CRO LVX CTX PIP"
  ok = ok & getPredicted(a,p) == "AMC GEN CIP TZP MFX FEP CAZ"
  ok = ok & getNotPredicted(a,p) == "OFX TOB"
  ok
}

getNotPredicted <- function(a,p){
  p <- unlist(strsplit(p,split=" "))
  a <- unlist(strsplit(a,split=" "))
  allP <- substr(p,1,3)
  allA <- substr(a,1,3)
  paste(allA[!(allA %in% allP)],collapse=" ")
}


getVME <- function(a,p){
  p <- unlist(strsplit(p,split=" "))
  a <- unlist(strsplit(a,split=" "))
  # for vme a should be R and p should be S
  sensP <- substr(unlist(regmatches(p, gregexpr(PATTERN_S, p))),1,3)
  resA <- substr(unlist(regmatches(a, gregexpr(PATTERN_R, a))),1,3)
  paste(intersect(sensP,resA),collapse=" ")
}

getME <- function(a,p){
  p <- unlist(strsplit(p,split=" "))
  a <- unlist(strsplit(a,split=" "))
  # for vme a should be S and p should be R
  resP <- substr(unlist(regmatches(p, gregexpr(PATTERN_R, p))),1,3)
  sensA <- substr(unlist(regmatches(a, gregexpr(PATTERN_S, a))),1,3)
  paste(intersect(resP,sensA),collapse=" ")
}

getAmbiguous <- function(a,p){
  p <- unlist(strsplit(p,split=" "))
  # all that are _SR
  ambP <- substr(unlist(regmatches(p, gregexpr(PATTERN_SR, p))),1,3)
  paste(ambP,collapse=" ")
}

getCorrect <- function(a,p){
  p <- unlist(strsplit(p,split=" "))
  a <- unlist(strsplit(a,split=" "))
  # all that are not _SR and are the same
  sensP <- substr(unlist(regmatches(p, gregexpr(PATTERN_S, p))),1,3)
  resA <- substr(unlist(regmatches(a, gregexpr(PATTERN_R, a))),1,3)
  resP <- substr(unlist(regmatches(p, gregexpr(PATTERN_R, p))),1,3)
  sensA <- substr(unlist(regmatches(a, gregexpr(PATTERN_S, a))),1,3)

  same <- union(intersect(resA,resP),intersect(sensA,sensP))
  paste(same,collapse=" ")
}

getPredicted <- function(a,p){
  p <- unlist(strsplit(p,split=" "))
  # all that are not _SR and are the same
  sensP <- substr(unlist(regmatches(p, gregexpr(PATTERN_S, p))),1,3)
  resP <- substr(unlist(regmatches(p, gregexpr(PATTERN_R, p))),1,3)
  paste(union(sensP,resP),collapse=" ")
}

getPredictedOrAmbiguous <- function(a,p){
  p <- unlist(strsplit(p,split=" "))
  # all that are not _SR and are the same
  sensP <- substr(unlist(regmatches(p, gregexpr(PATTERN_S, p))),1,3)
  resP <- substr(unlist(regmatches(p, gregexpr(PATTERN_R, p))),1,3)
  sensresP <- substr(unlist(regmatches(p, gregexpr(PATTERN_SR, p))),1,3)
  paste(union(sensP,union(resP,sensresP)),collapse=" ")
}



COMPARE_CORRECT <- "correct"
COMPARE_VME <- "VME"
COMPARE_ME <- "ME"
COMPARE_PREDICTED <- "predicted"
COMPARE_AMBIGUOUS <- "ambiguous"
COMPARE_PREDICTED_OR_AMBIGUOUS <- "predorambiguous"
COMPARE_NOT_PREDICTED <- "notpredicted"

COMPARE_METRICS_VECTOR <-c(
  COMPARE_CORRECT,
  COMPARE_VME,
  COMPARE_ME,
  COMPARE_PREDICTED,
  COMPARE_AMBIGUOUS,
  COMPARE_PREDICTED_OR_AMBIGUOUS,
  COMPARE_NOT_PREDICTED
)

#COMPARE_FUNCTIONS <- c(COMPARE_CORRECT=getCorrect, COMPARE_VME = getVME,COMPARE_ME=getME,COMPARE_PREDICTED=getPredicted,COMPARE_AMBIGUOUS=getAmbiguous,COMPARE_PREDICTED_OR_AMBIGUOUS=getPredictedOrAmbiguous)
COMPARE_FUNCTIONS <- c(getCorrect, getVME,getME,getPredicted,getAmbiguous,getPredictedOrAmbiguous,getNotPredicted)
names(COMPARE_FUNCTIONS)<-COMPARE_METRICS_VECTOR



GENERATE_COMPARES <- function()
{
  for(index in range){
    predsFiles <- makeShortFileName("preds",index)
    for(signi in signis){
      predsFiles <- c(predsFiles,list.files(PREDICTION_FOLDER,pattern = makeShortFileName("confpreds",index,signi)))
    }
    answerFile <-  makeShortFileName("answer",index)
    for(predsFile in predsFiles){
      predictionDataFrame <- readOutputFrameShort(short=predsFile,folder = PREDICTION_FOLDER)
      answerDataFrame <- readOutputFrameShort(short=answerFile,folder = PREDICTION_FOLDER)
      for ( name in names(COMPARE_FUNCTIONS)){
        f=COMPARE_FUNCTIONS[[name]]
        frame <- applyElementWise(x=answerDataFrame,y=predictionDataFrame,aFunction = f)
        newShortName <- paste(substring(predsFile,0,nchar(predsFile)-4),"_",name,".csv",sep="")
        writeOutputFrameShort(frame,shortfile=newShortName,folder=COMPARE_FOLDER)
      }
    }
  }
}

GENERATE_STATISTICS <- function()
{
  allStats <- data.frame(
    noinputab = integer(),
    metric = character(),
    count = integer(),
    antibiotic = character(),
    significanslevel = numeric()
  )

  for(index in range){
    predsFiles <- makeShortFileName("preds",index)
    for(signi in signis){
      predsFiles <- c(predsFiles,makeShortFileName("confpreds",index,signi))
    }
    names(predsFiles) <- c("crude",signis)
    for(predsName in names(predsFiles)){
      predsFile <- predsFiles[[predsName]]
      for ( name in names(COMPARE_FUNCTIONS)){
        newShortName <- paste(substring(predsFile,0,nchar(predsFile)-4),"_",name,".csv",sep="")
        frame <- readOutputFrameShort(shortfile=newShortName,folder=COMPARE_FOLDER)
        frame[is.na(frame)] <- ""
        for(ab in ANTIBIOTICS){
          count = countWords(frame,ab)
          if(predsName=="crude"){
            signi = NA
          }
          else {
            signi <- fixSigni(predsName)
          }
          row <- c (
            as.integer(index),
            name,
            count,
            ab,
            signi)
          allStats[nrow(allStats) + 1, ] <- row
        }
      }
    }
    writeStatisticsFrame(frame = allStats, name = paste("stats_temp",as.integer(index),sep="_"),folder = STATISTICS_FOLDER)
  }
  writeStatisticsFrame(allStats, name = "allstats",STATISTICS_FOLDER)
}


generate_signi_statistics <- function(technical=FALSE)
{
  allStats <- readStatisticsFrame("allstats",STATISTICS_FOLDER)
  range <- unique(allStats$noinputab)

  rows <- unique(allStats$significanslevel)
  rows[is.na(rows)]<-"crude"
  
  percentageFrame <- data.frame(row.names = rows)

  # Conformal predictions clinical & technical
  for(index in range) {
    percentageRow <- c()
    levels <- unique(allStats$significanslevel)
    subsetIndex <- allStats[allStats$noinputab==index,]
    for (level in levels){
      if(is.na(level)) {
        subset <- subsetIndex[is.na(subsetIndex$significanslevel),]
      }
      else {
        subset <- subsetIndex[!is.na(subsetIndex$significanslevel) & subsetIndex$significanslevel==level,]
      }
      #      print(sum(subset$count))
      if(technical){
        correctTechnical <- subset[subset$metric==COMPARE_CORRECT | subset$metric==COMPARE_AMBIGUOUS,] #   <- FIXME spelling 
        predictedTechnical <- subset[subset$metric==COMPARE_PREDICTED_OR_AMBIGUOUS,]
        percentageTechnical <- sum(correctTechnical$count)/sum(predictedTechnical$count)
        percentage <- percentageTechnical
      }
      else {
        correctClinical <- subset[subset$metric==COMPARE_CORRECT,]
        predictedClinical <- subset[subset$metric==COMPARE_PREDICTED,]
        percentageClinical <- sum(correctClinical$count)/sum(predictedClinical$count)
        percentage <- percentageClinical
      }
      percentageRow <- c(percentageRow,percentage)
#      print(paste(index,level,percentage,sep=" "))
    }
    percentageFrame <- cbind(percentageFrame,percentageRow)
  }
  colnames(percentageFrame) <- range
  significanceLevel <- rownames(percentageFrame)
  percentageFrame <- cbind(significanceLevel,percentageFrame)
  percentageFrame
}


generate_signi_statistics_details <- function(level = "0.025")
{
  allStats <- readStatisticsFrame("allstats",STATISTICS_FOLDER)
  range <- unique(allStats$noinputab)
  
  percentageFrame <- data.frame(row.names = c("correct","amb","vme","me","nopred"))

  for(index in range) {
    percentageRow <- c()
    subsetIndex <- allStats[allStats$noinputab==index,]
    if(is.na(level)) {
      subset <- subsetIndex[is.na(subsetIndex$significanslevel),]
    }
    else {
      subset <- subsetIndex[!is.na(subsetIndex$significanslevel) & subsetIndex$significanslevel==level,]
    }
    correct <- sum(subset[subset$metric==COMPARE_CORRECT,]$count)
    amb <- sum(subset[subset$metric==COMPARE_AMBIGUOUS,]$count)
    vme <- sum(subset[subset$metric==COMPARE_VME,]$count)
    me <- sum(subset[subset$metric==COMPARE_ME,]$count)
    nopred <- sum(subset[subset$metric==COMPARE_NOT_PREDICTED,]$count)
      
    percentageRow <- c(correct,amb,vme,me,nopred)
    
    ### columns X rows X number_of_possible_predictions
    maxPredictions = PASCAL[[index]]*ROWS*(NUMBER_OF_ANTIBIOTICS-index)
    sum(percentageRow)
    percentageRow <- percentageRow/maxPredictions
      #      print(paste(index,level,percentage,sep=" "))
    percentageFrame <- cbind(percentageFrame,percentageRow)
  }
  colnames(percentageFrame) <- range
  significanceLevel <- rownames(percentageFrame)
  percentageFrame <- cbind(significanceLevel,percentageFrame)
  percentageFrame
}



generate_ab_statistics <- function(metric = COMPARE_CORRECT)
{
  allStats <- readStatisticsFrame("allstats",STATISTICS_FOLDER)
  range <- unique(allStats$noinputab)
  abs = unique(allStats$antibiotic)

  rows <- unique(allStats$antibiotic)

  abFrame <- data.frame(row.names = rows)
  
  for(index in range) {
    abRow <- c()
    for(ab in abs) {
      correct <- allStats[allStats$noinputab==index & allStats$metric==metric & is.na(allStats$significanslevel) & allStats$antibiotic==ab,]
      predicted <- allStats[allStats$noinputab==index & allStats$metric==COMPARE_PREDICTED & is.na(allStats$significanslevel) & allStats$antibiotic==ab ,]
      percentage <- sum(correct$count)/sum(predicted$count)
      abRow <- c(abRow,percentage)
    }
    abFrame <- cbind(abFrame,abRow)
  }
  colnames(abFrame) <- range
  antibiotic <- rownames(abFrame)
  abFrame <- cbind(antibiotic,abFrame)
  abFrame
}
  

generate_ab_statistics_for_index <- function(index)
{
  allStats <- readStatisticsFrame("allstats",STATISTICS_FOLDER)
  range <- unique(allStats$noinputab)
  abs = unique(allStats$antibiotic)
  
  rows <- unique(allStats$antibiotic)
  
  abFrame <- data.frame(row.names = rows)
  
  metrics <- c(COMPARE_CORRECT,COMPARE_ME,COMPARE_VME)
  for(metric in metrics) {
    abRow <- c()
    for(ab in abs) {
      correct <- allStats[allStats$noinputab==index & allStats$metric==metric & is.na(allStats$significanslevel) & allStats$antibiotic==ab,]
      predicted <- allStats[allStats$noinputab==index & allStats$metric==COMPARE_PREDICTED & is.na(allStats$significanslevel) & allStats$antibiotic==ab ,]
      percentage <- sum(correct$count)/sum(predicted$count)
      abRow <- c(abRow,percentage)
    }
    abFrame <- cbind(abFrame,abRow)
  }
  colnames(abFrame) <- metrics
  antibiotic <- rownames(abFrame)
  abFrame <- cbind(antibiotic,abFrame)
  abFrame
}




countWordsScalar <-function(scalar,split=" ")
{
  length(unlist(strsplit(scalar,split=split)))
}

countWordsPerColumnFrame <-function(aFrame,split=" ")
{
  lapply(aFrame,function(x) { sum(countWordsScalar(x,split))})
}

countWordsPerRowFrame <-function(aFrame,split=" ")
{
  apply(aFrame,1,function(x) { sum(countWordsScalar(x,split))})
}





countWordsFrame <-function(aFrame,split=" ")
{
  resultMatrix <- mapply(function(u){mapply(function(x){countWordsScalar(x,split)},u)},aFrame)
  resultFrame <- as.data.frame(resultMatrix)
  rownames(resultFrame) <- rownames(aFrame)
  resultFrame
}

#aFrame <- frames[["VME"]]

sieveWordsFrame <-function(aFrame,ab)
{
  aFunction <- function(x){
    rv <- ""
    if(length(grep(ab,x))>0){
      rv <- ab
    }
    rv}
  
  resultMatrix <- mapply(function(u){mapply(aFunction,u)},aFrame)
  resultFrame <- as.data.frame(resultMatrix)
  rownames(resultFrame) <- rownames(aFrame)
  resultFrame
}



getCompareFrame <-function(comparetype,index,signi)
{
  if(is.na(signi)){
    predsFile <- makeShortFileName("preds",index)
  } else {
    predsFile <- makeShortFileName("confpreds",index,signi)
    
  }
  shortName <- paste(substring(predsFile,0,nchar(predsFile)-4),"_",comparetype,".csv",sep="")
  frame <- readOutputFrameShort(shortfile=shortName,folder=COMPARE_FOLDER)
  frame[is.na(frame)]<-""
  frame  
}

getCounts <-function(comparetype,index,signi)
{
  frame <- getCompareFrame(comparetype,index,signi)
  counts <- countWordsPerColumnFrame(frame)
  counts
}



findBestCombinations<- function(signi=NA,fromBestCrude=FALSE)
{
   
  bestFrame <- data.frame()
  
  noinputab <- c()
  percentage <- c()
  combinations <- c()
  correct <- c()
  predictions <- c()
  maxpredictions <-c()
  fractionCorrect <- c()
  fractionPredictions <-c()
  for(index in range){
    
    correctCounts <- getCounts(COMPARE_CORRECT,index,signi)
    predsCounts <- getCounts(COMPARE_PREDICTED,index,signi)

    correctCountsForBest <- correctCounts
    predsCountsForBest <- predsCounts
    if(!is.null(signi) & fromBestCrude){
      correctCountsForBest <- getCounts(COMPARE_CORRECT,index,NA)
      predsCountsForBest <- getCounts(COMPARE_PREDICTED,index,NA)
    }
    
    # 2:  1188*   99*(14-2)
    # maxPerCol  <- maxPredictions(correctFrame,index)/ncol(correctFrame)

    ## First find best
    percentages <- as.double(correctCountsForBest)/as.double(predsCountsForBest)
    names(percentages) <- names(correctCountsForBest)
    bestNames <- names(percentages[percentages==max(percentages)])

    anotherPercentages <- as.double(correctCounts)/as.double(predsCounts)
    names(anotherPercentages)<-names(correctCounts)
    maxPerColumn <- ROWS * (NUMBER_OF_ANTIBIOTICS-index) 
    for (aBestName in bestNames) {
      aCorrect  <- as.integer(correctCounts[aBestName])
      anotherMaxP <- as.double(anotherPercentages[aBestName])
      aPrediction <- as.integer(predsCounts[aBestName])
      aCorrectFraction = as.double(correctCounts[aBestName])/maxPerColumn
      aPredictionFraction = as.double(predsCounts[aBestName])/maxPerColumn
    
      correct <- c(correct,aCorrect)
      predictions <- c(predictions,aPrediction)
      percentage <- c(percentage,anotherMaxP)
      combinations <- c(combinations,aBestName)
      maxpredictions <- c(maxpredictions,maxPerColumn)
      fractionCorrect <-c(fractionCorrect,aCorrectFraction)
      fractionPredictions <-c(fractionPredictions,aPredictionFraction)
      noinputab <- c(noinputab,index)
    }
  }
  bestFrame <- cbind(noinputab,as.character(percentage),combinations,as.character(correct),as.character(predictions),as.character(maxpredictions),as.character(fractionCorrect),as.character(fractionPredictions))
    
  colnames(bestFrame) <- c("noinputab","percentage","combinations","correct","predictions","maxpredictions","fraction_correct","fraction_predictions")
  
  bestFrame
}
  
  

GENERATE_SUB_STATISTICS <- function()
{
  technicalSigniStatistics <- generate_signi_statistics(TRUE)
  clinicalSigniStatistics <- generate_signi_statistics(FALSE)
  writeStatisticsFrame(technicalSigniStatistics,name = "technicalSigniStatistics",STATISTICS_FOLDER)
  writeStatisticsFrame(clinicalSigniStatistics,name = "clinicalSigniStatistics",STATISTICS_FOLDER)

  for(signi in signis){
    statistics <- generate_signi_statistics_details(fixSigni(signi))
    writeStatisticsFrame(statistics,name = paste("signiStatisticsDetails",signi,sep="-"),STATISTICS_FOLDER)
  }
  
  
  
  abStatistics <- generate_ab_statistics(metric = COMPARE_CORRECT)
  writeStatisticsFrame(abStatistics,name = "abStatistics-correct",STATISTICS_FOLDER)

  abStatistics <- generate_ab_statistics(metric = COMPARE_VME)
  writeStatisticsFrame(abStatistics,name = "abStatistics-vme",STATISTICS_FOLDER)
  
  abStatistics <- generate_ab_statistics(metric = COMPARE_ME)
  writeStatisticsFrame(abStatistics,name = "abStatistics-me",STATISTICS_FOLDER)
  
  for(index in range){
    abStatistics <- generate_ab_statistics_for_index(index = index)
    writeStatisticsFrame(abStatistics,name = paste("abStatistics",index,sep="-"),STATISTICS_FOLDER)
  }
  
}


GENERATE_SPECIFIC_STATISTICS <- function()
{
  generate_best_statistics()
  generate_sample_statistics()
  generate_ab_vs_metric("selected")
  generate_ab_vs_metric("best")

  generate_sample_vs_metric()
}



generate_best_statistics <- function()
{
  bestCombinationStatistics <- findBestCombinations()
  writeStatisticsFrame(bestCombinationStatistics,name = "bestStatistics",STATISTICS_BEST_FOLDER)
  
  bestWords <- bestCombinationStatistics[,"combinations"]   
  dumpBestWordsFrame(bestWords,"bestWordFrame")
  for(signi in signis){
    bestCombinationStatistics <- findBestCombinations(signi = signi,fromBestCrude = FALSE)
    bestWords <- c(bestWords,bestCombinationStatistics[,"combinations"])
    writeStatisticsFrame(bestCombinationStatistics,name = paste("bestStatistics",signi,sep="_"),STATISTICS_BEST_FOLDER)
    
    bestCombinationStatistics <- findBestCombinations(signi = signi,fromBestCrude = TRUE)
    writeStatisticsFrame(bestCombinationStatistics,name = paste("bestStatistics_frombestnonconformal",signi,sep="_"),STATISTICS_BEST_FOLDER)
  }
  dumpBestWordsFrame(bestWords,"bestWordFrame_frombestnonconformal")
  
}

getComplementAntibiotics <-function(word)
{
  ANTIBIOTICS
  remainingAntibiotics <- lapply( ANTIBIOTICS, function(x){
    length(grep(x,word))>0
  }
  )
  ANTIBIOTICS[remainingAntibiotics!=TRUE]
}

generate_sample_vs_metric <-function()
{
  svm_range <- 4:8
  
  for (index in svm_range) {
    for (signi in c(NA,signis)){
      for(ab in c(NA,ANTIBIOTICS)){
        dump_sample_vs_metrics(index,signi,ab)
      }
    }
  }  
  for(index in range){
    dump_sample_vs_metrics(index,NA,"AMC")
  }
}

dump_sample_vs_metrics <- function(index,signi=NA,ab=NA)
{
  if(is.na(signi)){
    metrics <- c(COMPARE_VME,COMPARE_ME,COMPARE_CORRECT)
  } else {
    metrics <- c(COMPARE_VME,COMPARE_ME,COMPARE_CORRECT,COMPARE_AMBIGUOUS,COMPARE_NOT_PREDICTED)
  }
  numberOfRows <- ROWS
    numberOfColumns <- length(metrics)
    resultFrame <- data.frame(matrix(NA, nrow = numberOfRows, ncol = numberOfColumns))
    colnames(resultFrame) <- metrics
    for(metric in metrics) {
      compareFrame <- getCompareFrame(metric,index,signi)
      if(!is.na(ab)){
        compareFrame <- sieveWordsFrame(compareFrame,ab)
      }
      counts <- countWordsPerRowFrame(compareFrame)
      resultFrame[,metric] <- counts
      rownames(resultFrame) <- rownames(compareFrame)
    }
    
    if(is.na(ab)){
      maxPerRow <- PASCAL[[index]] * (NUMBER_OF_ANTIBIOTICS-index)
    } else {
      maxPerRow <- PASCAL[[index]] * (NUMBER_OF_ANTIBIOTICS-index)/NUMBER_OF_ANTIBIOTICS
    }
    apply(resultFrame,1,sum)
    
    resultFrame <- resultFrame/maxPerRow
    resultFrame <- resultFrame[order(resultFrame$VME+resultFrame$ME,decreasing = TRUE),]
    resultFrame <- resultFrame[order(resultFrame$correct,decreasing = FALSE),]
    resultFrame <- resultFrame[order(resultFrame$VME,decreasing = TRUE),]
    
    base <- "samplemetrics"
    if(!is.na(ab)){
      base <- paste(base,ab,sep="-")
    }
    writeStatisticsFrame(resultFrame,name=paste(base,fixSigni(signi),index,sep="-"),folder = STATISTICS_SAMPLE_METRICS_FOLDER,row.names = TRUE)

}

generate_ab_vs_metric <-function(tag="selected")
{

  for(signi in c(NA,signis)) {
    metrics <- c(COMPARE_CORRECT,COMPARE_ME,COMPARE_VME,COMPARE_AMBIGUOUS)
    frames <- list(metrics)
    for(metric in metrics){
      name <- paste("samplewords",metric,tag,fixSigni(signi),"words",sep = "-")  ## FIXME use signi
      frame <- readStatisticsFrame(name,STATISTICS_SAMPLE_WORDS_FOLDER)
      frames[[metric]] <- frame
    }
  
    all <- colnames(frames[[COMPARE_CORRECT]])[2:ncol(frames[[COMPARE_CORRECT]])]

    for(word in all)  {
    #  word <- all[[5]]
      complement <- getComplementAntibiotics(word)
      numberOfRows <- length(complement)
      numberOfColumns <- length(metrics)
      resultFrame <- data.frame(matrix(NA, nrow = numberOfRows, ncol = numberOfColumns))
      row.names(resultFrame) <- complement
      colnames(resultFrame) <- metrics
        
      for (metric in metrics){
        aFrame <- frames[[metric]]
        for(ab in complement){
          sievedFrame <- sieveWordsFrame(aFrame,ab)
          countAbInFrame <- sum(sievedFrame[,word]==ab)
          resultFrame[ab,metric] <- countAbInFrame
        }
      }
      antibiotics <- row.names(resultFrame)
      resultFrame <- cbind(antibiotics,resultFrame)
      writeStatisticsFrame(resultFrame,paste("metrics",fixSigni(signi),word,sep="-"),STATISTICS_ANTIBIOTICS_METRICS_FOLDER,row.names = FALSE)
    }
  }  
}


dumpBestWordsFrame <- function(bestWords,name)
{
  words <- unique(sort(bestWords))
  noinputab <- unlist(lapply(words, function(word){ length(unlist(strsplit(word,split = "_")))}))
  bestWordsFrame <- data.frame(noinputab,words)
  bestWordsFrame <- bestWordsFrame[order(noinputab),]
  writeStatisticsFrame(bestWordsFrame,name,STATISTICS_BEST_FOLDER)  
}


listToWordsFrame <- function(sentences)
{
  aFrame <- data.frame(sentences)
  word <- sentences
  noinputab <- unlist(lapply(sentences,function(x){countWordsScalar(x,split="_")}) )
  data.frame(noinputab,word)
}


generate_sample_statistics <-function()
{
  
  bestWordsFrame <- readStatisticsFrame("bestWordFrame",STATISTICS_BEST_FOLDER)
  for(signi in c(NA,signis)){
    abs_and_counts_best <- inspect_samples(bestWordsFrame,signi)
    dumpSampleWords(abs_and_counts_best,paste("best",fixSigni(signi),sep="-"))
  }
  if(NUMBER_OF_ANTIBIOTICS==14){  
    sentences <- c ("AMP_CRO_CIP_TOB","AMP_AMC_CRO_CIP_TOB", "AMP_AMC_CRO_CIP_OFX_TOB","AMC_PIP_TZP_CRO_CIP_MFX_GEN","AMC_PIP_TZP_CAZ_CRO_CIP_OFX_TOB")
  }
  else {
    sentences <- c ("AMP_CRO_CIP_TOB",
    "AMP_CRO_OFX_LVX_TOB",
    "AMP_CAZ_CRO_OFX_LVX_TOB",
    "AMP_TZP_CAZ_CRO_OFX_LVX_TOB",
    "AMP_TZP_CAZ_CRO_FEP_CIP_OFX_GEN")
  }
  selectedtWordsFrame <- listToWordsFrame(sentences)
  for(signi in c(NA,signis)){
    abs_and_counts_selected <- inspect_samples(selectedtWordsFrame,signi)
    dumpSampleWords(abs_and_counts_selected,paste("selected",fixSigni(signi),sep="-"))
  }
  
}




dumpSampleWords <- function(abs_and_counts,name){
  for(metric in COMPARE_METRICS_VECTOR){
    name_base <- paste("samplewords",metric,sep="-")
    wordsList <- abs_and_counts[[1]]
    metricWordsFrame <- wordsList[[metric]]
    writeStatisticsFrame(metricWordsFrame,paste(name_base,name,"words",sep="-"),STATISTICS_SAMPLE_WORDS_FOLDER,row.names = TRUE)
    
    countsList <- abs_and_counts[[2]]
    metricCountsFrame <- countsList[[metric]]
    writeStatisticsFrame(metricCountsFrame,paste(name_base,name,"counts",sep="-"),STATISTICS_SAMPLE_WORDS_FOLDER,row.names = TRUE)
  }
}

inspect_samples <- function(wordsFrame,signi=NA){
  allnoinputab  <- as.integer(unlist(unique(wordsFrame$noinputab)))
  countList <- list(COMPARE_METRICS_VECTOR)
  wordsList <- list(COMPARE_METRICS_VECTOR)
  noinputab <- allnoinputab[[1]]
  comparemetrics <- COMPARE_CORRECT
  for(comparemetrics in COMPARE_METRICS_VECTOR) {
    outputFrameWords <- NULL
    for(noinputab in allnoinputab){
        words <- wordsFrame[wordsFrame$noinputab==noinputab,]$word
        frame <- getCompareFrame(comparemetrics,noinputab,signi)
        selectedFrame <- frame[words]
        if(is.null(outputFrameWords)){
          outputFrameWords <- data.frame(row.names=rownames(selectedFrame))
        }
        outputFrameWords <- cbind(outputFrameWords,selectedFrame)
    }
    outputFrameCount <- countWordsFrame(outputFrameWords)
    
    item <- list(outputFrameWords)
    names(item)[1]<-comparemetrics
    wordsList[comparemetrics] <- item
    
    item <- list(outputFrameCount)
    names(item)[1]<-comparemetrics
    countList[comparemetrics] <- item
  }
  list(wordsList,countList)
} 


ALL <- function()
{
  stopifnot(dir.exists(PREDICTION_BASE))
  
  
  if(!dir.exists(COMPARE_FOLDER)){
    dir.create(COMPARE_FOLDER)
  }
  if(!dir.exists(STATISTICS_FOLDER)){
    dir.create(STATISTICS_FOLDER)
  }

  if(!dir.exists(STATISTICS_BEST_FOLDER)){
    dir.create(STATISTICS_BEST_FOLDER)
  }
  
  if(!dir.exists(STATISTICS_SAMPLE_WORDS_FOLDER)){
    dir.create(STATISTICS_SAMPLE_WORDS_FOLDER)
  }
  
  if(!dir.exists(STATISTICS_METRICS_FOLDER)){
    dir.create(STATISTICS_METRICS_FOLDER)
  }
  if(!dir.exists(STATISTICS_SAMPLE_METRICS_FOLDER)){
    dir.create(STATISTICS_SAMPLE_METRICS_FOLDER)
  }
  if(!dir.exists(STATISTICS_ANTIBIOTICS_METRICS_FOLDER)){
    dir.create(STATISTICS_ANTIBIOTICS_METRICS_FOLDER)
  }
  
  GENERATE_COMPARES()
  GENERATE_STATISTICS()

  GENERATE_SUB_STATISTICS()
  GENERATE_SPECIFIC_STATISTICS()
}






