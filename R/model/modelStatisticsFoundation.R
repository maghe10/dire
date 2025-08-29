source(file = 'model/modelcommon.R')
library(stringr)

ROWS <- 99

RANGE <- 1:13


PENICILLINS <- c("AMP","AMC","PIP","TZP")
CEPHALOSPORINS <- c("CAZ","CRO","CTX","FEP")
FLOUROQUINOLONS <- c("CIP" ,"OFX" ,"LVX" ,"MFX")
AMINOGLYCOSIDES <- c("GEN" ,"TOB")
AB_GROUPS <- list(PENICILLINS,CEPHALOSPORINS,FLOUROQUINOLONS,AMINOGLYCOSIDES)


allWordsWithAtLeastOneForEachGroup <- function()
{
  unlist(lapply(4:13,atLeastoneForEachGroup))
}

atLeastOneForEachGroup <- function(index)
{
  shortFile <- paste("modelOutput_answer_",index,".csv",sep="")
  columns <- colnames(readOutputFrameShort(shortfile=shortFile,folder=getPredictionFolder(MODE)))
  allGroupsInColumns <- unlist(lapply(columns,function(word){ all(unlist(lapply(AB_GROUPS,function(group){any(unlist(lapply(group, function(ab){grepl(pattern=ab,x = word)})))})))}))
  columns[allGroupsInColumns]
}



predsfilesForIndex <-function(index)
{
  predsFiles <- makeShortFileName("preds",index)
  for(signi in signis){
    predsFiles <- c(predsFiles,makeShortFileName("confpreds",index,signi))
  }
  predsFiles
}

predsTocompareFile <- function(aFileName)
{
  stringr::str_replace(aFileName,"modelOutput","compare")
}


compareToStatisticsFile <- function(aFileName)
{
  stringr::str_replace(aFileName,"compare","statistics")
}

#Collapsed on word but keep sample
statisticsToSampleFile <- function(aFileName)
{
  stringr::str_replace(aFileName,"statistics","sample")
}

#Collapsed on sample but keep word
statisticsToWordFile <- function(aFileName)
{
  stringr::str_replace(aFileName,"statistics","word")
}

#Collapsed on sample and word
statisticsToSumFile <- function(aFileName)
{
  stringr::str_replace(aFileName,"statistics","statistics-sum")
}


comparefilesForIndex <-function(index)
{
  unlist(lapply(predsfilesForIndex(index),predsTocompareFile))
}

statisticsfilesForIndex <-function(index)
{
  unlist(lapply(comparefilesForIndex(index),compareToStatisticsFile))
}

samplefilesForIndex <-function(index)
{
  unlist(lapply(statisticsfilesForIndex(index),statisticsToSampleFile))
}

wordfilesForIndex <-function(index)
{
  unlist(lapply(statisticsfilesForIndex(index),statisticsToWordFile))
}

sumfilesForIndex <-function(index)
{
  unlist(lapply(statisticsfilesForIndex(index),statisticsToSumFile))
}


  
significanceLevelForStatisticsTable <- function()
{
  c(NA,unlist(lapply(signis,fixSigni)))
}

reorganizePredictions <- function(index)
{
  #index <- 1
  #
  #index <- 4
  #index <- 13

  #allAbs <- c("AMP","AMC","PIP","TZP","CAZ","CRO","CTX","FEP","CIP","OFX","LVX","MFX","GEN","TOB")
  #correctFrame <- readOutputFrameShort(shortfile="modelOutput_answer_13.csv",folder=getPredictionFolder(MODE))
  #colnames(correctFrame) <- rev(allAbs)
  
  predsFiles <- predsfilesForIndex(index)
  answerFile <- makeShortFileName("answer",index)

  predictionFrames <- lapply(predsFiles,function(x){ 
    aFrame <- readOutputFrameShort(shortfile=x,folder=getPredictionFolder(MODE))
    aFrame
  })
  names(predictionFrames) <- predsFiles
  
  answerFrame <- readOutputFrameShort(shortfile=answerFile,folder=getPredictionFolder(MODE))

  reorganizeFrame <- function(frame,tag)
  {
    words <- colnames(frame)
    frame[frame=="<empty>"] <-NA
    frameLong <- cbind(frame,sample=rownames(frame)) %>% 
      tidyr::pivot_longer(cols = all_of(words), names_to = "word", values_to = tag, values_drop_na = TRUE) %>% 
      tidyr::separate_longer_delim(cols = all_of(tag),delim = " ") %>% 
      tidyr::separate_wider_delim(cols = all_of(tag),delim = "_",names = c("antibiotic",paste(tag,"sr",sep = "_")))
    
    frameLong
  }
  
  answerFrameLong <-  reorganizeFrame(answerFrame,"answer")
  for(aFileName in names(predictionFrames)){
    #aFileName <- names(predictionFrames)[[4]]
    predictionFrame <- predictionFrames[[aFileName]]
    predictionFrameLong <- reorganizeFrame(predictionFrame,"prediction")
    combinedFrameLong <- dplyr::right_join(predictionFrameLong, answerFrameLong, by = c("sample", "word", "antibiotic"))
    
    
    outFileName <- predsTocompareFile(aFileName)
    outFileNameFullGz <-paste(getCompareFolder(MODE),"/",outFileName,".gz",sep="")
    data.table::fwrite(combinedFrameLong,outFileNameFullGz)
    #writeStatisticsFrame(combinedFrameLong,name = outFileName,folder = getCompareFolder(MODE))
  }
}


# COLAPSE_ON_SAMPLE <- function()
# {
#   inFileName <- "allStats_widest.csv"
#   inFileNameFullGz <-paste(getStatisticsFolder(MODE),"/",inFileName,".gz",sep="")
#   full <- data.table::fread(inFileNameFullGz)
#   
#   sum <-  full %>% 
#     dplyr::group_by(noinputab,word,metric,antibiotic,significanceLevel) %>%
#     dplyr::summarize(count=sum(count),.groups = "drop")
# 
#   
#   outFileName <- "allStats_wider.csv"
#   outFileNameFullGz <-paste(getStatisticsFolder(MODE),"/",outFileName,".gz",sep="")
#   data.table::fwrite(sum,outFileNameFullGz)
# }

# COLAPSE_ON_WORD <- function()
# {
#   inFileName <- "allStats_widest.csv"
#   inFileNameFullGz <-paste(getStatisticsFolder(MODE),"/",inFileName,".gz",sep="")
#   full <- data.table::fread(inFileNameFullGz)
#   
#   sum <-  full %>% 
#     dplyr::group_by(across(c(-word,-sample))) %>%
#     dplyr::summarize(count=sum(count),.groups = "drop")
#   
#   
#   outFileName <- "allStats_sample.csv"
#   outFileNameFullGz <-paste(getStatisticsFolder(MODE),"/",outFileName,".gz",sep="")
#   data.table::fwrite(sum,outFileNameFullGz)
# }


collapseOnWord <-function(table)
{
  table %>% 
    dplyr::group_by(across(c(-word,-count))) %>%
    dplyr::summarize(count=sum(count),.groups = "drop")
}

collapseOnSample <-function(table)
{
  table %>% 
    dplyr::group_by(across(c(-sample,-count))) %>%
    dplyr::summarize(count=sum(count),.groups = "drop")
}


MERGE_SUM <- function(prefix=NA,range=RANGE)
{
  allFiles <- unlist(lapply(range,function(index){
    statisticsFiles <- sumfilesForIndex(index)
    if(!is.na(prefix)){
      statisticsFiles <- paste(prefix,statisticsFiles,sep="_")
    }
    statisticsFiles
  }))
  file_list <- paste(paste(getStatisticsTmpFolder(MODE),"/",allFiles,".gz",sep=""))

  mergedFrame <- data.table::rbindlist(lapply(file_list,data.table::fread))
  outFileName <- "allStats_sum"
  if(!is.na(prefix)){
    outFileName <- paste(prefix,outFileName,sep="_")
  }
  #outFileNameFullGz <-paste(getStatisticsFolder(MODE),"/",outFileName,".gz",sep="")
  writeStatisticsFrame(mergedFrame,name=outFileName,getStatisticsFolder(MODE))
}

MERGE_WORD <- function(prefix=NA,range=RANGE)
{
  allFiles <- unlist(lapply(range,function(index){
    statisticsFiles <- wordfilesForIndex(index)
    if(!is.na(prefix)){
      statisticsFiles <- paste(prefix,statisticsFiles,sep="_")
    }
    statisticsFiles
  }))
  file_list <- paste(paste(getStatisticsTmpFolder(MODE),"/",allFiles,".gz",sep=""))
  
  mergedFrame <- data.table::rbindlist(lapply(file_list,data.table::fread))
  nrow(mergedFrame)
  outFileName <- "allStats_word"
  if(!is.na(prefix)){
    outFileName <- paste(prefix,outFileName,sep="_")
  }
  #outFileNameFullGz <-paste(getStatisticsFolder(MODE),"/",outFileName,".gz",sep="")
  writeStatisticsFrame(mergedFrame,name=outFileName,getStatisticsFolder(MODE))
}

MERGE_SAMPLE <- function(prefix=NA,range=RANGE)
{
  allFiles <- unlist(lapply(range,function(index){
    statisticsFiles <- samplefilesForIndex(index)
    if(!is.na(prefix)){
      statisticsFiles <- paste(prefix,statisticsFiles,sep="_")
    }
    statisticsFiles
  }))
  file_list <- paste(paste(getStatisticsTmpFolder(MODE),"/",allFiles,".gz",sep=""))
  
  mergedFrame <- data.table::rbindlist(lapply(file_list,data.table::fread))
  nrow(mergedFrame)
  outFileName <- "allStats_sample"
  if(!is.na(prefix)){
    outFileName <- paste(prefix,outFileName,sep="_")
  }
  #outFileNameFullGz <-paste(getStatisticsFolder(MODE),"/",outFileName,".gz",sep="")
  writeStatisticsFrame(mergedFrame,name=outFileName,getStatisticsFolder(MODE))
}

# MERGE_STATISTICS_METRICS <- function()
# {
# 
#   allFiles <- unlist(lapply(RANGE,function(index){
#     statisticsFiles <- statisticsfilesForIndex(index)
#     statisticsFiles
#   }))
#   
#   file_list <- paste(paste(getStatisticsTmpFolder(MODE),"/",allFiles,".gz",sep=""))
#   outFileName <- "allStats_widest.csv"
#   outFileNameFullGz <-paste(getStatisticsFolder(MODE),"/",outFileName,".gz",sep="")
#   
#   for (i in seq_along(file_list)) {
#     dt <- data.table::fread(file_list[i])
#     print(i)
#     # data.table::fwrite(
#     #   dt,
#     #   file = outFileNameFullGz,
#     #   append = (i != 1),  # Only write header once (for the first file)
#     #   compress = "gzip"
#     #)
#   }
#   
# }

GENERATE_STATISTICS_METRICS <- function()
{

  for(index in RANGE){
    #MODE <- "Mode-A"
    #index <- 2
    compareFiles <- comparefilesForIndex(index)
    names(compareFiles) <- c("crude",signis)
    compareFrames <- lapply(names(compareFiles),function(x){ 
      #x<-"01"
      compareFile <- compareFiles[[x]]
      aFrame <- data.table::fread(paste(getCompareFolder(MODE),"/",compareFile,".gz",sep=""))
      significanceLevel <- NA
      if(x!="crude"){
        significanceLevel <- fixSigni(x)
      }
      # print(head(aFrame))
      # print(tail(aFrame))
      statisticsFrame <- comparesToStatisticsFrame(aFrame,index,significanceLevel)
      outFileName <- compareToStatisticsFile(compareFile)
      outFileNameFullGz <-paste(getStatisticsTmpFolder(MODE),"/",outFileName,".gz",sep="")
      data.table::fwrite(statisticsFrame,outFileNameFullGz)
    })
  }
}


GENERATE_STATISTICS_SUB_METRICS <- function()
{
  
  for(index in RANGE){
    #MODE <- "Mode-A"
    #index <- 4
    statisticsfiles <- statisticsfilesForIndex(index)
    names(statisticsfiles) <- c("crude",signis)
    lapply(names(statisticsfiles),function(x){ 
      #x<-"01"
      file <- statisticsfiles[[x]]
      aFrame <- data.table::fread(paste(getStatisticsTmpFolder(MODE),"/",file,".gz",sep=""))
      significanceLevel <- NA
      if(x!="crude"){
        significanceLevel <- fixSigni(x)
      }
      #Fix correct order
      aFrame <- aFrame %>% dplyr::select(noinputab,metric,count,antibiotic,significanceLevel,word,sample)
      # print(head(aFrame))
      # print(tail(aFrame))
      # print(nrow(aFrame))
      
      
      word <- collapseOnSample(aFrame)
      # print(nrow(word))
      sample <- collapseOnWord(aFrame)
      # print(nrow(sample))
      sum <- collapseOnWord(word)
      # print(nrow(sum))
      if(index>=4){
          words <- atLeastOneForEachGroup(index)
          aFrame_oneperabgroup <- aFrame %>% dplyr::filter(word %in% words)
          
          word_oneperabgroup <-  collapseOnSample(aFrame_oneperabgroup)
          sample_oneperabgroup <-  collapseOnWord(aFrame_oneperabgroup)
          sum_oneperabgroup <-  collapseOnWord(word_oneperabgroup)
      }
      
      outFileName <- statisticsToWordFile(file)
      outFileNameFullGz <-paste(getStatisticsTmpFolder(MODE),"/",outFileName,".gz",sep="")
      data.table::fwrite(x=word,outFileNameFullGz)
      
      outFileName <- statisticsToSampleFile(file)
      outFileNameFullGz <-paste(getStatisticsTmpFolder(MODE),"/",outFileName,".gz",sep="")
      data.table::fwrite(x=sample,outFileNameFullGz)

      outFileName <- statisticsToSumFile(file)
      outFileNameFullGz <-paste(getStatisticsTmpFolder(MODE),"/",outFileName,".gz",sep="")
      data.table::fwrite(x=sum,outFileNameFullGz)
      if(index>=4){
        outFileName <- paste("oneperabgroup",statisticsToWordFile(file),sep="_")
        outFileNameFullGz <-paste(getStatisticsTmpFolder(MODE),"/",outFileName,".gz",sep="")
        data.table::fwrite(x=word_oneperabgroup,outFileNameFullGz)
        
        outFileName <- paste("oneperabgroup",statisticsToSampleFile(file),sep="_")
        outFileNameFullGz <-paste(getStatisticsTmpFolder(MODE),"/",outFileName,".gz",sep="")
        data.table::fwrite(x=sample_oneperabgroup,outFileNameFullGz)
        
        outFileName <- paste("oneperabgroup",statisticsToSumFile(file),sep="_")
        outFileNameFullGz <-paste(getStatisticsTmpFolder(MODE),"/",outFileName,".gz",sep="")
        data.table::fwrite(x=sum_oneperabgroup,outFileNameFullGz)
      }
      
    })
  }
}



comparesToStatisticsFrame <- function(compareFrame,index,significanceLevel)
{
  #index <- 1
  #significanceLevel <- NA
  head(compareFrame)
  #compareFrame$prediction_sr!="SR" & compareFrame$prediction_sr!=compareFrame$answer
  #compareFrame$answer_sr
  
  wideMetricsFrame <- compareFrame %>% 
    dplyr::mutate(correct = prediction_sr==answer_sr, false = (prediction_sr=="R" & answer_sr=="S")|(prediction_sr=="S" & answer_sr=="R") , ambiguous=prediction_sr=="SR", notpredicted=prediction_sr=="",S=answer_sr=="S",R=answer_sr=="R") %>% 
    dplyr::mutate(correctS = S & correct, correctR = R & correct, falseR = S & false,falseS = R & false, notpredictedS = S & notpredicted, notpredictedR = R & notpredicted, ambiguousS = S & ambiguous,ambiguousR = R & ambiguous, total = S | R, predicted = correct | false, predorambiguous = correct | false | ambiguous) %>%
    dplyr::select(-prediction_sr,-answer_sr)
  stopifnot(all(wideMetricsFrame$total))
  wideMetricsFrame %>% 
    tidyr::pivot_longer(cols <- all_of(colnames(wideMetricsFrame)[4:ncol(wideMetricsFrame)]),names_to = "metric",values_to = "count") %>% 
    dplyr::mutate(count = as.integer(count), noinputab=index,significanceLevel=significanceLevel) %>%
    dplyr::select(noinputab,metric,count,antibiotic,significanceLevel,word,sample)
}


GENERATE_REORGANIZED_PREDICTIONS <- function()
{
  for(i in RANGE){
    reorganizePredictions(i)
  }
}

ALL <- function()
{
  checkDirs()
  # make each prediction file into a tidy form with
  GENERATE_REORGANIZED_PREDICTIONS()
  GENERATE_STATISTICS_METRICS()
  MOST()
}

MOST <- function()
{
  checkDirs()
  
  GENERATE_STATISTICS_SUB_METRICS()
  MERGE_SUM()
  MERGE_WORD()
  MERGE_SAMPLE()
  #SANITY_CHECK_SUM()
  #SANITY_CHECK_WORD()
  MERGE_SUM(prefix = "oneperabgroup",range = 4:max(RANGE))
  MERGE_WORD(prefix = "oneperabgroup",range = 4:max(RANGE))
  MERGE_SAMPLE(prefix = "oneperabgroup",range = 4:max(RANGE))
}





# GENERATE_STATISTICS_ONE_PER_AB_GROUP <- function()
# {
#   altWider <- readStatisticsFrame(name = "allstats_word",getStatisticsFolder(MODE))
#   
#   words <- allWordsWithAtLeastOneForEachGroup()
#   allStats <- altWider %>% dplyr::filter(word %in% words)
#   writeStatisticsFrame(allStats, name = "allstats_alt_wider_oneperabgroup",getStatisticsFolder(MODE))
#   
#   sum_wider <-  allStats %>% 
#     dplyr::group_by(across(c(-word,-count))) %>%
#     dplyr::summarize(count=sum(count),.groups = "drop")
#   
#   sum_wider <- sum_wider %>% dplyr::select(noinputab, metric, count, antibiotic, significanceLevel)
#   sum_wider <- sum_wider  %>% dplyr::arrange_all()
#   
#   # alt <-readStatisticsFrame("allStats_sum",getStatisticsFolder(MODE))
#   # alt <- alt %>% dplyr::arrange_all()
#   # 
#   # all(alt==sum_wider,na.rm=TRUE)
#   # any(alt!=sum_wider,na.rm=TRUE)
#   
#   writeStatisticsFrame(sum_wider, name = "allstats_alt_oneperabgroup",getStatisticsFolder(MODE))
# }




SANITY_CHECK_SUM <- function()
{
  alt <-readStatisticsFrame("allstats_alt",getStatisticsFolder(MODE))
  altSorted <- alt %>% dplyr::select(noinputab,metric,antibiotic, significanceLevel=significanslevel,count) %>% dplyr::mutate_at(dplyr::vars(significanceLevel), as.numeric)  %>% dplyr::arrange(across(everything(), desc))
  altSorted <- altSorted %>% dplyr::filter(noinputab %in% RANGE)
  
  sum <-readStatisticsFrame("allstats_sum",getStatisticsFolder(MODE))
  sumSorted <- sum %>% dplyr::select(noinputab,metric, antibiotic, significanceLevel,count) %>% dplyr::arrange(across(everything(), desc))


    
  all(altSorted$noinputab==sumSorted$noinputab)
  all(altSorted$metric==sumSorted$metric)
  
  altdiffCount <- altSorted[which(altSorted$count!=sumSorted$count),]
  head(altdiffCount)
  sumdiffCount <- sumSorted[which(altSorted$count!=sumSorted$count),]
  head(sumdiffCount)
  
  stopifnot(!any(altSorted!=sumSorted,na.rm=TRUE))
  stopifnot(all(altSorted==sumSorted,na.rm=TRUE))
}
  
  

SANITY_CHECK_WORD <- function()
{
  altWider <-readStatisticsFrame("allstats_alt_wider",getStatisticsFolder(MODE))
  unique(altWider$noinputab)
  altWiderSorted <- altWider %>% dplyr::select(noinputab = noinputab ,metric = metric, antibiotic=antibiotic, significanceLevel=significanslevel,word=word, count=count) %>% dplyr::arrange_all()
  altWiderSorted <- altWiderSorted %>% dplyr::filter(noinputab %in% RANGE)

  toKeep <- unlist(lapply(1:nrow(altWiderSorted),function(i){
    !grepl(altWiderSorted[i,"antibiotic"],altWiderSorted[i,"word"])
  }))
  altWiderSorted <- altWiderSorted[toKeep,]
  altWiderSorted[altWiderSorted=="NA"] <- NA
  altWiderSorted$significanceLevel <- as.numeric(altWiderSorted$significanceLevel)
  
  
  wider <-readStatisticsFrame("allstats_word",getStatisticsFolder(MODE))
  widerSorted <- wider %>% dplyr::select(noinputab = noinputab ,metric = metric,antibiotic=antibiotic, significanceLevel=significanceLevel,word=word, count=count) %>% dplyr::arrange_all()

  # for debugging
  # differ <- which(altWiderSorted$count !=widerSorted$count)
  # unique(altWiderSorted[differ,2])
  # head(differ)
  # altWiderSorted[differ,]
  # widerSorted[differ,]
  # colnames(altWiderSorted)
  # 
  # widerSorted[differ & widerSorted$antibiotic=="CAZ" & widerSorted$word=="AMC",]
  

  #FIXME format
  #stopifnot(all(sort(unique(altWider$significanslevel))==sort(unique(wider$significanceLevel))))
  stopifnot(nrow(altWiderSorted)==nrow(widerSorted))
  
    

  stopifnot(all(altWiderSorted==widerSorted,na.rm=TRUE))
  stopifnot(!any(altWiderSorted!=widerSorted,na.rm=TRUE))
  
  # all(alt==sum_wider,na.rm=TRUE)
  # any(alt!=sum_wider,na.rm=TRUE)
}

checkDirs <-function()
{
  stopifnot(dir.exists(getPredictionBase(MODE)))
  
  
  if(!dir.exists(getStatisticsFolder(MODE))){
    dir.create(getStatisticsFolder(MODE))
  }
  if(!dir.exists(getCompareFolder(MODE))){
    dir.create(getCompareFolder(MODE))
  }
  
  if(!dir.exists(getStatisticsTmpFolder(MODE))){
    dir.create(getStatisticsTmpFolder(MODE))
  }
  
}



