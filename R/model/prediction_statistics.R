source('model/modelcommon.R')
library(dplyr)
library(stringr)


SUFFIXES <- c("strict","normal")
STRICT <- c(TRUE,FALSE)


readBigDt <- function(mode){
  output_dir <- getPredictionAltFolder(mode)
  outRDS <- file.path(output_dir,sprintf("modelOutput_long_merged.rds"))
  big_dt <- readRDS(outRDS)
}

calcMetric <- function(longDt){
  longDt %>%
    mutate(
      S = answer == "S",
      R = answer == "R",
      correctS  = S & pred == "S",
      correctR  = R & pred == "R",
      falseS = R & pred == "S", # note! other way around in notation for R and S
      falseR = S & pred == "R", # note! other way around in notation for R and S
      zerolabelS = S & pred == "", #answer S, pred empty
      zerolabelR = R & pred == "", #answer R, pred empty
      twolabelS = S & pred == "S/R",  #answer S, pred S/R
      twolabelR = R & pred == "S/R", #answer R, pred S/R
      total = R | S,          # counts all cases
    )
}

recreateAllstats <- function(bigDt,
                             by = c("word", "sample"),
                             strict = TRUE,
                             noinputab_keep = 4:8,
                             wordFilter = NULL) {
  by <- match.arg(by)
  
  # ---- Long + optional word filter -------------------------------------
  longDt <- bigDtToLong(bigDt, strict)
  
  if (!is.null(wordFilter)) {
    stopifnot(is.character(wordFilter))
    longDt <- longDt %>%
      dplyr::filter(.data$word %in% wordFilter)
  }
  
  # ---- Aggregate metrics ------------------------------------------------
  counts_fast <- longDt %>%
    calcMetric() %>%
    dplyr::group_by(
      significanceLevel, antibiotic, noinputab, .data[[by]]
    ) %>%
    dplyr::summarise(
      dplyr::across(METRICS, ~ sum(.x, na.rm = TRUE)),
      .groups = "drop"
    )
  
  # ---- Sanity checks ----------------------------------------------------
  sanity_check <- counts_fast %>%
    dplyr::mutate(
      ok_total = total == R + S,
      ok_R     = R == correctR + falseS + zerolabelR + twolabelR,
      ok_S     = S == correctS + falseR + zerolabelS + twolabelS
    )
  
  stopifnot(
    all(sanity_check %>%
          dplyr::summarise(
            dplyr::across(dplyr::starts_with("ok_"), all)
          ))
  )
  
  # ---- Long result ------------------------------------------------------
  result <- counts_fast %>%
    tidyr::pivot_longer(
      cols      = METRICS,
      names_to  = "metric",
      values_to = "count"
    ) %>%
    dplyr::filter(noinputab %in% noinputab_keep) %>%
    dplyr::select(
      noinputab,
      metric,
      antibiotic,
      significanceLevel,
      !!rlang::sym(by),
      count
    )
  
  result
}



# ---- Convenience wrappers (same names as before) ------------------------

recreateAllstatsWord <- function(bigDt, strict = TRUE, noinputab_keep = 4:8) {
  recreateAllstats(bigDt = bigDt, by = "word", strict = strict, noinputab_keep = noinputab_keep)
}

recreateAllstatsSample <- function(mode, strict = TRUE, noinputab_keep = 4:8) {
  recreateAllstats(bigDt = bigDt, by = "sample", strict = strict, noinputab_keep = noinputab_keep)
}

sumOnWords <- function(longDtWord)
{
  longDtWord %>%
    group_by(
      noinputab,
      metric,
      antibiotic,
      significanceLevel
    ) %>%
    summarise(
      count = sum(count, na.rm = TRUE),
      .groups = "drop"
    )
}
  


bigDtToLong <- function(pDt,strict = TRUE)
{
  columnToSignificanceLevel <- function(pSigni, pStrict){
    confLevelColName <- "pred"
    if(!is.na(pSigni) && !is.na(pSigni)){
      if(pStrict){
        confLevelColName <- sprintf("confpred_%s_strict",pSigni)
      } else {
        confLevelColName <- sprintf("confpred_%s",pSigni)
      }
    } 
    confLevelColName
  }

  columnToSignificanceLevelStrict <- function(pSigni)
  {
    columnToSignificanceLevel(pSigni,pStrict=TRUE)
  }
  columnToSignificanceLevelNonStrict <- function(pSigni)
  {
    columnToSignificanceLevel(pSigni,pStrict=FALSE)
  }
  
  if(strict){
    columns <- unlist(lapply(c(NA,signis),columnToSignificanceLevelStrict))
    columnsToRemove <- unlist(lapply(signis,columnToSignificanceLevelNonStrict))
  } else {
    columns <- unlist(lapply(c(NA,signis),columnToSignificanceLevelNonStrict))
    columnsToRemove <- unlist(lapply(signis,columnToSignificanceLevelStrict))
  }
  

  mapping <- c(NA, 0.025, 0.05, 0.1)
  names(mapping) <- columns
  
  pDt_long <- pDt %>%
    tidyr::pivot_longer(
      cols = all_of(columns),
      names_to = "column",
      values_to = "pred"
    ) %>%
    mutate(
      significanceLevel = mapping[column]
    ) %>%
    select(
     -column
    ) %>% 
    select(-all_of(columnsToRemove))
  
  pDt_long
}
  


subsetForStatistics <-function(pDt, pAntibiotic=NULL, pNoinputab=NULL, pWord=NULL,pSigni=NULL,strict = TRUE)
{
  # pDt <- readBigDt(MODE_A)
  # pAntibiotic <- "TOB"
  # pNoinputab <- 4
  # pSigni <- "01"
  # pWord <- "AMP_AMC_PIP_TZP"
  # strict <- TRUE
  
  confLevelColName <- "pred"
  if(!is.na(pSigni) && !is.na(pSigni)){
      if(strict){
        confLevelColName <- sprintf("confpred_%s_strict",pSigni)
      } else {
        confLevelColName <- sprintf("confpred_%s",signi)
      }
  } 
  
  frame <- pDt
  if(!is.null(pNoinputab)){
    frame <- frame %>%  filter(noinputab==pNoinputab) 
  }
  if(!is.null(pWord)){
    frame <- frame %>%  filter(word == pWord) 
  }
  if(!is.null(pAntibiotic)){
    frame <- frame %>%  filter(antibiotic == pAntibiotic) 
  }

  head(
    frame
  )
  frame <- frame %>% select(sample,word,antibiotic,answer,noinputab,pred = all_of(confLevelColName))
  frame
}


atLeastOneForEachGroup <- function(mode, index)
{
  shortFile <- paste("modelOutput_answer_",index,".csv",sep="")
  columns <- colnames(readOutputFrameShort(shortfile=shortFile,folder=getPredictionFolder(mode)))
  allGroupsInColumns <- unlist(lapply(columns,function(word){ all(unlist(lapply(AB_GROUPS,function(group){any(unlist(lapply(group, function(ab){grepl(pattern=ab,x = word)})))})))}))
  columns[allGroupsInColumns]
}

atLeastOneForEachGroupWords <- function()
{
  words <- list()
  for(index in 4:8){
    words <- append(words,atLeastOneForEachGroup(MODE_A,index))
  }
  unlist(words)
}




RUN <-function()
{
  oneperabgroup <- atLeastOneForEachGroupWords()

  for(mode in MODES){
    #mode <- MODE_A
    bigDt <- readBigDt(mode)
    
    for(i in 1:2) {
      strict <- STRICT[[i]]
      suffix <-   SUFFIXES[[i]]

      recreateAndDump <- function(filter,prefix){
        recreated_sample <- recreateAllstats(bigDt,by="sample",strict = strict,wordFilter = filter)
        recreated_word <- recreateAllstats(bigDt,by="word",strict = strict,wordFilter = filter)
        recreated_sum <- recreated_word %>% sumOnWords()
        writeStatisticsFrame(recreated_word, name = paste_non_na(prefix,"allStats_word",suffix,sep = "_"),getStatisticsFolder(mode))
        writeStatisticsFrame(recreated_sample, name = paste_non_na(prefix,"allStats_sample",suffix,sep = "_"),getStatisticsFolder(mode))
        writeStatisticsFrame(recreated_sum, name = paste_non_na(prefix,"allStats_sum",suffix,sep = "_"),getStatisticsFolder(mode))
      }
      
      recreateAndDump(filter = NULL,prefix = NA)
      recreateAndDump(filter = oneperabgroup,prefix = ONEPERABGROUP)
      for(name in names(BEST_SELECTION)){
        selection <- BEST_SELECTION[name]
        recreateAndDump(filter = selection,prefix <- name)
      }
    }
  }
}


MERGE <- function()
{
  if(!dir.exists(getCommonModelFolder())){
    dir.create(getCommonModelFolder())
  }
  mergeModes()
}




mergeModes <- function()
{
  prefixes <- c(NA,"oneperabgroup",names(BEST_SELECTION))
  names <- c("abErrorStatisticsCount","allStats_sum","allStats_sample","allStats_word")
  
  for(prefix in prefixes){
    print(prefix)
    for(name in names){
      print(name)
      #prefix <- NA
      #name <- "abErrorStatisticsCount"
      frameList <- lapply(SUFFIXES,FUN <- function(suffix){
        inFile <- paste_non_na(prefix,name,suffix,sep="_")
        print(inFile)
        outFile <- inFile
        frames <- lapply(MODES,FUN <- function(mode){readStatisticsFrame(name = inFile,getStatisticsFolder(mode))})
        names(frames) <- MODES
        merged <- bind_rows(frames, .id = "mode")
        #print(merged)
        #writeStatisticsFrame(merged,outFile,getCommonModelFolder())
        merged
      })
      names(frameList) <- SUFFIXES
      merged <- bind_rows(frameList, .id = "cpmode")
      outFile <- paste_non_na(prefix,name,sep="_")
      print(outFile)
      saveRDS(merged, file.path(getCommonModelFolder(),sprintf("%s.rds",outFile)))
      #writeStatisticsFrame(merged,outFile,getCommonModelFolder())
    }
  }
}


makeAbErrorStatisticCounts <- function()
{
  prefixes <- c(NA,"oneperabgroup",names(BEST_SELECTION))
  for(mode in MODES){
    print(mode)
    for(prefix in prefixes){
      for (suffix in SUFFIXES){
        inFile <- paste_non_na(prefix,sprintf("allStats_sum_%s",suffix),sep="_")
        outFile <- paste_non_na(prefix,sprintf("abErrorStatisticsCount_%s",suffix),sep="_")
        frame <- readStatisticsFrame(name = inFile,getStatisticsFolder(mode))
        wide <- frame %>% tidyr::pivot_wider(names_from = metric,values_from = count)
        writeStatisticsFrame(wide,outFile,getStatisticsFolder(mode))
      }
    }
  }
}

ALL <-function()
{
  RUN()
  makeAbErrorStatisticCounts()
  MERGE()
}