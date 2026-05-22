source('model/modelcommon.R')

readAllstatsWord <- function(mode){
  output_dir <- getStatisticsFolder(mode)
  outAllstatsWord <- file.path(output_dir,sprintf("allStats_word.csv"))
  allstatsWord <- read.csv2(outAllstatsWord)
}


compateWithASW <- function(mode = MODE_A)
{
  asw <- readAllstatsWord(mode)
  recreatedAsw <- recreateAllstatsWord(mode,strict = FALSE)
  CORE_METRICS <- c("total","R","S","falseR","falseS","correctS","correctR")
  aswLimited <- asw %>% filter(metric %in% CORE_METRICS)
  recreatedAswLimited <- recreatedAsw %>% filter(metric %in% CORE_METRICS)
  stopifnot(nrow(recreatedAswLimited)==nrow(aswLimited))
  
  keys <- c("noinputab","metric","antibiotic","significanceLevel","word","count")
  rec_sorted <- recreatedAswLimited %>%
    select(all_of(keys)) %>%
    arrange(noinputab, metric, antibiotic, significanceLevel, word, count)
  asw_sorted <- aswLimited %>%
    select(all_of(keys)) %>%
    arrange(noinputab, metric, antibiotic, significanceLevel, word, count)
  rec_sorted <- as_tibble(rec_sorted)
  asw_sorted <- as_tibble(asw_sorted)
  typeof(rec_sorted$significanceLevel)
  typeof(asw_sorted$significanceLevel)
  
  print(n=100,anti_join(rec_sorted, asw_sorted, by = keys))
  print(n=100,anti_join(asw_sorted, rec_sorted, by = keys))
  
}

COMPARE <- function()
{
  print("Mode A")
  compateWithASW(mode=MODE_A)
  print("Mode B")
  compateWithASW(mode=MODE_B)
  print("Mode C")
  compateWithASW(mode=MODE_C)
}



