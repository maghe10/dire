MODE <- MODE_A
source("model/modelcommon.R")
library(dplyr)

sanitycheckCompares <- function()
{
  compareFile <- "compare_preds_6.csv"
  aFrame <- data.table::fread(paste(getCompareFolder(MODE),"/",compareFile,".gz",sep=""))
  
  stopifnot(sum(aFrame$prediction_sr == "R") + sum(aFrame$prediction_sr == "S") 
         == 
           sum(aFrame$answer_sr == "R") + sum(aFrame$answer_sr == "S")  )
  
  if(sum(aFrame$prediction_sr == "R") > sum(aFrame$answer_sr == "R"))
  {
    print("most ME")
  } else {
    print("most VME")
  }
}


sanitycheckPreds <- function()
{
  answerFile <- "modelOutput_answer_6.csv"
  predsFile <- "modelOutput_preds_6.csv"
  answerFrame <- data.table::fread(paste(getPredictionFolder(MODE),"/",answerFile,sep="")) 
  predsFrame <- data.table::fread(paste(getPredictionFolder(MODE),"/",predsFile,sep="")) 

  total_answer_S <- answerFrame %>% select(-1) %>% summarise(across(everything(), ~ sum(str_count(.x, "_S\\b"), na.rm = TRUE))) %>% rowSums()
  total_answer_R <- answerFrame %>% select(-1) %>% summarise(across(everything(), ~ sum(str_count(.x, "_R\\b"), na.rm = TRUE))) %>% rowSums()

  total_preds_S <- predsFrame %>% select(-1) %>% summarise(across(everything(), ~ sum(str_count(.x, "_S\\b"), na.rm = TRUE))) %>% rowSums()
  total_preds_R <- predsFrame %>% select(-1) %>% summarise(across(everything(), ~ sum(str_count(.x, "_R\\b"), na.rm = TRUE))) %>% rowSums()
  
  stopifnot(total_preds_S + total_preds_R 
            == 
              total_answer_S+ total_answer_R )
  
  if(total_preds_R > total_answer_R)
  {
    print("most ME")
  } else {
    print("most VME")
  }
}

sanitycheckJuanDirectly <- function(){
 file <- "C:\\Users\\magnu\\OneDrive - Västra Götalandsregionen\\git\\dire\\model\\Confidence-based-Prediction-of-Antibiotic-Resistance\\data\\output_Mode-A.csv"
 frame <- read.csv(file)
 sum(frame$AST_prediction=="R" | frame$AST_prediction=="S")
 if(sum(frame$AST_true=="R" & frame$AST_prediction=="S") < sum(frame$AST_true=="S" & frame$AST_prediction=="R")){
   print("most ME")
} else {
  print("most VME")
}

}


