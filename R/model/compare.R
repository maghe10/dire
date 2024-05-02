source(file = 'model/modelcommon.R')


createComparison <- function(correctCount, predictionsCount)
{
  doubleCorrect <-
    apply(correctCount, c(1, 2), function(x) {
      as.double(x)
    })
  
  resultFrame = doubleCorrect / predictionsCount
  
  resultFrame
}

checkAmbiguous <- function (word,correct) {
  length(grep("_SR", word)) > 0
}

countAmbiguous <- function (answerDataFrame, predictionDataFrame,ab=NULL) {
  genericCount(answerDataFrame,
               predictionDataFrame,ab,checkAmbiguous) 
}

checkIncorrect <- function(word,correct,vme)
{
  match <- FALSE
  if (!length(grep("_SR", word)) > 0) {
    # Not ending with _SR => compare with _R or _S
    if (length(grep("_R", word)) > 0) {
      opposite = paste(substr(word, 1, 4), "S", sep = "")
    }
    else {
      opposite = paste(substr(word, 1, 4), "R", sep = "")
    }
    if (length(grep(opposite, correct)) > 0) {
      predR = length(grep("_R", word)) > 0
      if (xor(x = vme, y = predR)) {
        match <- TRUE
      }
    }
  }
  match  
}

countIncorrect <-
  function (answerDataFrame,
            predictionDataFrame,
            vme = FALSE, ab=NULL)
  {
    genericCount(answerDataFrame,
                 predictionDataFrame,ab,checkIncorrect,vme)
  }

checkCorrect <- function(word,correct)
{
  length(grep(word, correct) > 0)
}

countCorrect <-
  function (answerDataFrame,
            predictionDataFrame, ab=NULL)
  {
    genericCount(answerDataFrame,
                 predictionDataFrame,ab,checkCorrect)
  }

checkPredicted <- function(word,correct,technical)
{
  answer = FALSE
  if (!technical) {
    ## Count only those that are unambiguous 
    answer <- length(grep("_SR", word, invert = TRUE))>0
  }
  else {
    ## Count all predictions 
    answer = TRUE 
  }
  answer
}

countPredicted <-
  function (answerDataFrame,
            predictionDataFrame,technical = FALSE, ab=NULL)
  {
    genericCount(answerDataFrame,
                 predictionDataFrame,ab,checkPredicted,technical)
  }


STATISTICS  <- function()
{
  statistics(TRUE, FALSE)
  statistics(FALSE, FALSE)
  statistics(TRUE, TRUE)
  statistics(FALSE, TRUE)
}

statistics <- function(bySample, technical)
{
  if (bySample) {
    applyIndex = 1
  } else {
    applyIndex = 2
  }
  
  resultFrame <- NULL
  
  for (index in range) {
    correct <-
      readOutputFrame(
        k = index,
        type = "correct",
        folder = paste(modelDirectory, "temp", sep = "/")
      )
    vme <-
      readOutputFrame(
        k = index,
        type = "vme",
        folder = paste(modelDirectory, "temp", sep = "/")
      )
    me <-
      readOutputFrame(
        k = index,
        type = "me",
        folder = paste(modelDirectory, "temp", sep = "/")
      )
    predictions <-
      readOutputFrame(
        k = index,
        type = "predictions",
        folder = paste(modelDirectory, "temp", sep = "/")
      )
    correctSum <- apply(correct, applyIndex, sum)
    predictionsSum <- apply(predictions, applyIndex, sum)
    vmeSum <- apply(vme, applyIndex, sum)
    meSum <- apply(me, applyIndex, sum)
    percentage <- correctSum / predictionsSum
    if (bySample) {
      if (is.null(resultFrame)) {
        resultFrame <- data.frame(percentage)
      } else {
        resultFrame <- cbind(resultFrame, percentage)
      }
    } else {
      resultFrame <- data.frame(percentage)
    }
    colnames(resultFrame)[ncol(resultFrame)] <-
      paste("percentage_preds", index, sep = "_")
    resultFrame <- cbind(resultFrame, correctSum, predictionsSum)
    colnames(resultFrame)[ncol(resultFrame) - 1] <-
      paste("correct_preds", index, sep = "_")
    colnames(resultFrame)[ncol(resultFrame)] <-
      paste("all_preds", index, sep = "")
    resultFrame <- cbind(resultFrame, vmeSum, meSum)
    colnames(resultFrame)[ncol(resultFrame) - 1] <-
      paste("vme", index, sep = "_")
    colnames(resultFrame)[ncol(resultFrame)] <-
      paste("me", index, sep = "")
    
    
    best = percentage[which(percentage == max(percentage))]
    worst = percentage[which(percentage == min(percentage))]
    #    print(paste("Index", index, "Mean", mean(percentage), "Max",max(percentage),"Min", min(percentage),"Best",names(best),"Worst",names(worst),sep=" "))
    for (signi in signis) {
      correct_type <-  "correct_conformal"
      if (technical) {
        prediction_type <- "predictions_conformal_technical"
      } else {
        prediction_type <- "predictions_conformal_clinical"
      }
      correct <-
        readOutputFrame(
          k = index,
          type = correct_type,
          signi = signi,
          folder = paste(modelDirectory, "temp", sep = "/")
        )
      ambiguous <-
        readOutputFrame(
          k = index,
          type = "ambiguous_conformal",
          signi = signi,
          folder = paste(modelDirectory, "temp", sep = "/")
        )
      predictions <-
        readOutputFrame(
          k = index,
          type = prediction_type,
          signi = signi,
          folder = paste(modelDirectory, "temp", sep = "/")
        )
      correctSum <- apply(correct, applyIndex, sum)
      if(technical){
        correctSum <- correctSum + apply(ambiguous, applyIndex, sum)
      }
      predictionsSum <- apply(predictions, applyIndex, sum)
      percentage <- correctSum / predictionsSum
      
      vme_conformal <-
        readOutputFrame(
          k = index,
          type = "vme_conformal",
          signi = signi,
          folder = paste(modelDirectory, "temp", sep = "/")
        )
      me_conformal <-
        readOutputFrame(
          k = index,
          type = "me_conformal",
          signi = signi,
          folder = paste(modelDirectory, "temp", sep = "/")
        )
      vme_conformal_sum <- apply(vme_conformal, applyIndex, sum)
      me_conformal_sum <- apply(me_conformal, applyIndex, sum)
      
      
      resultFrame <-
        cbind(resultFrame, percentage, correctSum, predictionsSum)
      colnames(resultFrame)[ncol(resultFrame) - 2] <-
        paste(fixSigni(signi), "percentage", index, sep = "_")
      colnames(resultFrame)[ncol(resultFrame) - 1] <-
        paste(fixSigni(signi), "correct", index, sep = "_")
      colnames(resultFrame)[ncol(resultFrame)] <-
        paste(fixSigni(signi), "all", index, sep = "_")
      resultFrame <-
        cbind(resultFrame, vme_conformal_sum, me_conformal_sum)
      colnames(resultFrame)[ncol(resultFrame) - 1] <-
        paste(fixSigni(signi), "vme", index, sep = "_")
      colnames(resultFrame)[ncol(resultFrame)] <-
        paste(fixSigni(signi), "me", index, sep = "_")
      
      
      best = percentage[which(percentage == max(percentage))]
      worst = percentage[which(percentage == min(percentage))]
      #      print(paste("Index", index, "Signi",signi, "Mean", mean(percentage), "Max",max(percentage),"Min", min(percentage),"Best",names(best),"Worst",names(worst),sep=" "))
      if (!bySample) {
        #        print(resultFrame)
        if (technical) {
          writeOutputFrame(resultFrame, k = index, type = "stats_comb_technical")
        } else {
          writeOutputFrame(resultFrame, k = index, type = "stats_comb_clinical")
        }
        
      }
    }
    
  }
  if (bySample) {
    #    print(resultFrame)
    if (technical) {
      writeOutputFrame(resultFrame, k = "all", type = "stats_samp_technical")
    }
    else {
      writeOutputFrame(resultFrame, k = "all", type = "stats_samp_clinical")
    }
  }
}




COUNT <- function()
{
  for (index in range) {
    inFrame <- readInputframe(index)
    answerFrame <- readOutputFrame(k = index, type = "answer")
    
    for (signi in signis) {
      conformalPredictionFrame <-
        readOutputFrame(k = index, signi, type = "confpreds")
      vme <-
        countIncorrect(answerFrame, conformalPredictionFrame, TRUE)
      writeOutputFrame(vme, k = index, signi, type = "vme_conformal")
      me <-
        countIncorrect(answerFrame, conformalPredictionFrame, FALSE)
      writeOutputFrame(me, k = index, signi, type = "me_conformal")
      
      
      correct <-
        countCorrect(answerFrame, conformalPredictionFrame)
      ambiguous <-
        countAmbiguous(answerFrame, conformalPredictionFrame)
      predictionsTechnical <- countPredicted(answerFrame,conformalPredictionFrame,TRUE)

      
      comparisonTechnical <- createComparison(correct+ambiguous, predictionsTechnical)
      writeOutputFrame(comparisonTechnical, k = index, signi, type = "comparison_conformal_technical")
      writeOutputFrame(predictionsTechnical, k = index, signi, type = "predictions_conformal_technical")
      writeOutputFrame(correct, k = index, signi, type = "correct_conformal")
      writeOutputFrame(ambiguous, k = index, signi, type = "ambiguous_conformal")
      
            
      predictionsClinical <- countPredicted(answerFrame,conformalPredictionFrame, FALSE)
      comparisonClinical <- createComparison(correct, predictionsClinical)
      writeOutputFrame(comparisonClinical, k = index, signi, type = "comparison_conformal_clinical")
      writeOutputFrame(predictionsClinical, k = index, signi, type = "predictions_conformal_clinical")
      
    }
    predsFrame <- readOutputFrame(k = index, type = "preds")
    correct <- countCorrect(answerFrame, predsFrame)
    predictions <- countPredicted(answerFrame,predsFrame)
    comparison <- createComparison(correct, predictions)
    writeOutputFrame(comparison, k = index, type = "comparison")
    writeOutputFrame(predictions, k = index, type = "predictions")
    writeOutputFrame(correct, k = index, type = "correct")
    
    vme <- countIncorrect(answerFrame, predsFrame, TRUE)
    writeOutputFrame(vme, k = index, type = "vme")
    me <- countIncorrect(answerFrame, predsFrame, FALSE)
    writeOutputFrame(me, k = index, type = "me")
  }
}


COUNT_WITH_AB <- function(ab=NULL)
{
  for (index in range) {
    inFrame <- readInputframe(index)
    answerFrame <- readOutputFrame(k = index, type = "answer")
    
    for (signi in signis) {
      conformalPredictionFrame <-
        readOutputFrame(k = index, signi, type = "confpreds")
      vme <-
        countIncorrect(answerFrame, conformalPredictionFrame,ab,TRUE)
      writeOutputFrame(vme, k = index, signi, type = paste("vme_conformal",ab,sep = "_"))
      me <-
        countIncorrect(answerFrame, conformalPredictionFrame,ab,FALSE)
      writeOutputFrame(me, k = index, signi, type = paste("me_conformal",ab,sep = "_"))
      
      
      correct <-
        countCorrect(answerFrame, conformalPredictionFrame,ab)
      ambiguous <-
        countAmbiguous(answerFrame, conformalPredictionFrame,ab)
      predictionsTechnical <- countPredicted(answerFrame,conformalPredictionFrame,TRUE,ab)
      
      
      comparisonTechnical <- createComparison(correct+ambiguous, predictionsTechnical)
      writeOutputFrame(comparisonTechnical, k = index, signi, type = paste("comparison_conformal_technical",ab,sep = "_"))
      writeOutputFrame(predictionsTechnical, k = index, signi, type = paste("predictions_conformal_technical",ab,sep = "_"))
      writeOutputFrame(correct, k = index, signi, type = paste("correct_conformal",ab,sep = "_") )
      writeOutputFrame(ambiguous, k = index, signi, type = paste("ambiguous_conformal",ab,sep = "_") )
      
      
      predictionsClinical <- countPredicted(answerFrame,conformalPredictionFrame, FALSE,ab)
      comparisonClinical <- createComparison(correct, predictionsClinical)
      writeOutputFrame(comparisonClinical, k = index, signi, type = paste("comparison_conformal_clinical",ab,sep = "_"))
      writeOutputFrame(predictionsClinical, k = index, signi, type = paste("predictions_conformal_clinical",ab,sep = "_"))
      
    }
    predsFrame <- readOutputFrame(k = index, type = "preds")
    correct <- countCorrect(answerFrame, predsFrame,ab)
    predictions <- countPredicted(answerFrame,predsFrame,FALSE,ab)
    comparison <- createComparison(correct, predictions)
    writeOutputFrame(comparison, k = index, type = paste("comparison",ab,sep = "_"))
    writeOutputFrame(predictions, k = index, type = paste("prediction",ab,sep = "_"))
    writeOutputFrame(correct, k = index, type = paste("correct",ab,sep = "_"))
    
    vme <- countIncorrect(answerFrame, predsFrame, TRUE,ab)
    writeOutputFrame(vme, k = index, type = paste("vme",ab,sep = "_"))
    me <- countIncorrect(answerFrame, predsFrame, FALSE,ab)
    writeOutputFrame(me, k = index, type = paste("me",ab,sep = "_"))
  }
}






percentageFromFrame <- function(k, signi = NULL, type)
{
  frame = readOutputFrame(k, signi, type, folder = TEMP_FOLDER)
  sum(frame) / maxPredictions(frame, k)
}

remainingPercentage<- function(k, signi = NULL)
{
  frame = readOutputFrame(k, signi, "predictions_conformal_technical", folder = TEMP_FOLDER)
  countRemaining <- maxPredictions(frame, k) - sum(frame)
  
  countRemaining / maxPredictions(frame, k)
}




EXTRACT_SUMMARY <- function()
{
  allStats <- data.frame(
    antibiotics = integer(),
    correct = numeric(),
    ambiguous = numeric(),
    vme = numeric(),
    me = numeric(),
    nopredictions = numeric(),
    significanslevel = numeric()
  )
  for (index in range) {
    row <- c (
      as.integer(index),
      percentageFromFrame(k = index, type = "correct"),
      0,
      percentageFromFrame(k = index, type = "vme"),
      percentageFromFrame(k = index, type = "me")
    )
    row <- c(row, 1 - sum(row[2:5]))
    row <- c(row, "none")
    allStats[nrow(allStats) + 1, ] <- row
  }
  
  for (index in range) {
    for (signi in signis) {
      row <- c (
        as.integer(index),
        percentageFromFrame(k = index, signi, type = "correct_conformal"),
        percentageFromFrame(k = index, signi, type = "ambiguous_conformal"),
        percentageFromFrame(k = index, signi, type = "vme_conformal"),
        percentageFromFrame(k = index, signi, type = "me_conformal")
      )

      row <- c(row, remainingPercentage(k = index, signi))
      row <- c(row, as.numeric(fixSigni(signi)))
      allStats[nrow(allStats) + 1, ] <- row
    }
  }
  writeOutputFrame(allStats, k = NUMBER_OF_ANTIBIOTICS, type = "allstats")
  allStats
}




EXTRACT_SUMMARY_SPECIFIC <- function()
{
  # columns range
  predStatsClinical <- NULL
  predStatsTechnical <- NULL
  vmeStats <- NULL
  meStats <- NULL
  for (index in range) {
    predictions <-
      readOutputFrame(k = index,
                      type = "predictions",
                      folder = TEMP_FOLDER)
    correct <-
      readOutputFrame(k = index,
                      type = "correct",
                      folder = TEMP_FOLDER)
    vme <-
      readOutputFrame(k = index,
                      type = "vme",
                      folder = TEMP_FOLDER)
    me <- readOutputFrame(k = index,
                          type = "me",
                          folder = TEMP_FOLDER)
    max_predictions <- maxPredictions(predictions, index)
    
    vmeRow <- sum(vme) / max_predictions
    meRow <- sum(me) / max_predictions
    correctRowTechnical <- sum(correct) / sum(predictions)
    correctRowClinical <- sum(correct) / sum(predictions)
    
    
    for (signi in signis) {
      correct_type <-  "correct_conformal"
      prediction_type_technical <- "predictions_conformal_technical"
      prediction_type_clinical <- "predictions_conformal_clinical"

      predictionsConfTechnical <-
        readOutputFrame(k = index,
                        signi,
                        type = prediction_type_technical,
                        folder = TEMP_FOLDER)
      predictionsConfClinical <-
        readOutputFrame(k = index,
                        signi,
                        type = prediction_type_clinical,
                        folder = TEMP_FOLDER)
      correctConf <-
        readOutputFrame(k = index,
                        signi,
                        type = correct_type,
                        folder = TEMP_FOLDER)
      ambiguousConf <-
        readOutputFrame(k = index,
                        signi,
                        type = "ambiguous_conformal",
                        folder = TEMP_FOLDER)

      vmeConf <-
        readOutputFrame(k = index,
                        signi,
                        type = "vme_conformal",
                        folder = TEMP_FOLDER)
      meConf <-
        readOutputFrame(k = index,
                        signi,
                        type = "me_conformal",
                        folder = TEMP_FOLDER)
      correctRowTechnical <-
        c(correctRowTechnical, (sum(correctConf)+sum(ambiguousConf)) / sum(predictionsConfTechnical))
      correctRowClinical <-
        c(correctRowClinical, sum(correctConf)  / sum(predictionsConfClinical))
      vmeRow <- c(vmeRow, sum(vmeConf) / max_predictions)
      meRow <- c(meRow, sum(meConf) / max_predictions)
    }
    
    ###################
    if (is.null(predStatsClinical)) {
      predStatsClinical <- data.frame(correctRowClinical)
      rownames(predStatsClinical)[1] <- "preds"
      rownames(predStatsClinical)[2:(length(signis) + 1)] = fixSigni(signis)
    }
    else {
      predStatsClinical <- cbind(predStatsClinical, correctRowClinical)
    }
    colnames(predStatsClinical)[ncol(predStatsClinical)] <- index
    ############
    if (is.null(predStatsTechnical)) {
      predStatsTechnical <- data.frame(correctRowTechnical)
      rownames(predStatsTechnical)[1] <- "preds"
      rownames(predStatsTechnical)[2:(length(signis) + 1)] = fixSigni(signis)
    }
    else {
      predStatsTechnical <- cbind(predStatsTechnical, correctRowTechnical)
    }
    colnames(predStatsTechnical)[ncol(predStatsTechnical)] <- index
    ###############
        
    if (is.null(vmeStats)) {
      vmeStats <- data.frame(vmeRow)
      rownames(vmeStats)[1] <- "preds"
      rownames(vmeStats)[2:(length(signis) + 1)] = fixSigni(signis)
    }
    else {
      vmeStats <- cbind(vmeStats, vmeRow)
    }
    colnames(vmeStats)[ncol(vmeStats)] <- index
    ####################
    if (is.null(meStats)) {
      meStats <- data.frame(meRow)
      rownames(meStats)[1] <- "preds"
      rownames(meStats)[2:(length(signis) + 1)] = fixSigni(signis)
    }
    else {
      meStats <- cbind(meStats, meRow)
    }
    colnames(meStats)[ncol(meStats)] <- index
    
  }
  
  writeOutputFrame(predStatsTechnical, k = index, type = "predsstat_technical")
  writeOutputFrame(predStatsClinical, k = index, type = "predsstat_clinical")
  writeOutputFrame(vmeStats, k = index, type = "vmestat")
  writeOutputFrame(meStats, k = index, type = "mestat")
}



# uncomment to debug functions
#signis = c("0025")
#range <- 1:2
#range <- 10:13
#signi = "0025"
#index = 2
#index <- 3
#inputDataFrame <- inFrame
#predictionDataFrame  <- predictionFrame
#answerDataFrame <- answerFrame
#correctCount <- correct
#predictionsCount <- predictions
#bySample = FALSE



ALL <- function()
{
  COUNT()
  abs = c("AMP",	"AMC"	,"PIP"	,"TZP",	"CAZ",	"CRO",	"CTX"	,"FEP"	,"CIP"	,"OFX"	,"LVX"	,"MFX"	,"GEN"	,"TOB")
  for(ab in abs){
    print(ab)
    COUNT_WITH_AB(ab)
  }
  STATISTICS()
  EXTRACT_SUMMARY()
  EXTRACT_SUMMARY_SPECIFIC()
}
