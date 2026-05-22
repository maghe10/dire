source('model/modelcommon.R')
source('model/make_models_output_file.R')

library(data.table)
library(dplyr)
library(stringr)



add_strict_confpred <- function(df, conf_cols = NULL) {
  # If user does not supply a list of confpred columns, detect them
  if (is.null(conf_cols)) {
    conf_cols <- grep("^confpred_", names(df), value = TRUE)
  }
  
  df %>%
    mutate(across(
      all_of(conf_cols),
      ~ {
        conf <- .x
        p    <- pred
        
        single_conf <- nchar(conf) == 1 & !is.na(conf)
        single_pred <- nchar(p) == 1 & !is.na(p)
        
        ifelse(single_conf & single_pred & conf == p, conf, "")
      },
      .names = "{.col}_strict"
    ))
}


RUN <- function()
{
  modelFile <- paste("sirAntibioticsModel","_",sub(" ", "-", MODE_A),".csv",sep="")
  sirAntibioticsModel <- read.csv2(row.names=1,paste(modelDirectory,"input",modelFile,sep="/"))
  realSampleNumbers <- as.integer(rownames(sirAntibioticsModel))
  modes <- MODES
  #modes <- c(MODE_C)
  for(mode in modes)
    #mode <- "Mode-A"
  for(k in 4:13){
    #k <- 13
    template_file <- paste(modelDirectory,"input",sprintf("modelOutput%d.csv",k),sep="/")
    stopifnot(file.exists(template_file))
    fileName <- paste("sirAntibioticsModelWordsJuan_",k,"_1_cp.csv",sep="")
    file <- paste(getConformalPredictionFolder(mode),fileName,sep="/")
    stopifnot(file.exists(file))
    input_dir <- getConformalPredictionFolder(mode)
    output_dir <- getPredictionAltFolder(mode)
    print(input_dir)
    print(output_dir)
    print(k)
    print(template_file)
    make_model_output_files(input_dir = input_dir,output_dir = output_dir,k=k,sample_range = 1:99,true_sample_ids = realSampleNumbers,template_file = template_file)
  }  
}




merge_raw <- function()
{
  realSampleNumbers <- realSampleNumbers()
  allmodes <- lapply(MODES,function(mode){
    allnoinputab <- lapply (4:13,function(k){
      allSamplesFrames <- lapply(1:99,function (p){
        read.csv(file.path(getConformalPredictionFolder(mode), sprintf("sirAntibioticsModelWordsJuan_%d_%d_cp.csv",k,p)))[-1,-c(1,2)]
      })
      names(allSamplesFrames) <- realSampleNumbers
      bind_rows(allSamplesFrames, .id = "sample")
    })
    bind_rows(allnoinputab,.id = "noinputab")
  })
  names(allmodes) <- MODES
  bind_rows(allmodes,.id = "mode")
}


realSampleNumbers <- function()
{
  modelFile <- paste("sirAntibioticsModel","_",sub(" ", "-", MODE_A),".csv",sep="")
  sirAntibioticsModel <- read.csv2(row.names=1,paste(modelDirectory,"input",modelFile,sep="/"))
  as.integer(rownames(sirAntibioticsModel))
}

onn2Softmax_vec_scan <- function(x) {
  # Collapse all strings, strip brackets, parse all numbers in one scan
  txt <- gsub("\\[|\\]", "", paste(x, collapse = "\n"))
  nums <- scan(text = txt, sep = ",", quiet = TRUE)  # numeric vector length 2N
  
  mat <- matrix(nums, ncol = 2, byrow = TRUE)
  a <- mat[, 1]; b <- mat[, 2]
  p1 <- 1 / (1 + exp(b - a))
  cbind(p1, 1 - p1)
}  

applySoftmax_cached <- function(df) {
  x <- df$Output_neural_networks
  ux <- unique(x)
  Pu <- onn2Softmax_vec_scan(ux)          # compute once per unique string
  map <- match(x, ux)
  df$psSoftmax <- Pu[map, 1]
  df$prSoftMax <- Pu[map, 2]
  df
}

applySoftmax <-function(df)
{
  frame <- df
  #frame <- head(df)
  P <- onn2Softmax_vec(frame$Output_neural_networks)
  frame$psSoftmax <- P[,1]
  frame$prSoftmax <- P[,2]
  frame
}

MERGE_RAW <- function()
{
  df <- merge_raw()
  saveRDS(df,file.path(getCommonModelFolder(),"merged_cp_raw.rds"))
  #df <-readRDS(file.path(getCommonModelFolder(),"merged_cp_raw.rds"))
  df_with_softMax <- applySoftmax_cached(df)
  saveRDS(df_with_softMax,file.path(getCommonModelFolder(),"merged_cp_raw_with_softmax.rds"))
}


MERGE <- function()
{
  modelFile <- paste("sirAntibioticsModel","_",sub(" ", "-", MODE_A),".csv",sep="")
  sirAntibioticsModel <- read.csv2(row.names=1,paste(modelDirectory,"input",modelFile,sep="/"))
  realSampleNumbers <- as.integer(rownames(sirAntibioticsModel))
  
  modes <- MODES
  
  for (mode in modes) {
    #mode <- MODE_A
    input_dir  <- getConformalPredictionFolder(mode)
    output_dir <- getPredictionAltFolder(mode)
    
    ks <- 4:13
    
    big_dt <- rbindlist(
      lapply(ks, function(k) {
        message("mode = ", mode, ", k = ", k)
        
        dt <- fread(file.path(output_dir,
                              sprintf("modelOutput_long_%d.csv", k)))
        dt[, noinputab := k]
        dt
      }),
      use.names = TRUE,
      fill = TRUE
    )
    
    big_dt <- add_strict_confpred(big_dt)
    outRDS <- file.path(output_dir,
                        sprintf("modelOutput_long_merged.rds"))
    # R-object
    saveRDS(big_dt,
            outRDS)

    big_dt_read <- readRDS(outRDS)
    #verify
    stopifnot(all(big_dt_read==big_dt))
    #colnames(big_dt_read)
    #head(big_dt)
  }
}


ALL <- function()
{
  RUN()
  MERGE()
  MERGE_RAW()
}
