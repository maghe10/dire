source(file='common/imports.R')

OUTPUT_DIR <- paste(modelDirectory,"input",sep="/")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
INPUT_DIR <- OUTPUT_DIR


writeCommonInputTablesAndersStyle <- function()
{
  demographicsModel <- readDemographicsTable()
  
  mode <- MODE_A
  for(mode in MODES){
    sirDataFrame <- readSirTable(mode)
    #sample columns as rownames
    rownames(sirDataFrame) <- sirDataFrame$sample
    sirDataFrame <- sirDataFrame[-1]
    antibioticsNames <- colnames(sirDataFrame)
    inputAntibioticsNames <- antibioticsNames
    numberOfInputAntibiotics <- length(inputAntibioticsNames)
    demographicsDataframe <- demographicsModel
    sirDataFrameWords <- makeWordsDataFrame(sirDataFrame,demographicsModel)

    csvTable <- sampleAsColumns(sirDataFrameWords)
    file = paste(OUTPUT_DIR,paste("sirAntibioticsModelWords_",sub(" ", "-", mode),".csv",sep=""),  sep="/")
    write.csv2(x=csvTable,file=file,row.names = FALSE)
    
    demographicsModelNoDate <- demographicsModel
    demographicsModelNoDate$date <- "<unk>"
    sirDataFrameWords <- makeWordsDataFrame(sirDataFrame,demographicsModelNoDate)
    csvTable <- sampleAsColumns(sirDataFrameWords)
    file = paste(OUTPUT_DIR,paste("sirAntibioticsModelWordsNoDate_",sub(" ", "-", mode),".csv",sep=""),  sep="/")
    write.csv2(x=csvTable,file=file,row.names = FALSE)
  }
}

makeWordsDataFrame <- function(sirDataFrame,demographicsDataframe)
{
  tmp <- sirDataFrame
  for(column in colnames(sirDataFrame)){
    tmp[,column] <- paste(column,sirDataFrame[,column],sep="_")
  }
  
  antibiotics <- c()
  row <- "2"
  for(row in rownames(tmp)){
    antibiotic <- paste(tmp[row,],collapse = ',')    
    antibiotics <- c(antibiotics,antibiotic)
  }
  result <- data.frame(antibiotics)
  n <- rep(ncol(tmp))
  result <- cbind(result,n)  
  # x: ['SE 72 F 2018-10', 'SE 63 M 2016-03']}
  x <- paste("SE",demographicsDataframe$age,demographicsDataframe$sex,substr(demographicsDataframe$date,start=0,stop=7))
  result <- cbind(result,x)  
  
  colnames(result) = c("Antibiotic","n","x")
  result
}


abListToabHeading <- function(abList)
{
  result <- abList[[1]]
  if(length(abList)>1){
    for(index in 2:length(abList)){
      result <- paste(result,abList[[index]],sep="_")
    }
  }
  result
}

asColumnNameList <- function (comb){
  aList <- c()
  for(index in 1:ncol(comb)){
    aList <- c(aList,abListToabHeading(as.list(comb[,index])))
  }
  aList
}

createDataframeWithColumns <- function (k,inputAntibioticsNames,sirDataFrame)
{
  comb <- combn(inputAntibioticsNames,m=k)
  columnnames <- asColumnNameList(comb) 
  df <- data.frame()
  for(i in 1:nrow(sirDataFrame)){
    row <- rep(NA,length(columnnames)) 
#    for (j in 1:length(columnnames)) {
#      word <- sirToWord(sirDataFrame[i,],columnnames[[j]])
#      row[[j]] <- word
#    }
    df <- rbind(df,row)
  }
  colnames(df) <- columnnames
  rownames(df) <- rownames(sirDataFrame)
  df
}


createDataframeWithRows <- function (k,inputAntibioticsNames)
{
  comb <- t(combn(inputAntibioticsNames,m=k))
  comb
}


#
#
#

readAndersFrame <- function(fileNameIn= "sirAntibioticsModelWords_Mode-A.csv")
{
  read.csv2(file=paste(INPUT_DIR,fileNameIn,sep = "/"))
  
}

# Convert anders frame to juan frame
#
# patient information from the format in Anders implementation
# to the format in Juan implementation
#
# Example:
# Anders: SE 70 F 2014-10
# Juan: SE F 70 2014_10
# 
# Anders:
# Antibiotic: comma separated list of all antibiotic
# x: patient information in anders format
#
#
convertFrame <- function(frame,inputAb)
{
  antibioticsSplit <- frame  %>%
    mutate(
      tokens = str_split(Antibiotic, ","),
      in_inputAb  = map(tokens, ~ .x[str_extract(.x, "^[A-Z]+") %in% inputAb]),
      not_inputAb = map(tokens, ~ .x[!(str_extract(.x, "^[A-Z]+") %in% inputAb)]),
      
      in_inputAb  = map_chr(in_inputAb,  ~ paste(.x, collapse = " ")),
      not_inputAb = map_chr(not_inputAb, ~ paste(.x, collapse = " "))
    ) %>% select(in_inputAb,not_inputAb)
  
  convertedPatientID <- frame %>%
    separate(
      x,
      into = c("country", "age", "sex", "month"),
      sep = " ",
      remove = FALSE
    ) %>%
    mutate(
      age  = as.integer(age),
      month_mod = str_replace(month, "-", "_")  # "2022-09" → "2022_09"
    ) %>%
    mutate(
      month_mod = str_replace(month_mod, "_0", "_")  # "2022_09" → "2022_9"
    )  %>%
    #Remove when model is updated to include our dates. 
    mutate(
      month_mod = "2020_12"  # Last date juan handles
    ) %>%
    mutate(
      combined = paste(country, sex, age, month_mod, sep = " ")
    ) %>% select(combined)
  
  
  anotherFrame <- data.frame(rep("ESCCOL",nrow(frame)),convertedPatientID,antibioticsSplit)
  colnames(anotherFrame) <- c("bacteria","patientdata","in_ab","out_ab")
  anotherFrame
}

# Write a frame in the juan format, i.e. without row names and header.
# The columns, separated by comma are
#
# <species>,<patientid>, <input_ab, <output_ab>
#
# Example: 
# ESCCOL, HR F 70 2014_10, CRO_S AMX_R CIP_S GEN_S, CAZ_S AMP_S
#
writeJuanFrame <- function(frame,outputFile = fileNameOut,outputFolder = INPUT_DIR )
{
  if (!dir.exists(outputFolder)) {
    dir.create(outputFolder, recursive = TRUE)
  }
  
  
  write.table(frame,file=paste(outputFolder,outputFile,sep = "/"),
              sep = ",",
              col.names = FALSE,
              row.names = FALSE,
              quote = FALSE)
}



convertAndersJuanFrame <- function(name,
                                   mode,
                                   indices  = 4:13,
                                   parallel = FALSE,
                                   workers  = NULL) {
  # name  <- "sirAntibioticsModelWords"
  # mode  <- "Mode-A"
  # indices <- 12:13
  # parallel <- TRUE/FALSE
  
  # ------------------------------------------------------------------
  # 1) Read the Anders frame once
  # ------------------------------------------------------------------
  input_file <- sprintf("%s_%s.csv", name, mode)
  frame <- readAndersFrame(input_file)
  
  if (!"sample" %in% names(frame)) {
    stop("convertAndersJuanFrame: 'frame' must have a 'sample' column.")
  }
  
  # Split once per sample so we don't repeatedly subset by row index
  rows_by_sample <- split(frame, frame$sample)  # named list: names = sample IDs
  
  # Output directory for this mode
  out_dir_mode <- file.path(INPUT_DIR, mode)
  if (!dir.exists(out_dir_mode)) {
    dir.create(out_dir_mode, recursive = TRUE)
  }
  
  # ------------------------------------------------------------------
  # 2) Optionally set up parallel plan
  # ------------------------------------------------------------------
  if (parallel) {
    # If workers not specified, use all available minus one
    if (is.null(workers)) {
      workers <- max(1L, future::availableCores() - 1L)
    }
    future::plan(future::multisession, workers = workers)
  }
  
  # ------------------------------------------------------------------
  # 3) Loop over all requested indices (e.g. 12, 13)
  # ------------------------------------------------------------------
  for (i in indices) {
    # --------------------------------------------------------------
    # 3a) Read combinations for this i only once
    # --------------------------------------------------------------
    model_output_file <- file.path(INPUT_DIR, sprintf("modelOutput%d.csv", i))
    
    # Read just header to get column names (faster)
    combos_header <- read.csv2(file = model_output_file, nrows = 1)
    combinations  <- colnames(combos_header)[-1]   # drop first col (sample)
    
    # Pre-parse all combinations into parts once
    parts_list <- strsplit(combinations, "_", fixed = FALSE)
    
    message("Processing i = ", i, " with ", length(combinations), " combinations")
    message("Number of samples: ", length(rows_by_sample))
    
    # --------------------------------------------------------------
    # 3b) Define per-sample worker function
    # --------------------------------------------------------------
    process_one_sample <- function(row_df, sample_id) {
      # row_df is a 1-row data frame for this sample
      # sample_id is the name from split() (character)
      
      frames_list <- lapply(parts_list, function(parts) {
        convertFrame(row_df, parts)
      })
      
      juanFrame <- dplyr::bind_rows(frames_list)
      
      out_file <- sprintf("%sJuan_%d_%s.csv", name, i, sample_id)
      
      writeJuanFrame(
        frame        = juanFrame,
        outputFile   = out_file,
        outputFolder = out_dir_mode
      )
      
      invisible(NULL)
    }
    
    # --------------------------------------------------------------
    # 3c) Either parallel or sequential over samples
    # --------------------------------------------------------------
    if (parallel) {
      # Parallel over samples with furrr
      furrr::future_walk2(
        .x        = rows_by_sample,
        .y        = names(rows_by_sample),
        .f        = process_one_sample,
        .progress = TRUE,
        .options  = furrr::furrr_options(seed = TRUE)
      )
    } else {
      # Sequential over samples with nice logging using purrr::iwalk
      purrr::iwalk(
        .x = rows_by_sample,
        .f = function(row_df, sample_id) {
          message("  Sample ", sample_id, " for i = ", i)
          process_one_sample(row_df, sample_id)
        }
      )
    }
  }
  
  invisible(NULL)
}



ALL  <- function()
{
  writeCommonInputTablesAndersStyle()
  for(mode in MODES) {
    convertAndersJuanFrame("sirAntibioticsModelWords",mode,parallel = TRUE)
  }
}
