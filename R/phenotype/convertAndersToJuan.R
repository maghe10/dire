library(stats)
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

source(file='common.R')

INPUT_DIR <- paste(modelDirectory,"input",sep="/")
fileNameIn <- "sirAntibioticsModelWords_Mode-A.csv"
fileNameOut <- "input_Mode-A.csv"

readAndersFrame <- function(fileNameIn= "sirAntibioticsModelWords_Mode-A.csv")
{
  read.csv2(file=paste(INPUT_DIR,fileNameIn,sep = "/"))
  
}





# Convert patient information from the format in Anders implementation
# to the format in Juan implementation
#
# Example:
# Anders: SE 70 F 2014-10
# Juan: SE F 70 2014_10
#
# Operates on the x field
convertPatientID <-function(frame)
{
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
    ) %>%
  #Remove when model is updated to include our dates. 
  mutate(
    month_mod = "2020_12"  # Last date juan handles
   ) %>%
   mutate(
      combined = paste(country, sex, age, month_mod, sep = " ")
    ) %>% select(combined)
  frame$x <- convertedPatientID$combined
  frame
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

#
# Convert the input file in Anders format to a similar, but with Juan patient information to be used in my python script that calls juans implementation.
# 
#
convertinputFrame <-function(name,mode)
{
  frame <- readAndersFrame(paste(name,"_",mode,".csv",sep=""))
  frame <- convertPatientID(frame)
  write.csv2(frame,paste(INPUT_DIR,paste(name,"Juan2020_12","_",mode,".csv",sep=""),sep = "/"),row.names = FALSE)
}

# Make an example that can be run with
# directly with juans inplementaton
#
CONVERTINPUTFRAMES <- function()
{
  convertinputFrame("sirAntibioticsModelWords","Mode-A")
  convertinputFrame("sirAntibioticsModelWords","Mode-B")
  convertinputFrame("sirAntibioticsModelWords","Mode-C")
}



# Write a frame in the juan format, i.e. without row names and header.
# The columns, separated by comma are
#
# <species>,<patientid>, <input_ab, <output_ab>
#
# Example: 
# ESCCOL, HR F 70 2014_10, CRO_S AMX_R CIP_S GEN_S, CAZ_S AMP_S
#
WriteJuanFrame <- function(frame)
{
  write.table(frame,file=paste(INPUT_DIR,fileNameOut,sep = "/"),
              sep = ",",
              col.names = FALSE,
              row.names = FALSE,
              quote = FALSE)
}

# Make an example in juan format based on my samples to be called
# directly with juans inplementaton
#
# It outputs a random selection of all combinations with 4 input antibiotics
CONVERT_EXAMPLE <- function()
{

  frame <- readAndersFrame()
  
  combinations <- colnames(read.csv2(file=paste(INPUT_DIR,"modelOutput4.csv",sep = "/")))[-1]
  combinations %>% stringr::str_split("_", simplify = TRUE) %>% as.vector()


  frame2 <- tibble(text = combinations) %>%
    mutate(parts = strsplit(text, "_")) 

  frames_list <- lapply(seq_len(nrow(frame2)), function(i) {
    convertFrame(frame, frame2$parts[[i]])
  })
  
  returnFrame <- bind_rows(frames_list)
  returnFrame <- returnFrame %>% slice_sample(n = 1000)


  WriteJuanFrame(returnFrame)
}