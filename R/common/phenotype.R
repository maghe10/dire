source("common/imports.R")
library(tibble)
library(tidyr)


OUTPUT_DIR <- processedRootRcommon
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

ATU <- "atu"
BREAKPOINTS <- "breakpoints"
DEMOGRAPICS <- "demographicstable"
ESBL <- "esbl"
MILLIMETERS <- "millimetertable"
MODEL_MILLIMETERS <- "modelmillimetertable"
RESCALED_MILLIMETERS <- "millimetertableRescaled"
RESCALED_MODEL_MILLIMETERS <- "modelmillimetertableRescaled"
SIR <- "sirAntibiotics"
SIR_model <- "sirAntibioticsModel"
sir_mode <- function(mode){
  mode <- gsub(" ", "-", mode)
  paste(SIR_model,mode,sep="_") }
SIR_A = sir_mode(MODE_A)
SIR_B= sir_mode(MODE_B)
SIR_C= sir_mode(MODE_C)

ALL_NAMES <- c(
  ATU ,
  BREAKPOINTS ,
  DEMOGRAPICS,
  ESBL,
  MILLIMETERS,
  MODEL_MILLIMETERS,
  RESCALED_MILLIMETERS,
  RESCALED_MODEL_MILLIMETERS,
  SIR,
  SIR_model,
  SIR_A ,
  SIR_B,
  SIR_C
)

rescaleMillimeters <- function(millimeterTableToRescale)
{
  colMax <- function (colData) {
    apply(colData, MARGIN=c(2), max)
  }
  colMin <- function (colData) {
    apply(colData, MARGIN=c(2), min)
  }
  
  abMax <- colMax(millimeterTableToRescale)
  abMin <- colMin(millimeterTableToRescale)
  
  returnValue <- millimeterTableToRescale
  
  for(col in 1:ncol(returnValue)){
    returnValue[,col] <-  (millimeterTableToRescale[,col]-abMin[col])/(abMax[col]-abMin[col])
  }
  returnValue
}


#from generate_model_input.R
importFromExcel <- function()
{
  getIncludedSubset <-function(measurementTable)
  {
    exkludera <- measurementTable$Exkludera=="Ja"
    exkludera[is.na(exkludera)] <- FALSE
    includedMeasurementTable <- measurementTable[!exkludera,]
    rownames(includedMeasurementTable) <- includedMeasurementTable$Studienummer
    includedMeasurementTable
  }
  
  getMillimetertable <- function(includedMeasurementTable)
  {
    millimeterTable <- includedMeasurementTable[,allAntibioticsWithZonesInputColumns]
    rownames(millimeterTable) <- includedMeasurementTable$Studienummer
    millimeterTable
  }
  
  getDemographicsTable <- function(includedMeasurementTable)
  {
    demographicsModel <- includedMeasurementTable[,demographicsInputColumns]
    demographicsModel$Datum <- as.Date(as.integer(demographicsModel$Datum),origin=as.Date("1900-01-01"))
    demographicsModel$Kön[which(demographicsModel$Kön=="K")] <- "F"
    colnames(demographicsModel) <- c("date","sex","age")
    
    demographicsModel
  }
  
  write_csv2_sample <- function(x, name) {
    x_out <- data.frame(
      sample = rownames(x),
      x,
      row.names = NULL,
      check.names = FALSE
    )
    
    write.csv2(x_out, file=paste(OUTPUT_DIR,sprintf("%s.csv",name),sep = "/"), row.names = FALSE)
  }
  
  


  writeSirAntibiotics <- function(modelMillimeterTable,modelLimitTable)
  {
    modes <- colnames(modelLimitTable[2:ncol(modelLimitTable)])
    zones <- modelMillimeterTable
    for(mode in modes){
      limit <- modelLimitTable[,mode]
      sir <-data.frame(zones)
      for(abindex in  1:ncol(zones)){
        sir[zones[,abindex]<limit[abindex],abindex] <- "R"
        sir[zones[,abindex]>=limit[abindex],abindex] <- "S"
      }
      write_csv2_sample(sir,name = sir_mode(mode))
    }
  }
 
  phenotypeExcel =  paste(processedRootExcel, "ModelInput.xlsx", sep="/")
  
  SHEET_BREAKPOINTS_INVASIV <- "Breakpoints invasiv"
  SHEET_BREAKPOINTS_UTI <- "Breakpoints UTI"
  measurementTable <- sheetToDataFrame(read_xlsx(phenotypeExcel,sheet = SHEET_MEASUREMENTS_RAW))
  modelLimitTable <- sheetToDataFrame(read_xlsx(phenotypeExcel,sheet = SHEET_MODEL_LIMIT))
  headingToAntibioticTable <- sheetToDataFrame(read_xlsx(phenotypeExcel,sheet = SHEET_HEADING_TO_ANTIBIOTIC))

  breakpointsTableInvasive <- sheetToDataFrame(read_xlsx(phenotypeExcel,sheet = SHEET_BREAKPOINTS_INVASIV))
  breakpointsTableUti <- sheetToDataFrame(read_xlsx(phenotypeExcel,sheet = SHEET_BREAKPOINTS_UTI))
  
  #Missing from invasive
  diff <- setdiff(
    unique(breakpointsTableUti$Antibiotika),
    unique(breakpointsTableInvasive$Antibiotika)
  )
  
  breakpointsTable <- rbind(breakpointsTableInvasive,breakpointsTableUti[breakpointsTableUti$Antibiotika %in% diff,])
  
  breakpointsTable <- breakpointsTable %>% 
    filter_and_drop(Version,EUCAST_VERSION) %>% 
    rename(antibiotic = Antibiotika) %>%
    dplyr::mutate(
      antibiotic = dplyr::recode(
        antibiotic,
        "CN"  = "GEN",
        "LEV" = "LVX",
        "PRL" = "PIP"
      )
    ) %>% 
    filter(antibiotic %in% allAntibioticsWithBreakpointsNames)
  
  
  
  breakpointsTable %>% write.csv2(
    file.path(OUTPUT_DIR, sprintf("%s.csv",BREAKPOINTS)),
    na = "",
    row.names = FALSE
  )

  
  modelLimitTable <- modelLimitTable %>%
    mutate(Antibiotic = factor(Antibiotic, levels = allAntibioticsInModelModelNames)) %>%
    arrange(Antibiotic)

  includedMeasurementTable <- getIncludedSubset(measurementTable)
  millimeterTable <- getMillimetertable(includedMeasurementTable)
  
  # rename
  lookup <- setNames(
    headingToAntibioticTable$Antibiotic,
    headingToAntibioticTable$Heading
  )
  new_names <- lookup[names(millimeterTable)]
  new_names[is.na(new_names)] <- names(millimeterTable)[is.na(new_names)]
  names(millimeterTable) <- new_names
  
  #model subset
  modelMillimeterTable <- millimeterTable %>%
    select(all_of(allAntibioticsInModelModelNames))
  write_csv2_sample(millimeterTable,MILLIMETERS)
  write_csv2_sample(rescaleMillimeters(millimeterTable),RESCALED_MILLIMETERS)
  write_csv2_sample(modelMillimeterTable,MODEL_MILLIMETERS)
  write_csv2_sample(rescaleMillimeters(modelMillimeterTable),RESCALED_MODEL_MILLIMETERS)
  writeSirAntibiotics(modelMillimeterTable,modelLimitTable)
  
  get_SIR <- function(mm, S, R) {
    ifelse(
      is.na(mm), NA,
      ifelse(mm < R, "R",
             ifelse(mm >= S, "S", "I"))
    )
  }
  
  sirAntibiotics <- c(c("F", "MEL", "CFR", "W", "SXT","MEM"),colnames(modelMillimeterTable))
  
  sirWide <- millimeterTable[,unique(breakpointsTable$antibiotic)] %>%
    rownames_to_column("sample") %>%
    pivot_longer(-sample, names_to = "antibiotic", values_to = "mm") %>%
    left_join(breakpointsTable, by = c("antibiotic")) %>%
    mutate(
      SIR = get_SIR(mm,as.numeric(S),as.numeric(R))
    ) %>% select(sample, antibiotic, SIR) %>%
    pivot_wider(names_from = antibiotic, values_from = SIR)
  
  writeSir <- function(sirWide)
  {
    fileName <- sprintf("%s.csv",SIR)
    write.csv2(sirWide,file=paste(OUTPUT_DIR,fileName,sep = "/"),row.names = FALSE)
  }
  writeSir(sirWide)
  
  
    #SIR 
  sirWideModel <- modelMillimeterTable %>%
    rownames_to_column("sample") %>%
    pivot_longer(-sample, names_to = "antibiotic", values_to = "mm") %>%
    left_join(breakpointsTable, by = c("antibiotic")) %>%
    mutate(
      SIR = get_SIR(mm,as.numeric(S),as.numeric(R))
    ) %>% select(sample, antibiotic, SIR) %>%
    pivot_wider(names_from = antibiotic, values_from = SIR)
  
  writeSirModel <- function(sirWideModel)
  {
    fileName <- sprintf("%s.csv",SIR_model)
    write.csv2(sirWideModel,file=paste(OUTPUT_DIR,fileName,sep = "/"),row.names = FALSE)
  }
  writeSirModel(sirWideModel)
  
  # ATU 
  atuWide <- modelMillimeterTable %>%
    rownames_to_column("sample") %>%
    pivot_longer(-sample, names_to = "antibiotic", values_to = "mm") %>%
    left_join(breakpointsTable, by = c("antibiotic")) %>%
    mutate(in_atu = !is.na(ATU_LOW) & mm >= ATU_LOW & mm <= ATU_HIGH) %>%
    select(sample, antibiotic, in_atu) %>%
    pivot_wider(names_from = antibiotic, values_from = in_atu)
  writeAtu <- function(atuWide)
  {
    fileName <- sprintf("%s.csv",ATU)
    write.csv2(atuWide,file=paste(OUTPUT_DIR,fileName,sep = "/"),row.names = FALSE)
  }
  writeAtu(atuWide)
  
  demographicsTable <- getDemographicsTable(includedMeasurementTable)
  write_csv2_sample(demographicsTable,DEMOGRAPICS)

  # ESBL
  readESBL <- function(includedMeasurementTable)
  {
    esblSubset <- includedMeasurementTable %>% select(
      c("Studienummer","DDT(p/n)","AMPC(p/n)","KISS II-MEM")
    )    
    includedMeasurementTable %>%  
      mutate(
        ESBL_A = `DDT(p/n)` == "p",
        ESBL_M = `AMPC(p/n)` == "p",
        ESBL_CARBA = `KISS II-MEM` < 25)  %>%
      select(ESBL_A,ESBL_M,ESBL_CARBA)
  }
  
  esbl <- readESBL(includedMeasurementTable)
  write_csv2_sample(esbl,ESBL)

}


# DUMP_MODEL_MM_FRAME <-function()
# {
#   frame <- fetchTables()[["MM"]]
#   sample <- row.names(frame)
#   frame <- cbind(sample,frame)
#   # colnames(frame)
#   write.csv2(frame,file.path(getCommonModelFolder(),"modelzonemillimeters.csv"),row.names = FALSE)
# }

ALL <- function()
{
  importFromExcel()
  HEAD()
}

HEAD <- function()
{
  invisible(lapply(ALL_NAMES, function(name) {
    file <- file.path(OUTPUT_DIR, paste0(name, ".csv"))
    cat("\n===", paste0(name, ".csv"), "===\n")
    print(head(read.csv2(file)))
  }))
}
