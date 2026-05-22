library(readxl)
library(dplyr)

source(file='common.R')

dir = processedRootExcel

allAntibioticsWithZonesInputColumns <- c("KISS I-F",	"KISS I-MEL","KISS I-CFR",	"KISS I-W",	"KISS I-CIP",	"KISS I-AMC",	"KISS II-CTX",	"KISS II-CAZ"	,"KISS II-MEM",	"KISS II-TOB",	"KISS II-TZP",	"KISS II-SXT",	
                                         "DDT-FOX","DDT-FEP", "Studie-1-AMP"	,"Studie-1-PRL",	"Studie-1-CN"	,"Studie-1-CRO",	"Studie-2-LEV","Studie-2-MFX",	"Studie-2-OFX",	"Studie-2-NA")
allAntibioticsInModelModelNames <- c("AMP","AMC","PIP","TZP","CAZ","CRO","CTX","FEP","CIP","OFX","LVX","MFX","GEN","TOB")
allAntibioticsWithBreakpointsNames <- c("AMP","AMC","PIP","TZP","CAZ","CRO","CTX","FEP","CIP","OFX","LVX","MFX","GEN","TOB","MEM","SXT","W","F","MEL","CFR")

allMeasurableAntibioticsNames <- c("AMP","AMC","PIP","TZP","MEL","CFR","FOX","CAZ","CAZ30","CRO","CTX","CTX30","FEP","MEM","CIP","NAL","OFX","LVX","MFX","GEN","TOB","F","W","SXT")
length(allMeasurableAntibioticsNames)

demographicsInputColumns <-c("Datum", "Kön",	"Ålder")

phenotypeExcel <-  paste(processedRootExcel, "ModelInput.xlsx", sep="/")
SHEET_MEASUREMENTS_RAW <- "MeasurementsRaw"
SHEET_MODEL_LIMIT <- "ModelLimit"
SHEET_HEADING_TO_ANTIBIOTIC <- "HeadingToAntibiotic"

QCLimitsExcel <-  paste(processedRootExcel, "QC-limits.xlsx", sep="/")
SHEET_RANGE_HIGH <- "Range_high"
SHEET_RANGE_LOW <- "Range_low"
SHEET_TARGET_HIGH <- "Target_high"
SHEET_TARGET_LOW <- "Target_low"


sheetToDataFrame <- function(sheet)
{
  table <- data.frame(sheet)
  colnames(table) <- colnames(sheet)
  table
}

changeToAntibioticTypeOrder <- function(table,column,orderTo)
{
  reverse_order <- c()
  for(i in 1:length(orderTo)){
    match <- which(column == orderTo[i])
    #    print(match)
    reverse_order <- c(reverse_order,match)
  }
  table[reverse_order,]
}
readPhenotypeTable <- function()
{
  measurementTable <- sheetToDataFrame(read_xlsx(phenotypeExcel,sheet = SHEET_MEASUREMENTS_RAW))
  measurementTable
}

readHeadingToAntibioticTable <-function ()
{
  headingToAntibioticTable <- sheetToDataFrame(read_xlsx(phenotypeExcel,sheet = SHEET_HEADING_TO_ANTIBIOTIC))
  headingToAntibioticTable <- changeToAntibioticTypeOrder(headingToAntibioticTable,headingToAntibioticTable$Antibiotic,allAntibioticsInModelModelNames)
  headingToAntibioticTable
}

readHeadingToAntibioticTableAll <-function ()
{
  headingToAntibioticTable <- sheetToDataFrame(read_xlsx(QCLimitsExcel,sheet = SHEET_HEADING_TO_ANTIBIOTIC))
  headingToAntibioticTable <- changeToAntibioticTypeOrder(headingToAntibioticTable,headingToAntibioticTable$Antibiotic,allMeasurableAntibioticsNames)
  headingToAntibioticTable
}

recodeAntibiotic <- function(table)
{
  headingTable <- readHeadingToAntibioticTableAll()
  table %>%
    mutate(
      antibiotic = recode(
        antibiotic,
        !!!setNames(headingTable$Antibiotic, headingTable$Heading),
        .default = antibiotic
      )
    )
}

QCmillimetertable <- function()
{
  xxx <- readPhenotypeTable() %>% filter(Studienummer %in% c("CCUG_17620","CCUG_30600"))
  xxx_long <- xxx %>%
    select(Studienummer, where(is.numeric)) %>%
    pivot_longer(
      cols = -Studienummer,
      names_to = "antibiotic",
      values_to = "value",
      values_drop_na = TRUE
    )

  xxx_long %>% recodeAntibiotic()
}

readQCTargetTable <- function()
{
  limits <- c(SHEET_TARGET_HIGH,SHEET_TARGET_LOW)
  frames <- lapply(limits,FUN = function(sheet){
    xxx <- sheetToDataFrame(read_xlsx(QCLimitsExcel,sheet = sheet)) 
    xxx %>%  pivot_longer(colnames(xxx)[-1],names_to="antibiotic",values_to = "limit",values_drop_na = TRUE)
  })
  names(frames) <- c("high","low")
  bind_rows(frames,.id = "kind") %>% recodeAntibiotic()
}

readQCRangeTable <- function()
{
  limits <- c(SHEET_RANGE_HIGH,SHEET_RANGE_LOW)
  frames <- lapply(limits,FUN = function(sheet){
    xxx <- sheetToDataFrame(read_xlsx(QCLimitsExcel,sheet = sheet)) 
    xxx %>%  pivot_longer(colnames(xxx)[-1],names_to="antibiotic",values_to = "limit",values_drop_na = TRUE)
  })
  names(frames) <- c("high","low")
  bind_rows(frames,.id = "kind") %>% recodeAntibiotic()
}

antibioticsOrder <-function(){
  allMeasurableAntibioticsNames
}


