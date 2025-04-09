library(stats)
library(readxl)
library(dplyr)

source(file='common.R')

dir = processedRootExcel

allAntibioticsWithZonesInputColumns <- c("KISS I-F",	"KISS I-MEL","KISS I-CFR",	"KISS I-W",	"KISS I-CIP",	"KISS I-AMC",	"KISS II-CTX",	"KISS II-CAZ"	,"KISS II-MEM",	"KISS II-TOB",	"KISS II-TZP",	"KISS II-SXT",	
                                         "DDT-FOX","DDT-FEP", "Studie-1-AMP"	,"Studie-1-PRL",	"Studie-1-CN"	,"Studie-1-CRO",	"Studie-2-LEV","Studie-2-MFX",	"Studie-2-OFX",	"Studie-2-NA")

#allAntibioticsInModel <- c("AMP",	"AMC","PRL","TZP", "CAZ","CRO",	"CTX"	,"FEP", "CIP","OFX"	,"LEV","MFX","CN", "TOB" )

allAntibioticsInModelModelNames <- c("AMP","AMC","PIP","TZP","CAZ","CRO","CTX","FEP","CIP","OFX","LVX","MFX","GEN","TOB")



#allAntibioticsWithInvasiveSirInputColumns <- c("CIP",	"AMC",	"CTX",	"CAZ"	,"MEM",	"TOB",	"TZP",	"SXT",	
#                                       "FEP", "AMP"	,"PRL",	"CN"	,"CRO",	"LEV","MFX",	"OFX")

#allAntibioticsWithUTISirInputColumns <- c("F",	"MEL","CFR",	"W")

#studieNummerColumn <- c("Studienummer")
demographicsInputColumns <-c("Datum", "Kön",	"Ålder")


#phenotypeExcel =  paste(processedRootExcel, "Tolkad_data.xlsx", sep="/")
#zoneMillisSheet <- read_xlsx(phenotypeExcel,sheet = 'MeasureIncluded')
#sirSheetInvasiv <- read_xlsx(phenotypeExcel,sheet = 'Tolka Invasiv')
#sirSheetUTI <- read_xlsx(phenotypeExcel,sheet = 'Tolka UTI')
#measurementsRawSheet <- read_xlsx(phenotypeExcel,sheet = 'MeasurementsRaw')

SHEET_MEASUREMENTS_RAW = "MeasurementsRaw"
SHEET_MODEL_LIMIT = "ModelLimit"
SHEET_HEADING_TO_ANTIBIOTIC = "HeadingToAntibiotic"

OUTPUT_DIR <- paste(modelDirectory,"input",sep="/")

sheetToDataFrame <- function(sheet)
{
  table <- data.frame(sheet)
  colnames(table) <- colnames(sheet)
  table
}

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
  demographicsModel  <- cbind(includedMeasurementTable$Studienummer,demographicsModel)
  colnames(demographicsModel) <- c("sample","date","sex","age")
  
  demographicsModel
}


writeMillimeterTable <- function(includedMeasurementTable)
{
  millimeterTable <- getMillimetertable(includedMeasurementTable)

  millimeterTable <- cbind(rownames(millimeterTable),millimeterTable)
  colnames(millimeterTable)[1]<-"sample"
  fileName <- "millimeterTable.csv"
  write.csv2(millimeterTable,file=paste(OUTPUT_DIR,fileName,sep = "/"),row.names = FALSE)
}


writeSirAntibiotics <- function(includedMeasurementTable,headingToAntibioticTable,modelLimitTable)
{
  millimeterTable <- getMillimetertable(includedMeasurementTable)
  modes <- colnames(modelLimitTable[2:ncol(modelLimitTable)])
  zones <- millimeterTable[,headingToAntibioticTable$Heading]
  colnames(zones) <- headingToAntibioticTable$Antibiotic
  mode <- "Mode A"
  for(mode in modes){
      limit <- modelLimitTable[headingToAntibioticTable$Antibiotic==modelLimitTable$Antibiotic,mode]
      sir <-data.frame(zones)
      for(abindex in  1:ncol(zones)){
        sir[zones[,abindex]<limit[abindex],abindex] <- "R"
        sir[zones[,abindex]>=limit[abindex],abindex] <- "S"
      }
      fileName <- paste("sirAntibioticsModel","_",sub(" ", "-", mode),".csv",sep="")

      sir <- cbind(rownames(sir),sir)
      colnames(sir)[1]<-"sample"
      write.csv2(sir,file=paste(OUTPUT_DIR,fileName,sep = "/"),row.names = FALSE)
  }
}

writeDemographics <- function(includedMeasurementTable) {
  
  demographicsModel <- getDemographicsTable(includedMeasurementTable)
  write.csv2(demographicsModel,file=paste(OUTPUT_DIR,"demographicsModel.csv",sep = "/"),row.names = FALSE)
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

readMeasurementsTable <- function()
{
  phenotypeExcel =  paste(processedRootExcel, "ModelInput.xlsx", sep="/")
  sheetToDataFrame(read_xlsx(phenotypeExcel,sheet = SHEET_MEASUREMENTS_RAW))
}

writeSeparateInputTables <- function()
{
  phenotypeExcel =  paste(processedRootExcel, "ModelInput.xlsx", sep="/")
  
  measurementTable <- sheetToDataFrame(read_xlsx(phenotypeExcel,sheet = SHEET_MEASUREMENTS_RAW))
  modelLimitTable <- sheetToDataFrame(read_xlsx(phenotypeExcel,sheet = SHEET_MODEL_LIMIT))
  headingToAntibioticTable <- sheetToDataFrame(read_xlsx(phenotypeExcel,sheet = SHEET_HEADING_TO_ANTIBIOTIC))
  
  
  headingToAntibioticTable <- changeToAntibioticTypeOrder(headingToAntibioticTable,headingToAntibioticTable$Antibiotic,allAntibioticsInModelModelNames)
  modelLimitTable <- changeToAntibioticTypeOrder(modelLimitTable,modelLimitTable$Antibiotic,allAntibioticsInModelModelNames)
  
  includedMeasurementTable <- getIncludedSubset(measurementTable)

    
  writeSirAntibiotics(includedMeasurementTable,headingToAntibioticTable,modelLimitTable)
  writeMillimeterTable(includedMeasurementTable)
  writeDemographics(includedMeasurementTable)
  
}


writeCommonInputTables <- function()
{
  phenotypeExcel =  paste(processedRootExcel, "ModelInput.xlsx", sep="/")
  modelLimitTable <- sheetToDataFrame(read_xlsx(phenotypeExcel,sheet = SHEET_MODEL_LIMIT))

  modes <- colnames(modelLimitTable[2:ncol(modelLimitTable)])
  mode <- "Mode A"
  for(mode in modes){
    modelFile <- paste("sirAntibioticsModel","_",sub(" ", "-", mode),".csv",sep="")
    sirAntibioticsModel <- read.csv2(row.names=1,paste(modelDirectory,"input",modelFile,sep="/"))
    demographicsModel <- read.csv2(row.names=1,paste(modelDirectory,"input","demographicsModel.csv",sep="/"))
    antibioticsNames <- colnames(sirAntibioticsModel)
    inputAntibioticsNames <- antibioticsNames
    numberOfInputAntibiotics <- length(inputAntibioticsNames)
    sirDataFrame <- sirAntibioticsModel
    demographicsDataframe <- demographicsModel

    sirDataFrameWords <- makeWordsDataFrame(sirDataFrame,demographicsModel)
    csvTable <- sampleAsColumns(sirDataFrameWords)
    file = paste(modelDirectory,"input",paste("sirAntibioticsModelWords_",sub(" ", "-", mode),".csv",sep=""),  sep="/")
    write.csv2(x=csvTable,file=file,row.names = FALSE)
    
    demographicsModelNoDate <- demographicsModel
    demographicsModelNoDate$date <- "<unk>"
    sirDataFrameWords <- makeWordsDataFrame(sirDataFrame,demographicsModelNoDate)
    csvTable <- sampleAsColumns(sirDataFrameWords)
    file = paste(modelDirectory,"input",paste("sirAntibioticsModelWordsNoDate_",sub(" ", "-", mode),".csv",sep=""),  sep="/")
    write.csv2(x=csvTable,file=file,row.names = FALSE)
  }
}

countPredictableAntibiotics <- function()
{
  phenotypeExcel =  paste(processedRootExcel, "ModelInput.xlsx", sep="/")
  modelLimitTable <- sheetToDataFrame(read_xlsx(phenotypeExcel,sheet = SHEET_MODEL_LIMIT))
  
  modes <- colnames(modelLimitTable[2:ncol(modelLimitTable)])
  mode <- "Mode A"
  for(mode in modes){
    modelFile <- paste("sirAntibioticsModel","_",sub(" ", "-", mode),".csv",sep="")
    sirAntibioticsModel <- read.csv2(row.names=1,paste(modelDirectory,"input",modelFile,sep="/"))
    demographicsModel <- read.csv2(row.names=1,paste(modelDirectory,"input","demographicsModel.csv",sep="/"))
    antibioticsNames <- colnames(sirAntibioticsModel)
    inputAntibioticsNames <- antibioticsNames
    numberOfInputAntibiotics <- length(inputAntibioticsNames)

    sirDataFrame <- sirAntibioticsModel
    file = paste(modelDirectory,"input",paste("sirDataFrame_",sub(" ", "-", mode),".csv",sep=""),  sep="/")
    write.csv2(x=sirDataFrame,file=file,row.names = TRUE)

    R = apply(sirDataFrame, 2, \(x) sum(x=="R"))    
    S = apply(sirDataFrame, 2, \(x) sum(x=="S"))    
    df <- data.frame(S,R)
    file = paste(modelDirectory,"input",paste("sirStatsPerAntibiotic_",sub(" ", "-", mode),".csv",sep=""),  sep="/")
    write.csv2(x=df,file=file,row.names = TRUE)

    R = sum(sirDataFrame=="R")    
    S = sum(sirDataFrame=="S")    
    df <- data.frame(S,R)
    file = paste(modelDirectory,"input",paste("sirStats_",sub(" ", "-", mode),".csv",sep=""),  sep="/")
    write.csv2(x=df,file=file,row.names = TRUE)
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


# sirToWord <- function(sampleRow,abHeading)
# {
#   result <- ""
#   for(ab in colnames(sampleRow)){
#     if(str_detect(abHeading,ab)){
#       abString <- paste(ab,sampleRow[ab],sep="_")
#       if(result==""){
#         result <- abString 
#       } else {
#         result <- paste(result,abString,sep=" ")
#       }
#     }
#   }
#   result
# }



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


createOutputAndCompinationFiles <- function()
{
  phenotypeExcel =  paste(processedRootExcel, "ModelInput.xlsx", sep="/")
  modelLimitTable <- sheetToDataFrame(read_xlsx(phenotypeExcel,sheet = SHEET_MODEL_LIMIT))
  
  modes <- colnames(modelLimitTable[2:ncol(modelLimitTable)])
  mode <- modes[[1]] # Same for all moder
  modelFile <- paste("sirAntibioticsModel","_",sub(" ", "-", mode),".csv",sep="")
  sirAntibioticsModel <- read.csv2(row.names=1,paste(modelDirectory,"input",modelFile,sep="/"))
  antibioticsNames <- colnames(sirAntibioticsModel)
  inputAntibioticsNames <- antibioticsNames
  sirDataFrame <- sirAntibioticsModel

  for(k in 1:length(inputAntibioticsNames)){
    table <- createDataframeWithColumns(k,inputAntibioticsNames,sirDataFrame)
    emptyTable <- table
    emptyTable[] = ""
    csvTable <- sampleAsColumns(emptyTable)
    file = paste(modelDirectory,"input",paste("modelOutput",k,".csv",sep=""),  sep="/")
    write.csv2(x=csvTable,file=file,row.names = FALSE)
    
    csvTable <- createDataframeWithRows(k,inputAntibioticsNames)
    file = paste(modelDirectory,"input",paste("combinations",k,".csv",sep=""),  sep="/")
    write.table(x=csvTable,file=file,sep=";",row.names = FALSE, col.names = TRUE)
  }
  
}




ALL  <- function()
{
  if(!dir.exists(OUTPUT_DIR)){
    dir.create(OUTPUT_DIR)
  }
  writeSeparateInputTables()

  writeCommonInputTables()
  createOutputAndCompinationFiles()
}
