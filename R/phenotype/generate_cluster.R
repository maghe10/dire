library(stats)
library(readxl)
library(dplyr)
library(tibble)
library(pheatmap)
library(tidyverse) 
library(cluster)   
library(factoextra)
source(file='common.R')


#setwd("C:\\Users\\magnu\\OneDrive - Västra Götalandsregionen\\Documents\\KMIK\\Forskning\\CHAIR\\DIRE\\Resultat\\Arbetsmaterial\\R")

dir = processedRootExcel

allAntibioticsWithZonesInputColumns <- c("KISS I-F",	"KISS I-MEL","KISS I-CFR",	"KISS I-W",	"KISS I-CIP",	"KISS I-AMC",	"KISS II-CTX",	"KISS II-CAZ"	,"KISS II-MEM",	"KISS II-TOB",	"KISS II-TZP",	"KISS II-SXT",	
                                         "DDT-FOX","DDT-FEP", "Studie-1-AMP"	,"Studie-1-PRL",	"Studie-1-CN"	,"Studie-1-CRO",	"Studie-2-LEV","Studie-2-MFX",	"Studie-2-OFX",	"Studie-2-NA")

allAntibioticsInModel <- c("AMP",	"AMC","PRL","TZP", "CAZ","CRO",	"CTX"	,"FEP", "CIP","OFX"	,"LEV","MFX","CN", "TOB" )

allAntibioticsInModelModelNames <- c("AMP","AMC","PIP","TZP","CAZ","CRO","CTX","FEP","CIP","OFX","LVX","MFX","GEN","TOB")



allAntibioticsWithInvasiveSirInputColumns <- c("CIP",	"AMC",	"CTX",	"CAZ"	,"MEM",	"TOB",	"TZP",	"SXT",	
                                       "FEP", "AMP"	,"PRL",	"CN"	,"CRO",	"LEV","MFX",	"OFX")

allAntibioticsWithUTISirInputColumns <- c("F",	"MEL","CFR",	"W")

studieNummerColumn <- c("Studienummer")
demographicsInputColumns <-c("Datum", "Kön",	"Ålder")


phenotypeExcel =  paste(processedRootExcel, "Tolkad_data.xlsx", sep="/")
zoneMillisSheet <- read_xlsx(phenotypeExcel,sheet = 'MeasureIncluded')
sirSheetInvasiv <- read_xlsx(phenotypeExcel,sheet = 'Tolka Invasiv')
sirSheetUTI <- read_xlsx(phenotypeExcel,sheet = 'Tolka UTI')
measurementsRawSheet <- read_xlsx(phenotypeExcel,sheet = 'MeasurementsRaw')



readMillimetertable <- function()
{
  millimeterTable <- data.frame(zoneMillisSheet[1:99,c(allAntibioticsWithZonesInputColumns)])
  row.names(millimeterTable) <- zoneMillisSheet$Studienummer
  
  g = regexpr("\\.[^\\.]*$",colnames(millimeterTable))
  newNames <- substring(colnames(millimeterTable),g+1)
  colnames(millimeterTable) <- newNames
  

  millimeterTable[millimeterTable == "X"] <- NA
  for (col in colnames(millimeterTable)){
    millimeterTable[[col]] <- as.integer(millimeterTable[[col]])
  }
  millimeterTable
}

printAntibiotics <-function()
{
  modelInputAntibiotics4 <- c("AMP","CTX","CIP","TOB")
  modelInputColumns <- c(studieNummerColumn,demographicsInputColumns,modelInputAntibiotics4)
  modelInputTable <- cbind(zoneMillisSheet[,c(studieNummerColumn,demographicsInputColumns)], data.frame(sirSheet[1:100,modelInputAntibiotics4]))
  write.csv2(modelInputTable,file="modelInputAntibiotics4.csv",row.names = FALSE)
}



millimeterTable <- readMillimetertable()
modelMillimeterTable <- millimeterTable %>% select(allAntibioticsInModel)

g = regexpr("\\.[^\\.]*$",colnames(millimeterTable))
newNames <- substring(colnames(millimeterTable),g+1)
colnames(millimeterTable) <- newNames


write.csv2(millimeterTable,file="millimeterTable.csv",row.names = TRUE)



clusterOn  <- function(millimeterTable,centers) 
{
  df <- na.omit(millimeterTable)
  df <- scale(df)
  k <- kmeans(df, centers, nstart = 25)
  k
}


sirAntibiotics <-  cbind(data.frame(sirSheetInvasiv[1:99,allAntibioticsWithInvasiveSirInputColumns]),
                         data.frame(sirSheetUTI[1:99,allAntibioticsWithUTISirInputColumns]))

sirAntibioticsModel <- data.frame(sirSheetInvasiv[1:99,allAntibioticsInModel])
colnames(sirAntibioticsModel) <- allAntibioticsInModelModelNames
sirAntibioticsModel <- cbind(zoneMillisSheet$Studienummer,sirAntibioticsModel)
colnames(sirAntibioticsModel)[1] <- "sample"
#write.csv2(sirAntibioticsModel,file=paste(modelDirectory,"sirAntibioticsModel.csv",sep = "/"),row.names = FALSE)



measurementsRawSheetIncluded <- measurementsRawSheet[which(measurementsRawSheet$Exkludera %in% c(NA,"Nej")),]
demographicsModel <-  data.frame(measurementsRawSheetIncluded[,demographicsInputColumns])
demographicsModel$Datum <- as.Date(as.integer(demographicsModel$Datum),origin=as.Date("1900-01-01"))
demographicsModel$Kön[which(demographicsModel$Kön=="K")] <- "F"

demographicsModel  <- cbind(measurementsRawSheetIncluded$Studienummer,demographicsModel)
colnames(demographicsModel) <- c("sample","date","sex","age")
#write.csv2(demographicsModel,file=paste(modelDirectory,"demographicsModel.csv",sep = "/"),row.names = FALSE)







row.names(sirAntibiotics) <- zoneMillisSheet$Studienummer

k5 <- clusterOn(millimeterTable,5)
k10 <- clusterOn(millimeterTable,10)
k20 <- clusterOn(millimeterTable,20)
k30 <-clusterOn(millimeterTable,30)

k5model <- clusterOn(modelMillimeterTable,5)
k10model <- clusterOn(modelMillimeterTable,10)
k20model <- clusterOn(modelMillimeterTable,20)
k30model <-clusterOn(modelMillimeterTable,30)




clusteredSirAntibiotics <- sirAntibiotics
clusteredSirAntibiotics <- clusteredSirAntibiotics %>% add_column(kluster5=0)
clusteredSirAntibiotics[names(k5$cluster),'kluster5'] <- as.integer(k5$cluster)
clusteredSirAntibiotics <- clusteredSirAntibiotics %>% add_column(kluster10=0)
clusteredSirAntibiotics[names(k10$cluster),'kluster10'] <- as.integer(k10$cluster)
clusteredSirAntibiotics <- clusteredSirAntibiotics %>% add_column(kluster20=0)
clusteredSirAntibiotics[names(k20$cluster),'kluster20'] <- as.integer(k20$cluster)
clusteredSirAntibiotics <- clusteredSirAntibiotics %>% add_column(kluster30=0)
clusteredSirAntibiotics[names(k30$cluster),'kluster30'] <- as.integer(k30$cluster)

clusteredSirAntibiotics <- clusteredSirAntibiotics %>% add_column(kluster5model=0)
clusteredSirAntibiotics[names(k5model$cluster),'kluster5model'] <- as.integer(k5model$cluster)
clusteredSirAntibiotics <- clusteredSirAntibiotics %>% add_column(kluster10model=0)
clusteredSirAntibiotics[names(k5model$cluster),'kluster10model'] <- as.integer(k10model$cluster)
clusteredSirAntibiotics <- clusteredSirAntibiotics %>% add_column(kluster20model=0)
clusteredSirAntibiotics[names(k20model$cluster),'kluster20model'] <- as.integer(k20model$cluster)
clusteredSirAntibiotics <- clusteredSirAntibiotics %>% add_column(kluster30model=0)
clusteredSirAntibiotics[names(k30model$cluster),'kluster30model'] <- as.integer(k30model$cluster)

write.csv2(clusteredSirAntibiotics,file="clusteredSirAntibiotics.csv",row.names = TRUE)


pdf("k5.pdf")
fviz_cluster(k5,na.omit(millimeterTable))
dev.off()

pdf("k5model.pdf")
fviz_cluster(k5model,na.omit(modelMillimeterTable))
dev.off()

pdf("k10.pdf")
fviz_cluster(k10,na.omit(millimeterTable))
dev.off()

pdf("k10model.pdf")
fviz_cluster(k10model,na.omit(modelMillimeterTable))
dev.off()


pdf("k20.pdf")
fviz_cluster(k20,na.omit(millimeterTable))
dev.off()

pdf("k20model.pdf")
fviz_cluster(k20model,na.omit(modelMillimeterTable))
dev.off()


pdf("k30.pdf")
fviz_cluster(k30,na.omit(millimeterTable))
dev.off()

pdf("k30model.pdf")
fviz_cluster(k30model,na.omit(modelMillimeterTable))
dev.off()

FONT_SIZE <- 3
show_rownames <- FALSE
tiff("pheatmapAll.tiff", units="in", width=5, height=5, res=300)
pheatmap(millimeterTable,clustering_method = "ward.D2",fontsize_row = FONT_SIZE,show_rownames = show_rownames  )
dev.off()


pdf("pheatmapAll.pdf")
pheatmap(millimeterTable,clustering_method = "ward.D2",fontsize_row = FONT_SIZE )
dev.off()

jpeg("pheatmapAll.jpg")
pheatmap(millimeterTable,clustering_method = "ward.D2",fontsize_row = FONT_SIZE )
dev.off()



pdf("pheatmapAll_scaled.pdf")
pheatmap(scale(millimeterTable),clustering_method = "ward.D2",fontsize_row = FONT_SIZE )
dev.off()

pdf("pheatmapModel.pdf")
pheatmap(modelMillimeterTable,clustering_method = "ward.D2",fontsize_row = FONT_SIZE )
dev.off()

tiff("pheatmapModel.tiff", units="in", width=5, height=5, res=300)
pheatmap(modelMillimeterTable,clustering_method = "ward.D2",fontsize_row = FONT_SIZE,show_rownames = show_rownames )
dev.off()



pdf("pheatmapModel_scaled.pdf")
pheatmap(scale(modelMillimeterTable),clustering_method = "ward.D2",fontsize_row = FONT_SIZE )
dev.off()


print("done")

  

