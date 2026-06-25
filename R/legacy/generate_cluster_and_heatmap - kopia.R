# #library(stats)
# library(readxl)
# library(dplyr)
# #library(tibble)
# library(pheatmap)
# #library(tidyverse) 
# #library(cluster)   
# library(factoextra)
# #library(fpc)
# #library(gtable)
# library(ggplot2)
# library(patchwork)
# 

library(readxl)
library(dplyr)
library(pheatmap)
library(factoextra)
library(ggplot2)
library(patchwork)

source(file='common - kopia.R')
source(file='manuscript/manuscriptcommon - kopia.R')
MODE <- "Mode-A"
source(file='model/modelcommon - kopia.R')

GENOTYPE_DIR <- paste(processedRootR,assemblymethod,"genotype",sep="//")

OUTPUT_DIR <- paste(processedRootR,"cluster",sep="/")
PDF_OUTPUT_DIR <- paste(OUTPUT_DIR,"pdf",sep="/")
TIFF_OUTPUT_DIR <- paste(OUTPUT_DIR,"tiff",sep="/")
PNG_OUTPUT_DIR <- paste(OUTPUT_DIR,"png",sep="/")

FONTSIZE_COL <- 10

#CLUSTER_LIST <- list()
#PHEATMAP_LIST <- list()

allAntibioticsWithZonesInputColumns <- c("KISS I-F",	"KISS I-MEL","KISS I-CFR",	"KISS I-W",	"KISS I-CIP",	"KISS I-AMC",	"KISS II-CTX",	"KISS II-CAZ"	,"KISS II-MEM",	"KISS II-TOB",	"KISS II-TZP",	"KISS II-SXT",	
                                         "DDT-FOX","DDT-FEP", "Studie-1-AMP"	,"Studie-1-PRL",	"Studie-1-CN"	,"Studie-1-CRO",	"Studie-2-LEV","Studie-2-MFX",	"Studie-2-OFX",	"Studie-2-NA")
allAntibioticsWithZonesInputNames <- c("F",	"MEL","CFR",	"W",	"CIP",	"AMC",	"CTX",	"CAZ"	,"MEM",	"TOB",	"TZP",	"SXT",	
                                       "FOX","FEP", "AMP"	,"PIP",	"GEN"	,"CRO",	"LVX","MFX",	"OFX",	"NAL")


allAntibioticsInModel <- c("AMP",	"AMC","PRL","TZP", "CAZ","CRO",	"CTX"	,"FEP", "CIP","OFX"	,"LEV","MFX","CN", "TOB" )

allAntibioticsInModelModelNames <- ALL_ANTIBIOTICS_IN_MODEL

allAntibioticsWithBreakpointNames <- c("F",	"MEL","CFR",	"W",	"CIP",	"AMC",	"CTX",	"CAZ"	,"MEM",	"TOB",	"TZP",	"SXT",	
                                       "FEP", "AMP"	,"PIP",	"GEN"	,"CRO",	"LVX","MFX",	"OFX")



allAntibioticsWithInvasiveSirInputColumns <- c("CIP",	"AMC",	"CTX",	"CAZ"	,"MEM",	"TOB",	"TZP",	"SXT",	
                                               "FEP", "AMP"	,"PRL",	"CN"	,"CRO",	"LEV","MFX",	"OFX")

allAntibioticsWithUTISirInputColumns <- c("F",	"MEL","CFR",	"W")


INPUT_ANTIBIOTICS <- 6

export_plot_bundle_no_excel <- function(plot,
                                        file_stub,
                                        width = 6.5,
                                        height = 4,
                                        dpi = 600,
                                        export = TRUE) {
  if (!isTRUE(export)) return(invisible(NULL))
  if (!dir.exists(PNG_OUTPUT_DIR)) dir.create(PNG_OUTPUT_DIR, recursive = TRUE)
  
  ggplot2::ggsave(file.path(PNG_OUTPUT_DIR, paste0(file_stub, ".png")),
                  plot, width = width, height = height, dpi = dpi, bg = "white")
  invisible(NULL)
}

wrap_pheatmap <- function(ph) {
  patchwork::wrap_elements(full = ph$gtable)
}

panel_box_theme <- theme(
  plot.background = element_rect(fill = NA, colour = NA),
  plot.margin = margin(4, 4, 4, 4)
)

readESBL <- function()
{
  studieNummerColumn <- c("Studienummer")
  phenotypeExcel =  paste(processedRootExcel, "Tolkad_data.xlsx", sep="/")
  zoneMillisSheet <- read_xlsx(phenotypeExcel,sheet = 'MeasureIncluded')
  
  sirModel <- readSirAntibioticsModel()
  
  esblSubset <- as.data.frame(
    zoneMillisSheet  %>% select(c(studieNummerColumn,"DDT","AMPC","KISS II-MEM")))
  sirSubset <-   as.data.frame(sirModel %>% select(c("CAZ","CTX")))
  
  rownames(esblSubset) <- esblSubset[,1]
  esblSubset <- esblSubset[-1]
  
  #  oldESBL <- esblSubset$DDT!="n" | esblSubset$AMPC!="n"
  #  length(which(oldESBL))
  
  isM <- (esblSubset$AMPC=="p") & (sirModel$CTX=="R" | sirModel$CAZ=="R")
  isA <- (esblSubset$DDT=="p" & esblSubset$AMPC!="p") & (sirModel$CTX=="R" | sirModel$CAZ=="R")
  isCARBA <- (esblSubset$`KISS II-MEM`<25) & (sirModel$CTX=="R" | sirModel$CAZ=="R")
  ESBL <- as.character(isA | isM | isCARBA)
  ESBL[ESBL=="FALSE"]="no"
  ESBL[isA]="A"
  ESBL[isM]="M"
  ESBL[isCARBA]="C"
  names(ESBL) <-rownames(esblSubset)
  
  ESBL  
}



readMillimetertable <- function()
{
  xxx <- read.csv2(paste(modelDirectory,"input","millimetertable.csv",sep="/"),check.names=FALSE)
  colnames(xxx)[-1] <- allAntibioticsWithZonesInputNames
  xxx
}

readPredictionsTable <- function(compare)
{
  #compare = "ME"
  name <- paste(tolower(compare),"RateSampleVsAntibiotic","-",INPUT_ANTIBIOTICS,sep="")
  readStatisticsExcel(name)
 }


readSirAntibioticsModel <- function()
{
  phenotypeExcel =  paste(processedRootExcel, "Tolkad_data.xlsx", sep="/")
  zoneMillisSheet <- read_xlsx(phenotypeExcel,sheet = 'MeasureIncluded')
  sirSheetInvasiv <- read_xlsx(phenotypeExcel,sheet = 'Tolka Invasiv')

  sirAntibioticsModel <- data.frame(sirSheetInvasiv[1:99,allAntibioticsInModel])
  colnames(sirAntibioticsModel) <- allAntibioticsInModelModelNames
  
  sirAntibioticsModel <- cbind(zoneMillisSheet$Studienummer,sirAntibioticsModel)
  colnames(sirAntibioticsModel)[1] <- "sample"
  
  sirAntibioticsModel
}


readATUAntibioticsModel <- function()
{
  phenotypeExcel =  paste(processedRootExcel, "Tolkad_data.xlsx", sep="/")
  zoneMillisSheet <- read_xlsx(phenotypeExcel,sheet = 'MeasureIncluded')
  atuSheet <- read_xlsx(phenotypeExcel,sheet = 'Tolka ATU')

  atuSheetsModel <- data.frame(atuSheet[1:99,allAntibioticsInModel])
  
  atuSheetsModel[,] <- !is.na(atuSheetsModel) & atuSheetsModel=="Ja"
  
  colnames(atuSheetsModel) <- allAntibioticsInModelModelNames
  rownames(atuSheetsModel) <- atuSheet$ATU
  atuSheetsModel
}




greenYellowRed <- function()
{
  # Define the colors to transition between
  colors <- c("red", "yellow", "green")
  
  # Create a color ramp palette function
  color_ramp <- colorRampPalette(colors)
  
  # Generate a sequence of colors
  n <- 100  # Number of steps in the gradient
  gradient_colors <- color_ramp(n)
  gradient_colors
}

green <- function()
{
  # Define the colors to transition between
  colors <- c("white", "green")
  
  # Create a color ramp palette function
  color_ramp <- colorRampPalette(colors)
  
  # Generate a sequence of colors
  n <- 100  # Number of steps in the gradient
  gradient_colors <- color_ramp(n)
  gradient_colors
}

red <- function()
{
  # Define the colors to transition between
  colors <- c("white", "red")
  
  # Create a color ramp palette function
  color_ramp <- colorRampPalette(colors)
  
  # Generate a sequence of colors
  n <- 100  # Number of steps in the gradient
  gradient_colors <- color_ramp(n)
  gradient_colors
}

orange <- function()
{
  # Define the colors to transition between
  colors <- c("white", "orange")
  
  # Create a color ramp palette function
  color_ramp <- colorRampPalette(colors)
  
  # Generate a sequence of colors
  n <- 100  # Number of steps in the gradient
  gradient_colors <- color_ramp(n)
  gradient_colors
}



clustersAndAnnotation <- function(dataFrame,numberOfClusters)
{
  data <- as.matrix(dataFrame)
  rownames(data) <- rownames(dataFrame)
  colnames(data) <- colnames(dataFrame)
  
  NO_CLUSTERS <- numberOfClusters
  
  # Perform hierarchical clustering
  dist_matrix <- dist(data)
  hc <- hclust(dist_matrix, method = "ward.D2")
  

  
  # Cut the dendrogram into clusters
  #  clusters <- cutree(hc2, k = NO_CLUSTERS)
  clusters <- cutree(hc, k = NO_CLUSTERS)
  
  # Create a data frame for annotation with correct factor levels
  annotation <- data.frame(Cluster = factor(clusters))
  rownames(annotation) <- rownames(data)
  
  # Ensure the levels of the Cluster factor match the keys in annotation_colors
  cluster_levels <- as.character(1:NO_CLUSTERS)
  annotation$Cluster <- factor(annotation$Cluster, levels = cluster_levels)
  
  # Define colors for each of the clusters
  annotation_colors <- list(Cluster = setNames(c("red", "blue", "green", "purple", "orange", "pink", "brown","cyan","yellow","gray")[1:NO_CLUSTERS], cluster_levels))
  list ("data"=data,"clusters" = clusters, "annotation" = annotation,"annotation_colors"=annotation_colors)
}

decorateWithESBL <- function(caa)
{
  annotation <- caa$annotation
  annotation_colors <- caa$annotation_colors
  
  ESBL <- as.factor(readESBL())
  
  annotation <- cbind(annotation,ESBL)
  whiteYellowOrangeRed <- setNames(c("yellow", "orange", "red","white"), c("A","M","C","no"))
  colnames(annotation)[ncol(annotation)] <- "ESBL(pheno)" 
  
  annotation_colors <- append(annotation_colors, list(whiteYellowOrangeRed))
  names(annotation_colors)[length(annotation_colors)]<- "ESBL(pheno)" 
  
  caa$annotation_colors <- annotation_colors
  caa$annotation <- annotation
  caa
}

decorateWithBushJacobiAMPAMC <- function(caa)
{
  bushJacobi <- read.csv2(check.names = FALSE,file= paste(sep="",GENOTYPE_DIR, "\\" ,assemblymethod,"_bushJacobi_",amrfinderDatabase, ".csv"))
  
  annotation <- caa$annotation
  annotation_colors <- caa$annotation_colors
  
  SS_NAME <- "(none)"
  SR_NAME <- "2b"
  RR_NAME <- "1/2be/2br/2d/3"
  
  SS <- !(apply(bushJacobi[c("1","2b","2be","2br","2d","3")],1,any))
  RR <-  apply(bushJacobi[c("1","2be","2br","2d","3")],1,any)
  SR <- !(SS | RR)
  BushJacobi <- as.vector(SS)
  BushJacobi[which(SS)] <- SS_NAME
  BushJacobi[which(RR)] <- RR_NAME
  BushJacobi[which(SR)] <- SR_NAME
  BushJacobi <- as.factor(BushJacobi)
  
  whiteYellowRed <- setNames(c("white", "yellow","red"), c(SS_NAME,SR_NAME,RR_NAME))
  annotation <- cbind(annotation,BushJacobi)
  annotation_colors <- append(annotation_colors, list(whiteYellowRed))

  name <- "Bush Jacobi"
  colnames(annotation)[ncol(annotation)] <- name
  names(annotation_colors)[length(annotation_colors)]<- name 
  caa$annotation_colors <- annotation_colors
  caa$annotation <- annotation
  caa
}

decorateWithBushJacobiESBL <- function(caa)
{
  bushJacobi <- read.csv2(check.names = FALSE,file= paste(sep="",GENOTYPE_DIR, "\\" ,assemblymethod,"_bushJacobi_",amrfinderDatabase, ".csv"))
  
  annotation <- caa$annotation
  annotation_colors <- caa$annotation_colors
  name <-"ESBL(geno)" 
  
  NONE_NAME <- "(none)/2b/2br/2d/"
  ESBL_A_NAME <- "2be"
  ESBL_M_NAME <- "1"
  ESBL_CARBA_NAME <-"3"
  
  xxx <- which(apply(bushJacobi["2d"],1,any))
  yyy <- which(apply(bushJacobi[c("1","2be","3")],1,any))
  intersect(xxx,yyy)
    
  NONE <- !(apply(bushJacobi[c("1","2be","3")],1,any))
  ESBL_A <- !(apply(bushJacobi[c("1","3")],1,any)) & !NONE
  ESBL_M <-  !(apply(bushJacobi[c("3")],1,any)) & !(ESBL_A | NONE)
  ESBL_CARBA <- !(NONE | ESBL_A | ESBL_M)
  
  ESBL_geno <- as.vector(NONE)
  ESBL_geno[which(NONE)] <- NONE_NAME
  ESBL_geno[which(ESBL_A)] <- ESBL_A_NAME
  ESBL_geno[which(ESBL_M)] <- ESBL_M_NAME
  ESBL_geno[which(ESBL_CARBA)] <- ESBL_CARBA_NAME
  ESBL_geno <- as.factor(ESBL_geno)
  
  whiteYellowOrangeRed <- setNames(c("white", "yellow","orange","red"), c(NONE_NAME,ESBL_A_NAME,ESBL_M_NAME,ESBL_CARBA_NAME))
  annotation <- cbind(annotation,ESBL_geno)
  annotation_colors <- append(annotation_colors, list(whiteYellowOrangeRed))
  
  colnames(annotation)[ncol(annotation)] <- name
  names(annotation_colors)[length(annotation_colors)]<- name 
  
  caa$annotation_colors <- annotation_colors
  caa$annotation <- annotation
  caa
}


decorateWithEnzymeFamily <- function(caa)
{
  enzymeFamily <- read.csv2(check.names = FALSE,file= paste(sep="",GENOTYPE_DIR, "\\" ,assemblymethod,"_enzymeFamily_",amrfinderDatabase, ".csv"))

  typeof(enzymeFamily)
  ef = enzymeFamily
  ef <- lapply(ef, function(x) replace(x, x == "TRUE", "Yes"))
  ef <- lapply(ef, function(x) replace(x, x == "FALSE", "No"))
  typeof(ef)
  enzymeFamily <- ef
  
  annotation <- caa$annotation
  annotation_colors <- caa$annotation_colors
  
  orangeAndWhite = setNames(c("orange", "white"), c("Yes","No"))
  redAndWhite = setNames(c("red", "white"), c("Yes","No"))
  
  classes <- names(enzymeFamily[-1])
  classes <- rev(classes)
  for(class in classes) {
    # class <- "AMPC"
    # classes <- AMPC, CTX-M, LAP, NDM, OXA, SHV, TEM
    classVector <- enzymeFamily[[class]]
    names(classVector) <- enzymeFamily[[1]]
    classFactor <- as.factor(classVector) 
    
    annotation <- cbind(annotation,classFactor)
    colnames(annotation)[ncol(annotation)] <- class
    
    if(class %in% c("AMPC","NDM","OXA")) {
      annotation_colors <- append(annotation_colors, list(redAndWhite))
    } else {
      annotation_colors <- append(annotation_colors, list(orangeAndWhite))
    }
    name <- class
    colnames(annotation)[ncol(annotation)] <- name
    names(annotation_colors)[length(annotation_colors)]<- name 
    
  }
  caa$annotation_colors <- annotation_colors
  caa$annotation <- annotation
  caa
}




decorateWithBushJacobi <- function(caa)
{
  bushJacobi <- read.csv2(check.names = FALSE,file= paste(sep="",GENOTYPE_DIR, "\\" ,assemblymethod,"_bushJacobi_",amrfinderDatabase, ".csv"))
  
  annotation <- caa$annotation
  annotation_colors <- caa$annotation_colors
  
  yelowAndWhite = setNames(c("yellow", "white"), c("TRUE","FALSE"))
  redAndWhite = setNames(c("red", "white"), c("TRUE","FALSE"))
  
  classes <- colnames(bushJacobi[-1])
  classes <- rev(classes)
  for(class in classes) {  
    classVector <- bushJacobi[[class]]
    names(classVector) <- bushJacobi[[1]]
    classFactor <- as.factor(classVector) 
    
    annotation <- cbind(annotation,classFactor)
    colnames(annotation)[ncol(annotation)] <- class
    
    if(class=="2b") {
      annotation_colors <- append(annotation_colors, list(yelowAndWhite))
    } else {
      annotation_colors <- append(annotation_colors, list(redAndWhite))
    }
    #names(annotation_colors)[length(annotation_colors)]<- class
    
    name <- paste("BJ",class)
    colnames(annotation)[ncol(annotation)] <- name
    names(annotation_colors)[length(annotation_colors)]<- name 
    
  }
  caa$annotation_colors <- annotation_colors
  caa$annotation <- annotation
  caa
}

exampleAnnotate <- function()
{
  # Load necessary libraries
  library(pheatmap)
  library(factoextra)
  
  # Generate some example data
  set.seed(123)
  data <- matrix(rnorm(100), nrow = 10)
  rownames(data) <- paste0("row_", seq(nrow(data)))
  colnames(data) <- paste0("row_", seq(ncol(data)))
  
  # Perform hierarchical clustering
  dist_matrix <- dist(data)
  hc <- hclust(dist_matrix,method = "ward.D2")
  
  # Cut the dendrogram into 7 clusters
  clusters <- cutree(hc, k = 7)
  
  # Create a data frame for annotation with correct factor levels
  annotation <- data.frame(Cluster = factor(clusters))
  rownames(annotation) <- rownames(data)
  
  # Ensure the levels of the Cluster factor match the keys in annotation_colors
  cluster_levels <- as.character(1:7)
  annotation$Cluster <- factor(annotation$Cluster, levels = cluster_levels)
  
  # Define colors for each of the 7 clusters
  annotation_colors <- list(Cluster = setNames(c("red", "blue", "green", "purple", "orange", "pink", "brown"), cluster_levels))
  
  # Create heatmap with annotation
  pheatmap(data, annotation_row = annotation, annotation_colors = annotation_colors)
  
}




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




makeSirForHeatmap <- function(sirAntibiotics,modelMillimeterTable)
{
  sirToClusterOn <- sirAntibiotics[-1]
  sirToClusterOn[sirToClusterOn=="S"] <- 1
  sirToClusterOn[sirToClusterOn=="I"] <- 0.5
  sirToClusterOn[sirToClusterOn=="R"] <- 0
  
  # Make a table 
  
  for (row in 1:nrow(modelMillimeterTable)){
    for (col in 1:ncol(modelMillimeterTable)){
      modelMillimeterTable[row,col] = as.numeric(sirToClusterOn[row,col])
    }  
  }
  modelMillimeterTable
}

makeClusterPlot <- function(name,caa)
{
  # Visualize the clusters with fviz_cluster (which will perform PCA by default)
  p <- fviz_cluster(list(data = caa$data, cluster = caa$clusters),
                    palette = unlist(caa$annotation_colors$Cluster),
                    ellipse.type = "convex",  # Type of cluster ellipse
                    geom = "point",           # Point geometry
                    show.clust.cent = TRUE) + ggtitle(label='')  # Show cluster centers
  dumpFVizObject(name,p)
  p
}



fetchTables<- function(){
  millimeterTable <- readMillimetertable()
  rownames(millimeterTable) <- millimeterTable$sample
  modelMillimeterTable <- millimeterTable %>% select(all_of(allAntibioticsInModelModelNames))
  
  sirAntibiotics <-  readSirAntibioticsModel()
  
  sirTable <- makeSirForHeatmap(sirAntibiotics,modelMillimeterTable)
  rescaledMillimeters <- rescaleMillimeters(modelMillimeterTable) 
  millimetersTable <- modelMillimeterTable
  
  tablesList <- list (sirTable, rescaledMillimeters, millimetersTable )
  names(tablesList) <- c("SIR","RESCALEDMM","MM")
  tablesList
}

ESTIMATE_CLUSTERS <- function()
{
  library(cluster)
  tablesList <- fetchTables()
  for(tableName in names(tablesList) ){
    print(tableName)
    tableToClusterOn <- tablesList[[tableName]]
    estimateOptimalNumberofclusters(paste("estimate",tableName,sep = "_"),tableToClusterOn)
  }
}

fetchOptimalClusterSizeList <-function() {
  optimalClusterList <- list(10,8,4)
  names(optimalClusterList) <- c("SIR","RESCALEDMM","MM")
  optimalClusterList
}


fetchPredictionsTable <- function(compare)
{
  table <- readPredictionsTable(compare)
  sampleId <- table$sample
  table <- as.data.frame(table %>% select(all_of(allAntibioticsInModelModelNames)))
  rownames(table) <- sampleId
  table
}



PHEATMAPS_AND_CLUSTERS_SELECTED <- function()
{
  tablesList <- fetchTables()
  primaryTableName <- "RESCALEDMM"
  tableToClusterOn <- tablesList[[primaryTableName]]
  numberOfClusters <- fetchOptimalClusterSizeList()[[primaryTableName]]
  baseCaa <- clustersAndAnnotation(tableToClusterOn,numberOfClusters)

  
  FONT_SIZE <- 3
  SHOW_NAMES <- FALSE
  COLOR <- greenYellowRed()
  LEGEND = FALSE
  CLUSTER_ROWS = TRUE
  CLUSTER_COLS = TRUE
  row_order <- NA
  col_order <- NA
  
  
  ##############################  ONE cluster and pheatmap ###################################  
  #Primary cluster and pheatmap as one picture in a grid  
  cpName <- paste("cluster",primaryTableName,numberOfClusters,sep = "_")
  cp <- makeClusterPlot(cpName,baseCaa)  
  
  caa <- baseCaa
  ANNOTATION_ROW <- caa$annotation
  ANNOTATION_COLORS <- caa$annotation_colors
  
  COLOR <- rev(red())
  LEGEND = FALSE
  
  
  name <- paste("pheatmapWithCluster",primaryTableName,numberOfClusters,sep = "_")
  p <- pheatmap(tableToClusterOn,
                clustering_method = "ward.D2",
                fontsize_row = FONT_SIZE,
                fontsize_col = FONTSIZE_COL,
                show_rownames = SHOW_NAMES,
                color=COLOR,
                cluster_rows = CLUSTER_ROWS,
                cluster_cols = CLUSTER_COLS, 
                annotation_row = ANNOTATION_ROW, 
                annotation_colors =  
                ANNOTATION_COLORS,
                border_color = NA,
                legend=LEGEND)
  row_order <-p$tree_row$order
  col_order <-p$tree_col$order

  # cp is ggplot, p is pheatmap object
  pA <- cp + panel_box_theme
  pB <- wrap_pheatmap(p) & panel_box_theme   # `&` applies theme to wrapped element
  
  combined <- (pA + pB) +
    plot_layout(ncol = 2) +
    plot_annotation(tag_levels = "A")  # gives A, B automatically
  
  combined
  export_plot_bundle_no_excel(plot = combined,file_stub ="ClusterAndHeatmap",width = 12 ,height = 6, export = TRUE )

  #From now on we skip clustering as row_order and col_orders are fixed
  CLUSTER_ROWS = FALSE
  CLUSTER_COLS = FALSE

  ##############################  TWO SIR  ###################################  

  # VME and MEs as one picture in a grid
  caa <- baseCaa
  ANNOTATION_ROW <- caa$annotation
  ANNOTATION_COLORS <- caa$annotation_colors
  LEGEND <- FALSE
  
  predictions <- fetchPredictionsTable("VME")
  COLOR <-  red()
  vmePheatmap <- pheatmap(predictions[row_order,col_order],clustering_method = "ward.D2",fontsize_row = FONT_SIZE,fontsize_col = FONTSIZE_COL,show_rownames = SHOW_NAMES,color=COLOR,cluster_rows = CLUSTER_ROWS,cluster_cols = CLUSTER_COLS,annotation_row = ANNOTATION_ROW,annotation_legend=FALSE, annotation_colors =  ANNOTATION_COLORS,border_color = NA,legend=FALSE,main = "")
  
  predictions <- fetchPredictionsTable("ME")
  COLOR <-  orange()
  mePheatmap <- pheatmap(predictions[row_order,col_order],
                         clustering_method = "ward.D2",
                         fontsize_row = FONT_SIZE,
                         fontsize_col = FONTSIZE_COL,
                         show_rownames = SHOW_NAMES,
                         color= COLOR,
                         cluster_rows = CLUSTER_ROWS,
                         cluster_cols = CLUSTER_COLS,
                         annotation_row = ANNOTATION_ROW, 
                         annotation_legend = FALSE,
                         annotation_colors =ANNOTATION_COLORS,
                         legend=FALSE,
                         border_color = NA,
                         main = "")
  
  sirTable <- tablesList[["SIR"]]
  atuTable <- readATUAntibioticsModel()
  display_numbers <- matrix(ifelse(as.matrix(atuTable), "\u2217", ""), nrow(atuTable))
  
  COLOR <- greenYellowRed()
  sirPheatmap <- pheatmap(sirTable[row_order,col_order],
                          clustering_method = "ward.D2",
                          fontsize_row = FONT_SIZE,
                          fontsize_col = FONTSIZE_COL,
                          show_rownames = SHOW_NAMES,
                          color=COLOR,
                          cluster_rows = CLUSTER_ROWS,
                          cluster_cols = CLUSTER_COLS,
                          annotation_row = ANNOTATION_ROW, 
                          annotation_colors =  ANNOTATION_COLORS,
                          border_color = NA,
                          legend=FALSE,
                          main = "", 
                          display_numbers = display_numbers[row_order,col_order],
                          fontsize_number = 8)
  

  pA <- wrap_pheatmap(vmePheatmap) & panel_box_theme
  pB <- wrap_pheatmap(mePheatmap) & panel_box_theme
  pC <- wrap_pheatmap(sirPheatmap) & panel_box_theme
    
  combined <- (pA + pB +pC) +
    plot_layout(
      ncol = 3,
      widths = c(200, 200, 280),   # relative widths
      heights = c(480)
    ) + 
    plot_annotation(tag_levels = "A")  # gives A, B automatically
  
  combined
  export_plot_bundle_no_excel(plot = combined,file_stub ="HeatmapErrorSIR",width = 12 ,height = 6, export = TRUE )
  
  
  
  ##############################  THREE SIR  ###################################  
  caa <- baseCaa
  # remove cluster
  #caa$annotation <- caa$annotation[-1]
  #caa$annotation_colors <- caa$annotation_colors[-1]
  
  ANNOTATION_ROW <- caa$annotation
  ANNOTATION_COLORS <- caa$annotation_colors

  onlyBetalactams <- function(table){
    #FIXME
    table[,1:8]
    #table
  }
  
  selectedPrimaryTable <- onlyBetalactams(tablesList[[primaryTableName]])
  selectedSirTable <- onlyBetalactams(tablesList[["SIR"]])
  selectedAtuTable <- onlyBetalactams(readATUAntibioticsModel())
  selectedDisplay_numbers <- matrix(ifelse(as.matrix(selectedAtuTable), "\u2217", ""), nrow(selectedAtuTable))
  
  LEGEND <- FALSE
  COLOR <- greenYellowRed()
  selectedsirPheatmap <- pheatmap(selectedSirTable[row_order,],
                                  clustering_method = "ward.D2",
                                  fontsize_row = FONT_SIZE,
                                  fontsize_col = FONTSIZE_COL,
                                  show_rownames = SHOW_NAMES,
                                  color=COLOR,
                                  cluster_rows = CLUSTER_ROWS,
                                  cluster_cols = CLUSTER_COLS,
                                  annotation_row = ANNOTATION_ROW, 
                                  annotation_colors =  ANNOTATION_COLORS,
                                  border_color = NA,
                                  annotation_legend=TRUE,
                                  legend=FALSE,
                                  main = "", 
                                  display_numbers = selectedDisplay_numbers[row_order,],
                                  fontsize_number = 8)
  caa <- baseCaa
  # remove cluster
  caa$annotation <- caa$annotation[-1]
  caa$annotation_colors <- caa$annotation_colors[-1]
  caa <- decorateWithEnzymeFamily(caa)
  
  ANNOTATION_ROW <- caa$annotation
  ANNOTATION_COLORS <- caa$annotation_colors

  LEGEND <- FALSE
  COLOR <- rev(red())
  
  selectedPrimaryPheatmap <- pheatmap(selectedPrimaryTable[row_order,],
                                  clustering_method = "ward.D2",
                                  fontsize_row = FONT_SIZE,
                                  fontsize_col = FONTSIZE_COL,
                                  show_rownames = SHOW_NAMES,
                                  color=COLOR,
                                  cluster_rows = CLUSTER_ROWS,
                                  cluster_cols = CLUSTER_COLS,
                                  annotation_row = ANNOTATION_ROW, 
                                  annotation_colors =  ANNOTATION_COLORS,
                                  annotation_legend= TRUE, 
                                  #annotation_legend= list(TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),
                                  legend=FALSE,
                                  main = "", 
                                  fontsize_number = 8)
  
  
  pA <- wrap_pheatmap(selectedPrimaryPheatmap) & panel_box_theme
  pB <- wrap_pheatmap(selectedsirPheatmap) & panel_box_theme

  combined <- (pA + pB) +
    plot_layout(
      ncol = 2,
      widths = c(260,220),
      heights = c(480)
    ) +
  plot_annotation(tag_levels = "A")  # gives A, B automatically
  
  combined
  export_plot_bundle_no_excel(plot = combined,file_stub ="HeatmapBetalactamase",width = 12 ,height = 6, export = TRUE )
}


DUMP_MODEL_MM_FRAME <-function()
{
  frame <- fetchTables()[["MM"]]
  sample <- row.names(frame)
  frame <- cbind(sample,frame)
  # colnames(frame)
  write.csv2(frame,file.path(getCommonModelFolder(),"modelzonemillimeters.csv"),row.names = FALSE)
}

ALL  <- function()
{
  
  if(!dir.exists(PDF_OUTPUT_DIR)){
    dir.create(PDF_OUTPUT_DIR)
  }
  if(!dir.exists(TIFF_OUTPUT_DIR)){
    dir.create(TIFF_OUTPUT_DIR)
  }
  if(!dir.exists(PNG_OUTPUT_DIR)){
    dir.create(PNG_OUTPUT_DIR)
  }
  ESTIMATE_CLUSTERS()
  PHEATMAPS_AND_CLUSTERS_SELECTED()
  DUMP_MODEL_MM_FRAME()
  
}



dumpFVizObject <-function(name,p,...){
  pdf(paste(PDF_OUTPUT_DIR,paste(name,"pdf",sep="."),sep="/"))
  plot(p)
  dev.off()
  tiff(paste(TIFF_OUTPUT_DIR,paste(name,"tiff",sep="."),sep="/"),...)
  plot(p)
  dev.off()
}

estimateOptimalNumberofclusters <- function(name , dataMatrix)
{
  paste("gap",name,sep = "_")
  
  df <- scale(dataMatrix)
  df <- dataMatrix
  
  # Elbow method
  p = fviz_nbclust(df, hcut, method = "wss",hc_method = "ward.D2") +
    labs(subtitle = "Elbow method")
  dumpFVizObject(paste("wss",name,sep = "_"),p)
  
  # Silhouette method
  p = fviz_nbclust(df, hcut, method = "silhouette",hc_method = "ward.D2")+
    labs(subtitle = "Silhouette method")
  dumpFVizObject(paste("silhouette",name,sep = "_"),p)
  
  # Gap statistic
  # nboot = 50 to keep the function speedy. 
  # recommended value: nboot= 500 for your analysis.
  # Use verbose = FALSE to hide computing progression.
  set.seed(123)
  p = fviz_nbclust(df, hcut, nstart = 25,  method = "gap_stat", nboot = 50,hc_method = "ward.D2")+
    labs(subtitle = "Gap statistic method")
  dumpFVizObject(paste("gap_stat",name,sep = "_"),p)
  
  # Gap statistic
  set.seed(123)
  gap_stat <- clusGap(df, FUN = hcut, nstart = 25,
                      K.max = 10, B = 10)
  print(gap_stat, method = "firstmax")
  
  p = fviz_gap_stat(gap_stat)
  dumpFVizObject(paste("gap_stat_clusGap",name,sep = "_"),p)
  
}

