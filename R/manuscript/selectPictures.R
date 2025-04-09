source(file = 'common.R')

OUTPUT_DIR <- paste(processedRootR,"cluster",sep="/")
PDF_OUTPUT_DIR <- paste(OUTPUT_DIR,"pdf",sep="/")
TIFF_OUTPUT_DIR <- paste(OUTPUT_DIR,"tiff",sep="/")


selectHeatmapsAndClusters <- function()
{
  files <- c("ClusterAndHeatmap.tiff",
             "HeatmapErrorSIR.tiff",
             "HeatmapBetalactamase.tiff")
  files <- paste(TIFF_OUTPUT_DIR,files,sep="/")
  file.copy(files, manuscriptDirectory,overwrite = TRUE)
}

ALL <-function()
{
  selectHeatmapsAndClusters()
}