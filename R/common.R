library(readr)

oneDriveRoot <- paste(Sys.getenv("USERPROFILE"),"OneDrive - Västra Götalandsregionen",sep="\\")
direRoot <- paste(oneDriveRoot,"DIRE",sep="\\")


# Different directories for code

rRoot <- paste(oneDriveRoot,"git","dire","R",sep="\\")

workingDirectory <- rRoot


assemblymethod <- 'spades_standard'
#assemblymethod <- 'unicycler_normal'
#assemblymethod <- 'spades_optimized'
#assemblymethod <- 'unicycler_conservative'

assemblyDirectory <- paste(direRoot,assemblymethod,"assembly",sep="/")
#amrfinderDatabase <- "231115.1"
amrfinderDatabase <- "2024-12-18.1"
amrfinderDirectory <- paste(direRoot,assemblymethod,"amrfinder",amrfinderDatabase,sep="/")

resfinderDatabase <- "v460"
resfinderDirectory <- paste(direRoot,assemblymethod,"resfinder",resfinderDatabase,sep="/")

qualityDirectory <- paste(direRoot,assemblymethod,"quality",sep="/")
tygsDirectory <- paste(qualityDirectory,"TYGS",sep="/")
jspecieswsDirectory <- paste(qualityDirectory,"jspeciesws",sep="/")
quastDirectory <- paste(qualityDirectory,"multiqc",sep="/")
confindrtrimmedDirectory <- paste(qualityDirectory,"confindrtrimmed",sep="/")
confindrrawDirectory <- paste(qualityDirectory,"confindrraw",sep="/")
checkmDirectory <- paste(qualityDirectory,"checkm",sep="/")



# Different directories for data

processedRoot <- paste(direRoot,"Analyser","processed", sep="/")

processedRootR <-  paste(processedRoot,"R", sep="/")
processedRootRassembly <- paste(processedRootR,assemblymethod, sep="/")
modelDirectory <-  paste(processedRootR,"model", sep="/")
manuscriptDirectory <-  paste(processedRootR,"manuscript", sep="/")

processedRootExcel <-  paste(processedRoot,"Excel", sep="/")

# working directory is always R root
setwd(workingDirectory)

filesInroot <- list.files(direRoot)
filesInWd <- list.files(workingDirectory)
filesInPprocessedRootR <- list.files(processedRootR)

sampleAsColumns <- function(dataframe)
{
  csvTable <- cbind(rownames(dataframe),dataframe)
  colnames(csvTable)[1] = "sample"
  csvTable
}
