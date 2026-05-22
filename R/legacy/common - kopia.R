library(readr)


EUCAST_VERSION  <- "13.0"

# Anders last model
#modelVersion <- "250410"

#Juans first model
#modelVersion <- "251111"

#Juans first model with correct patient data
#modelVersion <- "251119"

#Juans first model with conformalprediction
modelVersion <- "251204"



MODE_A <- "Mode-A"
MODE_B <- "Mode-B"
MODE_C <- "Mode-C"

MODES <- c(MODE_A,MODE_B,MODE_C)


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
#amrfinderDatabase <- "2024-12-18.1"
amrfinderDatabase <- "2026-01-21.1"
amrfinderDirectory <- paste(direRoot,assemblymethod,"amrfinder",amrfinderDatabase,sep="/")

aribaDirectory <- paste(direRoot,"Illumina","ariba",sep="/")

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
manuscriptDirectory <-  paste(processedRootR,"manuscript", modelVersion, sep="/")
manuscriptPlotDirectory <- paste(manuscriptDirectory,"plot", sep="/")

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

normalize_sample_id <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA")] <- NA_character_
  
  out <- x |>
    stringr::str_remove("\\.tsv$") |>
    stringr::str_remove("_amrfinderplus$") |>
    stringr::str_remove("_amrfinder$") |>
    stringr::str_extract("[0-9]+")
  
  out_num <- suppressWarnings(as.integer(out))
  
  out_chr <- rep(NA_character_, length(out_num))
  ok <- !is.na(out_num)
  out_chr[ok] <- sprintf("%03d", out_num[ok])
  
  out_chr
}




