library(readr)


EUCAST_VERSION  <- "13.0"

NUMBER_OF_ANTIBIOTICS <- 14
ANTIBIOTICS = c("AMP",	"AMC"	,"PIP"	,"TZP",	"CAZ",	"CRO",	"CTX"	,"FEP"	,"CIP"	,"OFX"	,"LVX"	,"MFX"	,"GEN"	,"TOB")


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
#amrfinderDatabase <- "2026-01-21.1"
amrfinderDatabase <- "2026-05-15.1"
amrfinderDirectory <- paste(direRoot,assemblymethod,"amrfinder",amrfinderDatabase,sep="/")

#aribaDirectory <- paste(direRoot,"Illumina","ariba",sep="/")
#aribaResfinderDatabase <- "2026-03-09"
aribaResfinderDatabase <- "2026-05-25"
aribaDirectory <- paste(direRoot,assemblymethod,"ariba",aribaResfinderDatabase,sep="/")

#resfinderDatabase <- "v460"
#resfinderDirectory <- paste(direRoot,assemblymethod,"resfinder",resfinderDatabase,sep="/")

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
processedRootRcommon <-  paste(processedRootR,"common", sep="/")
processedRootRcluster <-  paste(processedRootR,"cluster", sep="/")


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

readModelMillimeterTable <- function()
{
    fileName <- "modelmillimetertable.csv"
    read.csv2(file.path(processedRootRcommon,fileName))
}

readMillimeterTable <- function()
{
  fileName <- "millimetertable.csv"
  read.csv2(file.path(processedRootRcommon,fileName))
}

readDemographicsTable <- function()
{
  fileName <- "demographicstable.csv"
  read.csv2(file.path(processedRootRcommon,fileName))
}

readSirTable <- function(mode = MODE_A)
{
  fileName <- sprintf("sirAntibioticsModel_%s.csv",mode)
  read.csv2(file.path(processedRootRcommon,fileName))
}

filter_and_drop <- function(df, col, values) {

  col_quo <- enquo(col)
  
  # dela upp värden i NA och icke-NA
  values_no_na <- values[!is.na(values)]
  want_na      <- any(is.na(values))
  
  df %>%
    dplyr::filter(
      (!!col_quo %in% values_no_na) | (want_na & is.na(!!col_quo))
    ) %>%
    select(-!!col_quo)
}

run_with_log <- function(
    calls,
    log_file,
    stop_on_error = FALSE,
    append = FALSE
) {
  if (!is.list(calls)) {
    stop("`calls` must be a named list of functions or expressions.")
  }
  
  if (is.null(names(calls)) || any(names(calls) == "")) {
    names(calls) <- paste0("step_", seq_along(calls))
  }
  
  if (!append && file.exists(log_file)) {
    file.remove(log_file)
  }
  
  write_log <- function(...) {
    cat(..., file = log_file, append = TRUE, sep = "")
  }
  
  timestamp <- function() {
    format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  }
  
  write_log(
    "============================================================\n",
    "LOG STARTED: ", timestamp(), "\n",
    "============================================================\n\n"
  )
  
  results <- vector("list", length(calls))
  names(results) <- names(calls)
  
  for (nm in names(calls)) {
    write_log(
      "\n\n",
      "============================================================\n",
      "SECTION: ", nm, "\n",
      "START:   ", timestamp(), "\n",
      "============================================================\n\n"
    )
    
    con <- file(log_file, open = "at")
    sink(con, type = "output", split = FALSE)
    sink(con, type = "message")
    
    result <- tryCatch(
      {
        if (is.function(calls[[nm]])) {
          calls[[nm]]()
        } else {
          eval(calls[[nm]], envir = parent.frame())
        }
      },
      error = function(e) {
        message("ERROR: ", conditionMessage(e))
        structure(list(error = e), class = "logged_error")
      },
      warning = function(w) {
        message("WARNING: ", conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    
    sink(type = "message")
    sink(type = "output")
    close(con)
    
    results[[nm]] <- result
    
    write_log(
      "\n",
      "------------------------------------------------------------\n",
      "END: ", timestamp(), "\n",
      "STATUS: ",
      if (inherits(result, "logged_error")) "ERROR" else "OK",
      "\n",
      "------------------------------------------------------------\n"
    )
    
    if (inherits(result, "logged_error") && stop_on_error) {
      stop("Stopped after error in section: ", nm)
    }
  }
  
  write_log(
    "\n\n",
    "============================================================\n",
    "LOG FINISHED: ", timestamp(), "\n",
    "============================================================\n"
  )
  
  invisible(results)
}


