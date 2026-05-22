source("common/imports.R")
library(tibble)
library(tidyr)


OUTPUT_DIR <- processedRootRcommon
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

PATIENTS <- "patients"
SAMPLES <- "samples"

ALL_NAMES <- c(
  PATIENTS,
  SAMPLES
)

excelFile =  paste(processedRootExcel, "PatientsAndSamples.xlsx", sep="/")

SHEET_ESTIMATES_PATIENTS <- "EstimatedPatients"
SHEET_RESISTANCE <- "Resistance"




importFromExcel <- function()
{
  estimatedPatients <- sheetToDataFrame(read_xlsx(excelFile,sheet = SHEET_ESTIMATES_PATIENTS))
  resistance <- sheetToDataFrame(read_xlsx(excelFile,sheet = SHEET_RESISTANCE))

  demoIncluded <- read.csv2(file.path(processedRootRcommon,"demographicstable.csv"))
  
  
  demoIncluded2 <- demoIncluded %>%
    mutate(
      Age = case_when(
        age < 18 ~ "<18",
        age >= 18 & age < 65 ~ "18-64",
        age >= 65 ~ ">=65"
      )
    )
  
  included_summary <- demoIncluded2 %>%
    count(Age, sex) %>%
    group_by(Age) %>%
    mutate(prop = n / nrow(demoIncluded2)) %>%
    ungroup()
  
  estimated_long <- estimatedPatients %>%
    pivot_longer(
      cols = c(Male, Female),
      names_to = "sex",
      values_to = "prop"
    ) %>%
    rename(Age = Age)
  
  included_long <- included_summary %>%
    mutate(
      sex = recode(sex, "M" = "Male", "F" = "Female")
    ) %>%
    select(Age, sex, prop)
  
  combined <- bind_rows(
    estimated_long %>% mutate(source = "Estimated"),
    included_long %>% mutate(source = "Included")
  )
  
  combined %>% write.csv2(file.path(processedRootRcommon,paste(PATIENTS,".csv",sep="")),row.names = FALSE)

  resistance %>% select(-`Included(count)`) %>% write.csv2(file.path(processedRootRcommon,paste(SAMPLES,".csv",sep="")),row.names = FALSE)
    
  
  head(demoIncluded)
  estimatedPatients
}



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

