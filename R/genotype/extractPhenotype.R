
source("manuscript/manuscriptCommon.R")
source(file = "model/modelCommon.R")
outdir <- paste(processedRootRassembly, "phenotype", sep = "/")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

readSirTable <- function(mode = "Mode-A")
{                         
  sirTable <- read.csv2(paste(processedRootRcommon,paste("sirAntibioticsModel_", mode ,".csv",sep=""),sep="/"))
  sirTable
} 

readModelMillimeterTable <- function()
{
  fileName <- "modelmillimetertable.csv"
  millimeterTable <- read.csv2(file.path(processedRootRcommon,fileName))
  millimeterTable %>% mutate(sample=normalize_sample_id(sample))
}
  

readBreakpointTable <-function()
{
  fileName <- "breakpoints.csv"
  breakpoints <- read.csv2(file.path(processedRootRcommon,fileName))
  breakpoints %>% 
    filter(antibiotic %in% ANTIBIOTICS)
}

readBreakpointTableOld <- function()
{
  phenotypeExcel =  paste(processedRootExcel, "Tolkad_data.xlsx", sep="/")
  breakpointsSheet <- read_xlsx(phenotypeExcel,sheet = 'Breakpoints invasiv')
  
  breakpointTable <- breakpointsSheet %>% 
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
    filter(antibiotic %in% ANTIBIOTICS)
  breakpointTable
}


combineSir <- function(modeA, modeC, sample_col = "sample") {
  out <- modeA
  cols <- setdiff(names(modeA), sample_col)
  
  out[cols] <- Map(
    function(a, c) ifelse(a == c, a, "U"),
    modeA[cols],
    modeC[cols]
  )
  
  out
}


modeA <- readSirTable(MODE_A)
modeC <- readSirTable(MODE_C)
combined <- combineSir(modeA, modeC)
combined <- combined %>%
  mutate(sample = normalize_sample_id(sample))

millis <- readModelMillimeterTable()


write.csv2(
  combined,
  file.path(outdir, "phenotype_SUR_wide.csv"),
  row.names = FALSE
)

write.csv2(
  millis,
  file.path(outdir, "phenotype_millimeters.csv"),
  row.names = FALSE
)

write.csv2(
  readBreakpointTable(),
  file.path(outdir, "phenotype_breakpoints.csv"),
  na = "",
  row.names = FALSE
)


