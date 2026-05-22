source(file='phenotype/phenotypeCommon.R')
source(file='manuscript/QC-toolkit.R')

outdir <- file.path(processedRootR,"qc")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

QC <- function()
{

  qc <- qc_compare(
    meas_tbl   = QCmillimetertable(),
    target_tbl = readQCTargetTable(),
    range_tbl  = readQCRangeTable(),
    ab_order = antibioticsOrder()
  )
  
  # Tables
  qc_detail  <- qc_table_detail(qc)
  qc_summary <- qc_table_summary(qc)
  
  qc_detail
  qc_summary
  
  # Plots
  qc_plot_bands(qc,filename = file.path(outdir,"qc_bands.png"))

}

ALL <- function()
{
  QC()
}

