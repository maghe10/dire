source(file='common.R')

allSamples <- function(){
  folder <- amrfinderDirectory
  files <- list.files(folder, pattern = "\\.tsv$", full.names = FALSE)
  ids <- normalize_sample_id(basename(files))
  ids[!ids %in% c("014", "038")]
}


# getGenotypeGroupTable <-function()
# {
#   aFile <- file.path(processedRootRassembly,"genotype","sample_genotype_grouping_simple.csv")
#   stopifnot(file.exists(aFile))
#   read.csv2(aFile)
# }

getGenotypeGroupTable <-function()
{
  aFile <- file.path(processedRootRassembly,"genotype","sample_genotype_grouping_all_classes.csv")
  stopifnot(file.exists(aFile))
  read.csv2(aFile)
}





getGeneMappingTable <-function()
{
  aFile <- file.path(processedRootRassembly,"genotype","resolved_amr_genes_manuscript_table.csv")
  stopifnot(file.exists(aFile))
  read.csv2(aFile)
}
