#pure genotype
source("clear.R")

source("genotype/localBLDB.R")
source("clear.R")

source("genotype/analyze_amrfinder_tsvs.R")
source("clear.R")

source("genotype/analyze_ariba_tsvs.R")
source("clear.R")

source("genotype/resolve_extract_genotype_all_genes.R")
source("clear.R")

source("genotype/sample_genotype_grouping_all_classes.R")
source("clear.R")


#genotype and model
source("genotype/model_genotype_all_classes.R")
RUN()
source("clear.R")

#model => move to
source("genotype/word_input_antibiotic_mmc_analysis.R")
source("clear.R")


#phenotype SIR, mm, and breakpointtable.
source("genotype/extractPhenotype.R")

#phenotype versus genotype
source("clear.R")
source("genotype/predictedPhenotypeFromGenotype_all_classes.R")
source("clear.R")
source("genotype/compareSUR.R")
source("clear.R")

source("genotype/plotSURcompare.R")
source("clear.R")


