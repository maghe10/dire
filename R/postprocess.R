# @deprecated => use other genotype scripts
# source(file = 'clear.R')
# source(file = 'genotype/legacy/amrfinder.R')
# ALL()

# @deprecated => use other genotype scripts
# source(file = 'clear.R')
# source(file = 'genotype/legacy/compareGenotypePhenotype.R')
# ALL()

# @deprecated generated tables that was used for figures
# TODO: Create tables to manuscript
# source(file = 'clear.R')
# source(file = 'manuscript/selectExcelStatisticsTables.R')
# ALL()

# @deprecated cluster/cluster.R now does the same
# TODO. More heatmaps. 
#source(file = 'clear.R')
#source(file = 'phenotype/generate_cluster_and_heatmap-kopia.R')
#ALL()

# @deprecated use manuscript/manuscript_figures.R
# source(file = 'clear.R')
# source(file = 'manuscript/manuscript.R')
# ALL()

source('clear.R')
source('genotype/all.R')


source('clear.R')
source('cluster/cluster.R')
ALL()

source('clear.R')
source('manuscript/manuscript.R')

