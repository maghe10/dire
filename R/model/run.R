source(file = 'clear.R')

MODE <- "Mode-A"
source(file = 'model/modelStatisticsFoundation.R')
print("starting A modelStatisticsFoundation.R")
ALL()
source(file = 'clear.R')
MODE <- "Mode-A"
source(file = 'model/modelStatistics.R')
print("starting A modelStatistics.R")
ALL()

source(file = 'clear.R')
MODE <- "Mode-B"
source(file = 'model/modelStatisticsFoundation.R')
print("starting B modelStatisticsFoundation.R")
ALL()
source(file = 'clear.R')
MODE <- "Mode-B"
source(file = 'model/modelStatistics.R')
print("starting B modelStatistics.R")
ALL()

source(file = 'clear.R')
MODE <- "Mode-C"
source(file = 'model/modelStatisticsFoundation.R')
print("starting C modelStatisticsFoundation.R")
ALL()
source(file = 'clear.R')
MODE <- "Mode-C"
source(file = 'model/modelStatistics.R')
print("starting C modelStatistics.R")
ALL()

source(file = 'clear.R')
source(file = 'model/poststatistics.R')
