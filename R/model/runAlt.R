source(file = 'clear.R')

MODE <- "Mode-A"
source(file = 'model/modelStatistics.R')
print("starting A")
ALL()


source(file = 'clear.R')
MODE <- "Mode-B"
source(file = 'model/modelStatistics.R')
print("starting B")
ALL()


source(file = 'clear.R')
MODE <- "Mode-C"
source(file = 'model/modelStatistics.R')
print("starting C")
ALL()


source(file = 'clear.R')
source(file = 'model/poststatistics.R')
