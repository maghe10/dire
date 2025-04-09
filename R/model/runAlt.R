source(file = 'clear.R')

MODE <- "Mode-A"
source(file = 'model/modelStatistics.R')
print("starting A")
GENERATE_RISK_STRATIFIED_STATISTICS()



source(file = 'clear.R')
MODE <- "Mode-B"
source(file = 'model/modelStatistics.R')
print("starting B")
GENERATE_RISK_STRATIFIED_STATISTICS()


source(file = 'clear.R')
MODE <- "Mode-C"
source(file = 'model/modelStatistics.R')
print("starting C")
GENERATE_RISK_STRATIFIED_STATISTICS()

source(file = 'clear.R')
source(file = 'model/poststatistics.R')
