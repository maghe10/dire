source(file = 'model/modelcommon.R')
library(stringr)

suspectedWrongAMCS <- c(13,19,21,22,23,27,34,40,43,54,55,61,63,66,70,75,76,77,83,84,89,91,94,100,103,109,110,111,116,125)
length(suspectedWrongAMCS)

shortFile <- "amrfinderblas.csv"
betalactamaseDataFrame <- read.csv2(file = paste(modelDirectory,shortFile,sep="/"),row.names = "Gene.symbol") 
betalactamaseDataFrame <- t(betalactamaseDataFrame)
rownames <- rownames(betalactamaseDataFrame)
sampleNumbers <- lapply(rownames, function(x) {str_extract(x, "(\\d+)")})
rownames(betalactamaseDataFrame) <- as.integer(sampleNumbers)

betalactamaseDataFrame <- betalactamaseDataFrame[order(as.integer(sampleNumbers)),]

subFrame <- betalactamaseDataFrame[rownames(betalactamaseDataFrame) %in% suspectedWrongAMCS,]
