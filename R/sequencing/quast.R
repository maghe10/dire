source(file='common.R')



dir <- quastDirectory
outdir <- paste(processedRootR, "quast", sep="/")






samplesReport <- list.files(dir,"transposed_report.tsv",recursive = TRUE)

sample <- samplesReport[1]

allSamplesReport <- data.frame()

for (sample in samplesReport) {
  aSplit<- strsplit(sample, split = "/")[[1]]  
  aSampleName <- aSplit[[1]]
  reportfilename <- aSplit[[2]]
  print(aSampleName)
  aTransposedTable <- read.csv(file = paste(sep="",dir, "\\" , sample),sep="\t")
#  colnames(aTransposedTable) <- aTransposedTable[1,]
#    aTransposedTable <- aTransposedTable[-1,]
    if(nrow(allSamplesReport) == 0) {
      allSamplesReport <- aTransposedTable
    }
    else {
      allSamplesReport <- rbind(allSamplesReport,aTransposedTable)
    }
}


#################

#samplesReport <- list.files(dir,"report.tsv",recursive = TRUE)

#allSamplesReport <- data.frame()

#for (sample in samplesReport) {
#  aSplit<- strsplit(sample, split = "/")[[1]]  
#  aSampleName <- aSplit[[1]]
#  reportfilename <- aSplit[[2]]
#  if(reportfilename == "report.tsv"){
#      print(aSampleName)
#      aTable <- read.csv(file = paste(sep="",dir, "\\" , sample),sep="\t")
#      aTransposedTable <- data.frame(t(aTable))
#      colnames(aTransposedTable) <- aTransposedTable[1,]
#      aTransposedTable <- aTransposedTable[-1,]
#      if(nrow(allSamplesReport) == 0) {
#        allSamplesReport <- aTransposedTable
#      }
#      else {
#        allSamplesReport <- rbind(allSamplesReport,aTransposedTable)
#      }
#  }
#}

rownames(allSamplesReport) <- allSamplesReport[,1]
allSamplesReport <- allSamplesReport[,-1]

allSamplesReport <- format(allSamplesReport, decimal.mark = ',')

columns <- c("# contigs (>= 0 bp)"	,"# contigs (>= 1000 bp)",	"# contigs (>= 5000 bp)",	"# contigs (>= 10000 bp)",	"# contigs (>= 25000 bp)",	"# contigs (>= 50000 bp)",	"Total length (>= 0 bp)",	"Total length (>= 1000 bp)",	"Total length (>= 5000 bp)",	"Total length (>= 10000 bp)",	"Total length (>= 25000 bp)",	"Total length (>= 50000 bp)",	"# contigs",	"Largest contig",	"Total length"	,"GC (%)"	,"N50",	"N90",	"auN",	"L50",	"L90",	"# N's per 100 kbp")

colnames(allSamplesReport) <- columns

write.csv2(x=allSamplesReport,row.names = TRUE,file= paste(sep="",outdir, "\\" ,"allSamplesReportQuast.csv"))

