source(file='common.R')


outdir <- paste(processedRootR, "checkm", sep="/")

qaFile <- paste(checkmDirectory,"qa.tsv",sep="/")
qaTable <- read.csv(file = qaFile,sep="\t")

columns <- c("Bin Id",	"Marker lineage",	"# genomes"	,"# markers"	,"# marker sets"	,"0"	,"1",	"2"	,"3"	,"4"	,"5+"	,"Completeness",	"Contamination",	"Strain heterogeneity")

colnames(qaTable) <- columns
write.csv2(x=qaTable,row.names = FALSE,file= paste(sep="",outdir, "\\" ,"qa.csv"))


