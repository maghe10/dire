source(file='common.R')


outdir <- paste(processedRootR, "confindr", sep="/")

fileTrimmed <- paste(confindrtrimmedDirectory,"confindr_report.csv",sep="/")
tableTrimmed <- read.csv(comment.char = "#",file = fileTrimmed ,sep=",")
tableTrimmed <- tableTrimmed[tableTrimmed$Sample != "sample14",]
write.csv2(x=tableTrimmed,row.names = FALSE,file= paste(sep="",outdir, "\\" ,"confindr_report_trimmed.csv"))

fileRaw <- paste(confindrrawDirectory,"confindr_report.csv",sep="/")
tableRaw <- read.csv(comment.char = "#",file = fileRaw ,sep=",")
tableRaw <- tableRaw[tableRaw$Sample != "sample14",]
write.csv2(x=tableRaw,row.names = FALSE,file= paste(sep="",outdir, "\\" ,"confindr_report_raw.csv"))

