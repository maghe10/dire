for(letter in c("A","B","C")){
  for(i in 1:13){
  fileConn<-file(paste("run_",i,"_Mode-",letter,".py",sep=""))
  writeLines(c("import os",
               "os.chdir(\"/root/dire/program/model/bin\")",
               "signi = \"0.025 0.05 0.10\"",
                paste("sir_csv = \"sirAntibioticsModelWordsNoDate_Mode-",letter,".csv\"",sep=""),
                paste("output_folder = \"/root/dire/data/Analyser/processed/R/model/tmp/Mode-", letter , "/preds\"",sep=""),
                "os.system(\"mkdir -p \" + output_folder)",
                paste("command = \"python runmodel.py -choices \" + str(",i,") + \" -sir \" + sir_csv + \" -outfolder \" + output_folder + \" --significant_level \" + signi",sep=""),
                "print(\"Running: \" + command)",
                "os.system(command)"
#               paste("print(\"Run \"+ str(",i,") + \" \" + signi)",sep=""),
#               paste("os.system(\"python runmodel.py -choices \" + str(",i,") + \" -sir \" + sir_csv + \" -outfolder \" + output_folder + \" --significant_level \" + signi)")
               ),fileConn)
              
  close(fileConn)
  }
}