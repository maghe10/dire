import os
## This directory should be in ~/sequencing/dire/model/bin directory.
## It could be achieved with symbolic links on your unix system
## 
## Nex to the ~/sequencing/dire/model/bin directory there should be a ~/sequencing/dire/model/dire directory with dire source code
## from https://github.com/FraunhoferChalmersCentre/dire
## 
os.chdir("~" +"/sequencing/dire/model/bin")

for i in range(1,14):
  print("Run "+ str(i) + " 0.025")
  os.system("python runmodel.py -choices " + str(i) +" -signi 0.025")

