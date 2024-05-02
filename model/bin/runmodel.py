import os
import sys
import subprocess # just to call an arbitrary command e.g. 'ls'


BIN_PATH = "/root/sequencing/dire/model/bin"
DIRE_PATH = "/root/sequencing/dire/model/dire"
SIR_FILE = "input/sirAntibioticsModelWordsNoDate_Mode-C.csv"
PROCESSED_PATH = "/root/sequencing/dire/Analyser/processed/R/model/tmp"

class cd:
  def __init__(self, newPath):
    self.newPath = newPath
    sys.path.append(newPath)

  def __enter__(self):
    self.savedPath = os.getcwd()
    os.chdir(self.newPath)

  def __exit__(self, etype, value, traceback):
    os.chdir(self.savedPath)
    sys.path.remove(self.newPath)

# Now you can enter the directory like this:
with cd(DIRE_PATH):
  # we are in ~/Library
  importFile = BIN_PATH + "/" + "imports.py"
  exec(open(importFile).read())
  execFile = BIN_PATH + "/" + "exec_model.py"
  exec(open(execFile).read())



parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-choices", "--number_of_choices", default=None, type=int,
                    help='Give a value 1 and 14')
parser.add_argument("-signi", "--significant_level", default=None, nargs='+', type=float,
                    help='Give a value between 0 and 1. If none is given, then conformal prediction is not applied')
args = parser.parse_args()


#pPrediction = conformal_pred
#pModelOutput = tmpModelOutput
#pColumn = myCombinations

def populateModelOutput(pPrediction,pModelOutput,pColumn):
  for i in range(0,len(pPrediction)):
    sentence = ""
    for j in range(0,len(pPrediction.columns)):
      value = pPrediction.iat[i,j]
      if(value=="S" or value=="R" or value=="S & R"):  
        pred = pPrediction.iat[i,j]
        if (pred == "S & R"):
          pred = "SR"
        word = pPrediction.columns[j] + "_" + pred
        if(sentence == ""):
          sentence = word
        else:
          sentence = sentence + " " + word
    if(sentence == ""):
      sentence = "<empty>"
    pModelOutput.at[i,pColumn] = sentence

if __name__ == '__main__':
  mySignificant_levels = args.significant_level
  myChoices = args.number_of_choices
  with cd(PROCESSED_PATH):
    modelOutput = pd.read_csv('input/modelOutput'+ str(myChoices) + '.csv',sep = ";",dtype=str)
    combinationsInput = pd.read_csv('input/combinations'+ str(myChoices) + '.csv',sep = ";",dtype=str)
    sirAntibioticsWords = pd.read_csv(SIR_FILE,sep = ";")
  
    tmpModelOutputPreds = modelOutput.copy()
    tmpModelOutputAnswer = modelOutput.copy()
    tmpModelOutputConfPredsList = list()
    for level in mySignificant_levels:
      tmpModelOutputConfPredsList.append(modelOutput.copy())

    for i in range(0,len(combinationsInput)):
      with cd(DIRE_PATH):
        # Init with the antibiotics for the column index
        # default myCombinations = ('AMC', 'TZP', 'PIP', 'CRO', 'CTX')
        myCombinations = combinationsInput.loc[i]
        myModel,myCombinations = INIT(myCombinations)
        readData =  pd.DataFrame.drop(sirAntibioticsWords,columns="sample")
        myData = readData
#        print(myData)
#        print(mySignificant_levels)
#        print(myCombinations)
        given, preds, answer, conformal_preds = RUN(myData,mySignificant_levels,myModel,myCombinations)
      comb_word = ''
      for ele in myCombinations:
        comb_word += str(ele)
        comb_word += "_"
      comb_word = comb_word[0:len(comb_word)-1]
      populateModelOutput(preds,tmpModelOutputPreds,comb_word)
      index = 0
      for conformal_pred in conformal_preds:
        populateModelOutput(conformal_pred,tmpModelOutputConfPredsList[index],comb_word)
        index = index + 1

      populateModelOutput(answer,tmpModelOutputAnswer,comb_word)
      with cd(PROCESSED_PATH):
        if(i%100==0 or i==len(combinationsInput)-1):
          tmpModelOutputPreds.to_csv('modelOutput_preds_'+ str(myChoices) + '_tmp.csv',sep = ";",index=False)
          tmpModelOutputAnswer.to_csv('modelOutput_answer_'+ str(myChoices) + '_tmp.csv',sep = ";",index=False)
          index = 0
          for level in mySignificant_levels:
            tmpModelOutputConfPredsList[index].to_csv('modelOutput_confpreds_'+ str(myChoices) + "_" + str(level).replace(".","") + '_tmp.csv',sep = ";",index=False)
            index = index + 1
  with cd(PROCESSED_PATH):
    tmpModelOutputPreds.to_csv('modelOutput_preds_'+ str(myChoices) + '.csv',sep = ";",index=False)
    tmpModelOutputAnswer.to_csv('modelOutput_answer_'+ str(myChoices) + '.csv',sep = ";",index=False)  
    index = 0
    for level in mySignificant_levels:
      tmpModelOutputConfPredsList[index].to_csv('modelOutput_confpreds_'+ str(myChoices) + "_" + str(level).replace(".","")  + '.csv',sep = ";",index=False)
      index = index + 1

