import os
import sys
import subprocess # just to call an arbitrary command e.g. 'ls'


BIN_PATH = "/root/dire/program/model/bin"
DIRE_PATH = "/root/dire/program/model/Confidence-based-Prediction-of-Antibiotic-Resistance"
PROCESSED_PATH = "/root/dire/data//Analyser/processed/R/model"

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
  importFile = BIN_PATH + "/" + "imports_juan.py"
  exec(open(importFile).read())
  execFile = BIN_PATH + "/" + "exec_model_juan.py"
  exec(open(execFile).read())



parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-choices", "--number_of_choices", default=None, type=int,
                    help='Give a value 1 and 14')
parser.add_argument("-signi", "--significant_level", default=None, nargs='+', type=float,
                    help='Give a value between 0 and 1. If none is given, then conformal prediction is not applied')

parser.add_argument("-outfolder", "--output_folder",default=PROCESSED_PATH)


#SIR_FILE = "input/sirAntibioticsModelWordsNoDate_Mode-C.csv"
parser.add_argument("-sir", "--load_sir_csvfile", default=None)
args = parser.parse_args()


#pPrediction = conformal_pred
#pModelOutput = tmpModelOutput
#pColumn = myCombinations


def asMatrices_fast(df):
    # Work on a minimal copy
    cols_needed = ["id", "Antibiotic", "AST_prediction", "AST_true", "Softmax"]
    d = df[cols_needed].copy()

    # Stack Softmax once -> (n, 2) array of [p_S, p_R]
    P = np.vstack(d["Softmax"].to_numpy()).astype(float)

    # Vectorized threshold → '', 'S', 'R', or 'SR'
    def confpred_from_P(P, tau):
        m = P >= float(tau)          # shape (n, 2)
        out = np.empty(len(P), dtype=object)
        out[:] = ""
        both  = m[:, 0] & m[:, 1]
        onlyS = m[:, 0] & ~m[:, 1]
        onlyR = ~m[:, 0] & m[:, 1]
        out[both]  = "SR"
        out[onlyS] = "S"
        out[onlyR] = "R"
        return out

    d["confpred_01"]   = confpred_from_P(P, 0.90)
    d["confpred_005"]  = confpred_from_P(P, 0.95)
    d["confpred_0025"] = confpred_from_P(P, 0.975)

    # If duplicates exist per (id, Antibiotic), keep first → enables fast pivot
    d = d.drop_duplicates(["id", "Antibiotic"], keep="first")

    # Fast path: pivot (fails if duplicates remain for any reason)
    try:
        wide = d.pivot(index="id",
                       columns="Antibiotic",
                       values=["AST_prediction", "AST_true",
                               "confpred_01", "confpred_005", "confpred_0025"])
    except ValueError:
        # Safe fallback (slower): pivot_table with first
        wide = d.pivot_table(index="id",
                             columns="Antibiotic",
                             values=["AST_prediction", "AST_true",
                                     "confpred_01", "confpred_005", "confpred_0025"],
                             aggfunc="first")  # you used this pattern earlier :contentReference[oaicite:0]{index=0}

    wide = wide.sort_index(axis=1)
    wide.columns.name = None

    preds         = wide["AST_prediction"].reset_index()
    answer        = wide["AST_true"].reset_index()
    confpred_01   = wide["confpred_01"].reset_index()
    confpred_005  = wide["confpred_005"].reset_index()
    confpred_0025 = wide["confpred_0025"].reset_index()

    return preds, answer, confpred_01, confpred_005, confpred_0025


#def asMatrices(df):
#
#    labels = np.array(["S", "R"])
#
#    def conf_set_from_probs(p, tau):
#      p = np.asarray(p)
#      return labels[p >= tau]           # returns array([], dtype='<U1'), ['S'], ['R'], or ['S','R']
#
#    def confpred_string(p, tau):
#      s = conf_set_from_probs(p, tau)
#      return "".join(s) if len(s) > 0 else ""  # None for empty set
#
#    df["confpred_01"] = df["Softmax"].apply(lambda p: confpred_string(p, 0.9))
#    df["confpred_005"] = df["Softmax"].apply(lambda p: confpred_string(p, 0.95))
#    df["confpred_0025"] = df["Softmax"].apply(lambda p: confpred_string(p, 0.975))
#
#    wide = (df.pivot_table(index="id",
#                           columns=["Antibiotic"],
#                           values=["AST_prediction","AST_true","confpred_01","confpred_005","confpred_0025"],
#                           aggfunc="first")
#              .sort_index(axis=1))
#    wide.columns.name = None
#    # Optionally split:
#    preds  = wide["AST_prediction"].reset_index()
#    answer = wide["AST_true"].reset_index()
#    confpred_01 = wide["confpred_01"].reset_index()
#    confpred_005 = wide["confpred_005"].reset_index()
#    confpred_0025 = wide["confpred_0025"].reset_index()
#    #print(wide)
#    #print(softmax)
#    return preds, answer,confpred_01,confpred_005,confpred_0025

def populateModelOutput(pred_df, model_df, col_name):
    cols = pred_df.columns.tolist()
    # build "ANTIBIOTIC_VALUE" or skip NaN/empty
    def joiner(row):
        vals = []
        for a, v in zip(cols, row):
            if isinstance(v, str) and v:           # 'S' or 'R' (or 'S & R')
                v = v.replace(" & ", "")           # SR
                vals.append(f"{a}_{v}")
        return " ".join(vals) if vals else "<empty>"

    model_df[col_name] = pred_df.apply(joiner, axis=1)


if __name__ == '__main__':
  mySignificant_levels = args.significant_level
  myChoices = args.number_of_choices
  with cd(DIRE_PATH):
    data_models,model_antibiotics, model_patient_info, classification_model, device, rank  = INIT()
  with cd(PROCESSED_PATH):
    modelOutput = pd.read_csv('input/modelOutput'+ str(myChoices) + '.csv',sep = ";",dtype=str)
    combinationsInput = pd.read_csv('input/combinations'+ str(myChoices) + '.csv',sep = ";",dtype=str)
    sirAntibioticsWords = pd.read_csv('input/'+args.load_sir_csvfile,sep = ";")
  
    tmpModelOutputPreds = modelOutput.copy()
    tmpModelOutputAnswer = modelOutput.copy()
    tmpModelOutputConfPredsList = list()
    mySignificant_levels = (0.1,0.05,0.025)
    for level in mySignificant_levels:
      tmpModelOutputConfPredsList.append(modelOutput.copy())

    for i in range(0,len(combinationsInput)):
      with cd(DIRE_PATH):
        # Init with the antibiotics for the column index
        # default myCombinations = ('AMC', 'TZP', 'PIP', 'CRO', 'CTX')
        myCombinations = combinationsInput.loc[i]
        readData =  pd.DataFrame.drop(sirAntibioticsWords,columns="sample")
        myData = readData
        df_test_f = RUN(myData,mySignificant_levels,myCombinations,data_models,model_antibiotics, model_patient_info, classification_model, device, rank)
        preds, answer,confpred_01,confpred_005,confpred_0025 = asMatrices_fast(df_test_f)
      comb_word = ''
      for ele in myCombinations:
        comb_word += str(ele)
        comb_word += "_"
      comb_word = comb_word[0:len(comb_word)-1]
      populateModelOutput(preds,tmpModelOutputPreds,comb_word)
      index = 0
      print(mySignificant_levels)
      conformal_preds = (confpred_01,confpred_005,confpred_0025)
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
  with cd(args.output_folder):
    tmpModelOutputPreds.to_csv('modelOutput_preds_'+ str(myChoices) + '.csv',sep = ";",index=False)
    tmpModelOutputAnswer.to_csv('modelOutput_answer_'+ str(myChoices) + '.csv',sep = ";",index=False)  
    index = 0
    for level in mySignificant_levels:
      tmpModelOutputConfPredsList[index].to_csv('modelOutput_confpreds_'+ str(myChoices) + "_" + str(level).replace(".","")  + '.csv',sep = ";",index=False)
      index = index + 1

