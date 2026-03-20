import os
import sys
import subprocess # just to call an arbitrary command e.g. 'ls'
from pathlib import Path
import numpy as np
import pandas as pd

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


# ---------------------------------------------------------------------
# 1. Load flattened calibration
# ---------------------------------------------------------------------

def load_calibration_flat(path="CICP/calibration_flat.csv"):
    """
    Load the flattened calibration table produced by flatten_calibration.R.

    Columns: model, pathogen, antibiotic, pheno, score, n, id
    """
    df = pd.read_csv(path)
    return df


# ---------------------------------------------------------------------
# 2. Core math
# ---------------------------------------------------------------------

def _p_from_calibration(cal_df: pd.DataFrame, score: float) -> float:
    """
    Same logic as before, but cal_df is already a pandas DataFrame for
    a single (model, pathogen, antibiotic, pheno) combination.
    """
    scores = cal_df["score"].to_numpy()

    scores_le = np.where(scores <= score)[0]
    scores_eq = np.where(scores == score)[0]

    if len(scores_le) > 0:
        max_le = scores_le.max()
        if len(scores_eq) > 0:
            min_eq = scores_eq.min()
            idx = min(max_le, min_eq)
        else:
            idx = max_le
    else:
        idx = np.inf

    if np.isinf(idx):
        return 1.0 / (cal_df["n"].iloc[0] + 1.0)
    else:
        idx = int(idx)
        return (cal_df["id"].iloc[idx] + 1.0) / (cal_df["n"].iloc[idx] + 1.0)


def _return_prediction_set(
    calibration_flat: pd.DataFrame,
    soft0: float,
    soft1: float,
    model: str,
    antibiotic: str,
    pheno: str,
) -> float:
    """
    Use the flattened calibration table instead of nested lists.

    Assumes pathogen = 'E. coli'.
    """
    pathogen = "E. coli"

    cal_df = calibration_flat[
        (calibration_flat["model"] == model)
        & (calibration_flat["pathogen"] == pathogen)
        & (calibration_flat["antibiotic"] == antibiotic)
        & (calibration_flat["pheno"] == pheno)
    ].copy()

    if cal_df.empty:
        raise ValueError(
            f"No calibration rows for model={model}, pathogen={pathogen}, "
            f"antibiotic={antibiotic}, pheno={pheno}"
        )

    if pheno == "susceptible":
        return _p_from_calibration(cal_df, soft0)
    else:
        return _p_from_calibration(cal_df, soft1)


def prediction_sets_for_confidences(
    calibration_flat: pd.DataFrame,
    soft0: float,
    soft1: float,
    antibiotic: str,
    model: str = "full",
):
    """
    Same API as before; uses flattened calibration.
    """
    ps = _return_prediction_set(
        calibration_flat, soft0, soft1, model, antibiotic, "susceptible"
    )
    pr = _return_prediction_set(
        calibration_flat, soft0, soft1, model, antibiotic, "resistant"
    )
    print(f"[antibiotic={antibiotic}] [ps={ps}] [pr={pr}] [soft0={soft0}] [soft1={soft1}]")

    def label_for_conf(conf: float) -> str:
        threshold = 1.0 - conf
        if ps > threshold and pr > threshold:
            return "SR"
        elif ps > threshold and pr <= threshold:
            return "S"
        elif ps <= threshold and pr > threshold:
            return "R"
        else:
            return ""

    return {
        "cp_90":  label_for_conf(0.90),
        "cp_95":  label_for_conf(0.95),
        "cp_975": label_for_conf(0.975),
    }

def confpred_from_Soft(soft: np.ndarray,
                       antibiotics: np.ndarray,
                       calibration_flat,
                       model_name: str = "full"):
    """
    soft: (n, 2) array with columns [p_S, p_R]
    antibiotics: length-n array of antibiotic codes (e.g. 'CIP')
    calibration_flat: the calibration DataFrame

    Returns:
        dict with keys 'cp_90', 'cp_95', 'cp_975',
        each a length-n array of labels ('S', 'R', 'SR', '').
    """
    n = soft.shape[0]
    cp_90 = np.empty(n, dtype=object)
    cp_95 = np.empty(n, dtype=object)
    cp_975 = np.empty(n, dtype=object)

    for i in range(n):
        soft0 = float(soft[i, 0])
        soft1 = float(soft[i, 1])

        # make sure this is a plain string, not a 1-element array
        ab = antibiotics[i]
        if isinstance(ab, (np.ndarray, list)):
            ab = ab[0]
        ab = str(ab)

        cp = prediction_sets_for_confidences(
            calibration_flat=calibration_flat,
            soft0=soft0,
            soft1=soft1,
            antibiotic=ab,
            model=model_name,
        )

        cp_90[i] = cp["cp_90"]
        cp_95[i] = cp["cp_95"]
        cp_975[i] = cp["cp_975"]

    return {
        "cp_90": cp_90,
        "cp_95": cp_95,
        "cp_975": cp_975,
    }



def softmax_2class(logit_s: float, logit_r: float):
    logits = np.array([logit_s, logit_r], dtype=float)
    exps = np.exp(logits - logits.max())
    probs = exps / exps.sum()
    return probs[0], probs[1]
  
def asMatrices_fast(df):
    # Work on a minimal copy
    cols_needed = ["id", "Antibiotic", "AST_prediction", "AST_true", "Output_neural_networks"]
    d = df[cols_needed].copy()
    
    # Stack soft once -> (n, 2) array of [p_S, p_R]
    logit = np.vstack(d["Output_neural_networks"].to_numpy()).astype(float)
    n = logit.shape[0]
    soft = np.empty((n, 2), dtype=float)
    for i in range(n):
        logit0 = float(logit[i, 0])
        logit1 = float(logit[i, 1])
        soft0,soft1 = softmax_2class(logit0,logit1)
        soft[i, 0] = soft0
        soft[i, 1] = soft1
        
        

    antibiotics = d["Antibiotic"].astype(str).to_numpy()
    #print(antibiotics)

    cp_all = confpred_from_Soft(soft, antibiotics, calibration_flat)
    print(cp_all)

    # Optionally add to DataFrame:
    d["confpred_01"] = cp_all["cp_90"]
    d["confpred_005"] = cp_all["cp_95"]
    d["confpred_0025"] = cp_all["cp_975"]
    
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
    calibration_flat = load_calibration_flat("CICP/calibration_flat.csv")
    print("Models:", calibration_flat["model"].unique())
    print("\nPathogens:", calibration_flat["pathogen"].unique())
    print("\nPhenotypes:", calibration_flat["pheno"].unique())
    print("\nAntibiotics")
    print(sorted(calibration_flat["antibiotic"].unique()))

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

