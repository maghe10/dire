import os
import sys
import subprocess  # just to call an arbitrary command e.g. 'ls'
from pathlib import Path
import argparse

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


# Load project-specific code
with cd(DIRE_PATH):
    importFile = BIN_PATH + "/" + "imports_juan.py"
    exec(open(importFile).read())
    execFile = BIN_PATH + "/" + "exec_model_juan.py"
    exec(open(execFile).read())


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument(
    "-choices", "--number_of_choices", default=None, type=int,
    help="Give a value 1 and 14",
)
parser.add_argument(
    "-signi", "--significant_level", default=None, nargs="+", type=float,
    help="Give a value between 0 and 1. If none is given, then conformal prediction is not applied",
)
parser.add_argument(
    "-outfolder", "--output_folder", default=PROCESSED_PATH,
    help="Output folder for final CSV files",
)
parser.add_argument(
    "-sir", "--load_sir_csvfile", default=None,
    help="Input SIR CSV file name in 'input/'",
)
args = parser.parse_args()


# ---------------------------------------------------------------------
# 1. Load flattened calibration
# ---------------------------------------------------------------------

def load_calibration_flat(path="CICP/calibration_flat.csv") -> pd.DataFrame:
    """
    Load the flattened calibration table produced by flatten_calibration.R.

    Columns: model, pathogen, antibiotic, pheno, score, n, id
    """
    return pd.read_csv(path)


# ---------------------------------------------------------------------
# 2. Fast conformal prediction utilities
# ---------------------------------------------------------------------

def build_calibration_cache(
    calibration_flat: pd.DataFrame,
    model: str = "full",
    pathogen: str = "E. coli",
):
    """
    Precompute calibration arrays for each antibiotic for a given model/pathogen.

    Returns:
        dict: antibiotic -> (cal_sus_vals, cal_res_vals)
        where each is a NumPy array with columns [score, n, id].
    """
    subset = calibration_flat[
        (calibration_flat["model"] == model)
        & (calibration_flat["pathogen"] == pathogen)
    ]

    cache = {}
    for ab, sub in subset.groupby("antibiotic"):
        cal_sus = sub[sub["pheno"] == "susceptible"][["score", "n", "id"]].to_numpy()
        cal_res = sub[sub["pheno"] == "resistant"][["score", "n", "id"]].to_numpy()

        if cal_sus.size == 0 or cal_res.size == 0:
            continue

        # ensure scores are sorted ascending
        idx_sus = np.argsort(cal_sus[:, 0])
        idx_res = np.argsort(cal_res[:, 0])
        cal_sus = cal_sus[idx_sus]
        cal_res = cal_res[idx_res]

        cache[ab] = (cal_sus, cal_res)

    if not cache:
        raise ValueError(f"No calibration entries found for model={model}, pathogen={pathogen}")

    return cache


def fast_p_values_vectorized(cal_vals: np.ndarray, scores_in: np.ndarray) -> np.ndarray:
    """
    Vectorized conformal p-values for many scores at once.

    cal_vals: shape (m, 3) with columns [score, n, id], scores sorted ascending.
    scores_in: shape (k,) array of scores to evaluate.

    Returns:
        p-values of shape (k,).
    """
    scores_cal = cal_vals[:, 0]
    n_vals = cal_vals[:, 1]
    id_vals = cal_vals[:, 2]

    scores_in = np.asarray(scores_in, dtype=float)
    k = scores_in.shape[0]
    p = np.empty(k, dtype=float)

    # index of last cal score <= s
    idx_right = np.searchsorted(scores_cal, scores_in, side="right") - 1

    # Case: no calibration scores <= s
    base_p = 1.0 / (n_vals[0] + 1.0)
    mask_no_le = idx_right < 0
    p[mask_no_le] = base_p

    # Case: at least one <= s
    mask_valid = ~mask_no_le
    if np.any(mask_valid):
        s_valid = scores_in[mask_valid]
        idx_r_valid = idx_right[mask_valid]

        # leftmost index >= s
        idx_left = np.searchsorted(scores_cal, s_valid, side="left")

        # SAFE equality mask – avoid out-of-bounds access
        valid_eq = idx_left < scores_cal.shape[0]
        eq_mask = np.zeros_like(valid_eq, dtype=bool)
        eq_mask[valid_eq] = (scores_cal[idx_left[valid_eq]] == s_valid[valid_eq])

        # Compute final index
        idx = idx_r_valid.copy()
        idx[eq_mask] = np.minimum(idx_r_valid[eq_mask], idx_left[eq_mask])

        n_sel = n_vals[idx]
        id_sel = id_vals[idx]
        p_valid = (id_sel + 1.0) / (n_sel + 1.0)

        p[mask_valid] = p_valid

    return p


def confpred_from_Soft(
    soft: np.ndarray,
    antibiotics: np.ndarray,
    calibration_cache: dict,
    model_name: str = "full",
):
    """
    soft: (n, 2) array with columns [p_S, p_R]
    antibiotics: length-n array of antibiotic codes (e.g. 'CIP')
    calibration_cache: dict antibiotic -> (cal_sus_vals, cal_res_vals)

    Returns:
        dict with keys 'cp_90', 'cp_95', 'cp_975',
        and internal pS, pR arrays (used only inside asMatrices_fast).
    """
    n = soft.shape[0]
    ps_arr = np.empty(n, dtype=float)
    pr_arr = np.empty(n, dtype=float)

    antibiotics = antibiotics.astype(str)
    unique_abs = np.unique(antibiotics)

    # 1) Compute ps/pr in bulk per antibiotic
    for ab in unique_abs:
        if ab not in calibration_cache:
            raise ValueError(f"No calibration in cache for antibiotic '{ab}' (model={model_name})")

        cal_sus_vals, cal_res_vals = calibration_cache[ab]
        idxs = np.where(antibiotics == ab)[0]

        scores_S = soft[idxs, 0]  # P(S)
        scores_R = soft[idxs, 1]  # P(R)

        ps_ab = fast_p_values_vectorized(cal_sus_vals, scores_S)
        pr_ab = fast_p_values_vectorized(cal_res_vals, scores_R)

        ps_arr[idxs] = ps_ab
        pr_arr[idxs] = pr_ab

    # 2) Compute labels vectorized for all rows
    t90 = 0.10
    t95 = 0.05
    t975 = 0.025

    cp_90 = np.full(n, "", dtype=object)
    cp_95 = np.full(n, "", dtype=object)
    cp_975 = np.full(n, "", dtype=object)

    # 90%
    s90 = ps_arr > t90
    r90 = pr_arr > t90
    cp_90[s90 & ~r90] = "S"
    cp_90[~s90 & r90] = "R"
    cp_90[s90 & r90] = "SR"

    # 95%
    s95 = ps_arr > t95
    r95 = pr_arr > t95
    cp_95[s95 & ~r95] = "S"
    cp_95[~s95 & r95] = "R"
    cp_95[s95 & r95] = "SR"

    # 97.5%
    s975 = ps_arr > t975
    r975 = pr_arr > t975
    cp_975[s975 & ~r975] = "S"
    cp_975[~s975 & r975] = "R"
    cp_975[s975 & r975] = "SR"

    return {
        "cp_90": cp_90,
        "cp_95": cp_95,
        "cp_975": cp_975,
    }


# ---------------------------------------------------------------------
# 3. asMatrices_fast — optimized, labels-only
# ---------------------------------------------------------------------

def asMatrices_fast(df: pd.DataFrame, calibration_cache: dict, model_name: str = "full"):
    # Work on a minimal copy
    cols_needed = ["id", "Antibiotic", "AST_prediction", "AST_true", "Output_neural_networks"]
    d = df[cols_needed].copy()

    # Stack logits once -> (n, 2) array
    logits = np.vstack(d["Output_neural_networks"].to_numpy()).astype(float)

    # Vectorized 2-class softmax -> (n, 2) array [p_S, p_R]
    logits_shifted = logits - logits.max(axis=1, keepdims=True)
    exps = np.exp(logits_shifted)
    soft = exps / exps.sum(axis=1, keepdims=True)

    antibiotics = d["Antibiotic"].astype(str).to_numpy()

    # Fast conformal prediction for all rows (returns only cp_* labels now)
    cp_all = confpred_from_Soft(soft, antibiotics, calibration_cache, model_name=model_name)

    # Add conformal prediction labels
    d["confpred_01"] = cp_all["cp_90"]
    d["confpred_005"] = cp_all["cp_95"]
    d["confpred_0025"] = cp_all["cp_975"]

    # If duplicates exist per (id, Antibiotic), keep first → enables fast pivot
    d = d.drop_duplicates(["id", "Antibiotic"], keep="first")

    # Pivot only what we need
    value_cols = [
        "AST_prediction", "AST_true",
        "confpred_01", "confpred_005", "confpred_0025",
    ]

    try:
        wide = d.pivot(
            index="id",
            columns="Antibiotic",
            values=value_cols,
        )
    except ValueError:
        wide = d.pivot_table(
            index="id",
            columns="Antibiotic",
            values=value_cols,
            aggfunc="first",
        )

    wide = wide.sort_index(axis=1)
    wide.columns.name = None

    preds         = wide["AST_prediction"].reset_index()
    answer        = wide["AST_true"].reset_index()
    confpred_01   = wide["confpred_01"].reset_index()
    confpred_005  = wide["confpred_005"].reset_index()
    confpred_0025 = wide["confpred_0025"].reset_index()

    return (
        preds,
        answer,
        confpred_01,
        confpred_005,
        confpred_0025,
    )


# ---------------------------------------------------------------------
# 4. ModelOutput population helper
# ---------------------------------------------------------------------

def populateModelOutput(pred_df: pd.DataFrame, model_df: pd.DataFrame, col_name: str):
    """
    For each row in pred_df, build strings like 'AMC_S CRO_R ...'
    and store them in a single column (col_name) of model_df.
    """
    cols = pred_df.columns.tolist()

    def joiner(row):
        vals = []
        for a, v in zip(cols, row):
            if a == "id":
                continue
            if isinstance(v, str) and v:
                v = v.replace(" & ", "")  # just in case
                vals.append(f"{a}_{v}")
        return " ".join(vals) if vals else "<empty>"

    model_df[col_name] = pred_df.apply(joiner, axis=1)


# ---------------------------------------------------------------------
# 5. Main
# ---------------------------------------------------------------------

if __name__ == "__main__":
    mySignificant_levels = args.significant_level
    myChoices = args.number_of_choices

    with cd(DIRE_PATH):
        calibration_flat = load_calibration_flat("CICP/calibration_flat.csv")
        print("Models:", calibration_flat["model"].unique())
        print("\nPathogens:", calibration_flat["pathogen"].unique())
        print("\nPhenotypes:", calibration_flat["pheno"].unique())
        print("\nAntibiotics")
        print(sorted(calibration_flat["antibiotic"].unique()))

        # Build calibration cache once for model "full" and pathogen "E. coli"
        calibration_cache = build_calibration_cache(calibration_flat, model="full", pathogen="E. coli")

        data_models, model_antibiotics, model_patient_info, classification_model, device, rank = INIT()

    with cd(PROCESSED_PATH):
        modelOutput = pd.read_csv(f"input/modelOutput{myChoices}.csv", sep=";", dtype=str)
        combinationsInput = pd.read_csv(f"input/combinations{myChoices}.csv", sep=";", dtype=str)
        sirAntibioticsWords = pd.read_csv(f"input/{args.load_sir_csvfile}", sep=";")
        MERGE_MODE = False
        if(MERGE_MODE):
            def toInputFile(pData,pCombinations):
                data = pData
                combinations = pCombinations
                combos = combinations.tolist()
                # Split antibiotics for each row
                rows = []
                for _, row in data.iterrows():
                    abx = [a.strip() for a in row["Antibiotic"].split(",")]
                    combo_group = [a for a in abx if any(a.startswith(c) for c in combos)]
                    other_group = [a for a in abx if a not in combo_group]
                    rows.append(["ESCCOL", row["x"], " ".join(combo_group), " ".join(other_group)])

                # Create final DataFrame (no headers)
                input_df = pd.DataFrame(rows) 
                return input_df

            #Create inputframe
            all_inputs = []  # collect all small frames here

            # Create inputframe
            for i in range(len(combinationsInput)):
                # Init with the antibiotics for the column index
                myCombinations = combinationsInput.loc[i]
                readData = pd.DataFrame.drop(sirAntibioticsWords, columns="sample")
                myData = readData

                inputFileData = toInputFile(myData, myCombinations)

                # Optional but very useful: label which combination this frame belongs to
                comb_word = "_".join(str(ele) for ele in myCombinations)
                inputFileData["combination"] = comb_word
                inputFileData["comb_index"] = i

                all_inputs.append(inputFileData)

            # Merge them all into one big frame
            big_input = pd.concat(all_inputs, ignore_index=True)
            print(big_input)
            with cd(DIRE_PATH):
                print("Start merged")
                df_test_f = RUN_INPUT_FRAME(
                        big_input,
                        data_models,
                        model_antibiotics,
                        model_patient_info,
                        classification_model,
                        device,
                        rank,
                    )
                (
                        preds,
                        answer,
                        confpred_01,
                        confpred_005,
                        confpred_0025,
                    ) = asMatrices_fast(
                        df_test_f,
                        calibration_cache=calibration_cache,
                        model_name="full",
                    )
                print("Stop merged")
            with cd(args.output_folder):
                # Core outputs
                preds.to_csv(
                    f"model_preds_{myChoices}.csv", sep=";", index=False
                )
        else :
            tmpModelOutputPreds = modelOutput.copy()
            tmpModelOutputAnswer = modelOutput.copy()
            tmpModelOutputConfPredsList = []

            # override / define significance levels for conformal prediction
            mySignificant_levels = (0.1, 0.05, 0.025)
            for _ in mySignificant_levels:
                tmpModelOutputConfPredsList.append(modelOutput.copy())

            for i in range(len(combinationsInput)):
                with cd(DIRE_PATH):
                    # Init with the antibiotics for the column index
                    myCombinations = combinationsInput.loc[i]
                    readData = pd.DataFrame.drop(sirAntibioticsWords, columns="sample")
                    myData = readData

                    df_test_f = RUN(
                        myData,
                        mySignificant_levels,
                        myCombinations,
                        data_models,
                        model_antibiotics,
                        model_patient_info,
                        classification_model,
                        device,
                        rank,
                    )

                    (
                        preds,
                        answer,
                        confpred_01,
                        confpred_005,
                        confpred_0025,
                    ) = asMatrices_fast(
                        df_test_f,
                        calibration_cache=calibration_cache,
                        model_name="full",
                    )

                comb_word = "_".join(str(ele) for ele in myCombinations)

                # predictions + confpred labels
                populateModelOutput(preds,        tmpModelOutputPreds, comb_word)
                populateModelOutput(answer,       tmpModelOutputAnswer, comb_word)

                conformal_preds = (confpred_01, confpred_005, confpred_0025)
                for idx, conformal_pred in enumerate(conformal_preds):
                    populateModelOutput(conformal_pred, tmpModelOutputConfPredsList[idx], comb_word)

        # Final writes
        with cd(args.output_folder):
            # Core outputs
            tmpModelOutputPreds.to_csv(
                f"modelOutput_preds_{myChoices}.csv", sep=";", index=False
            )
            tmpModelOutputAnswer.to_csv(
                f"modelOutput_answer_{myChoices}.csv", sep=";", index=False
            )

            # Conformal label outputs
            for idx, level in enumerate(mySignificant_levels):
                tmpModelOutputConfPredsList[idx].to_csv(
                    f"modelOutput_confpreds_{myChoices}_{str(level).replace('.', '')}.csv",
                    sep=";",
                    index=False,
                )
