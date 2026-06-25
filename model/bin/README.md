# AI susceptibility prediction scripts

This folder contains helper scripts used to run the AI-based antibiotic susceptibility prediction model and to generate conformal prediction outputs for the three susceptibility translation modes used in the manuscript.

The main entry point is:

```bash
python runABC_juan.py
```

Before running the Python scripts, the conformal prediction calibration file must be flattened once using:

```bash
Rscript flatten_calibration.R
```

## Script overview

| Script | Purpose |
|---|---|
| `flatten_calibration.R` | Converts the original nested conformal prediction calibration object (`CICP/calibration.rds`) into a flat CSV file (`CICP/calibration_flat.csv`). This must be run before the Python model scripts. |
| `runABC_juan.py` | Top-level wrapper script. Runs the model for translation modes A, B, and C and for a range of `number_of_choices` values. |
| `runmodel_juan_conformal.py` | Runs the prediction model for one mode/input size and writes standard and conformal prediction outputs. |
| `test_env.py` | Simple environment test for Python, PyTorch, CUDA availability, and GPU functionality. |
| `exec_model.py` | Model execution helper used in the original model workflow. In the current wrapper workflow, `runmodel_juan_conformal.py` loads project-specific model code from the model repository under `program/model/bin`. |

## Expected project structure

The scripts assume the following directory layout inside the DIRE project environment:

```text
/root/dire/
├── data/
│   └── Analyser/
│       └── processed/
│           └── R/
│               └── model/
│                   ├── input/
│                   │   ├── combinations4.csv
│                   │   ├── combinations5.csv
│                   │   ├── ...
│                   │   ├── combinations13.csv
│                   │   ├── modelOutput4.csv
│                   │   ├── modelOutput5.csv
│                   │   ├── ...
│                   │   ├── modelOutput13.csv
│                   │   ├── sirAntibioticsModelWordsJuan2020_12_Mode-A.csv
│                   │   ├── sirAntibioticsModelWordsJuan2020_12_Mode-B.csv
│                   │   └── sirAntibioticsModelWordsJuan2020_12_Mode-C.csv
│                   └── output/
└── program/
    └── model/
        ├── bin/
        │   ├── imports_juan.py
        │   └── exec_model_juan.py
        └── Confidence-based-Prediction-of-Antibiotic-Resistance/
            └── CICP/
                ├── calibration.rds
                └── calibration_flat.csv
```

The paths are currently hard-coded in the scripts:

```text
BIN_PATH       = /root/dire/program/model/bin
DIRE_PATH      = /root/dire/program/model/Confidence-based-Prediction-of-Antibiotic-Resistance
PROCESSED_PATH = /root/dire/data/Analyser/processed/R/model
```

Adjust these paths in the scripts if the project is run from another location.

## Required input files

For each `number_of_choices` value, the Python scripts expect two input files in:

```text
/root/dire/data/Analyser/processed/R/model/input/
```

| File | Description |
|---|---|
| `combinations<N>.csv` | AST input combinations for `N` antibiotics. |
| `modelOutput<N>.csv` | Output template used to store predictions for `N` AST results used as input. |

For each translation mode, the scripts expect one AST input file:

| File | Description |
|---|---|
| `sirAntibioticsModelWordsJuan2020_12_Mode-A.csv` | AST categories translated according to Mode A. |
| `sirAntibioticsModelWordsJuan2020_12_Mode-B.csv` | AST categories translated according to Mode B. |
| `sirAntibioticsModelWordsJuan2020_12_Mode-C.csv` | AST categories translated according to Mode C. |

The mode files are expected to contain a `sample` column and antibiotic susceptibility categories in the format required by the model workflow.

## Step 1: Flatten conformal prediction calibration data

Run this step once before running the Python scripts:

```bash
Rscript flatten_calibration.R
```

This script reads:

```text
CICP/calibration.rds
```

and writes:

```text
CICP/calibration_flat.csv
```

The generated `calibration_flat.csv` is required by `runmodel_juan_conformal.py`.

## Step 2: Test the Python environment

Optional but recommended:

```bash
python test_env.py
```

This prints the Python version, PyTorch version, CUDA availability, number of CUDA devices, and simple CPU/GPU matrix multiplication tests.

## Step 3: Run the model for modes A, B, and C

Run the top-level script:

```bash
python runABC_juan.py
```

By default, this runs:

```text
Modes: A, B, C
number_of_choices: 4 to 13
significant_level: 0.10
base output directory: /root/dire/data/Analyser/processed/R/model/output/temp
```

The wrapper calls `runmodel_juan_conformal.py` repeatedly, once for each mode and each `number_of_choices` value.

## Optional command-line arguments

`runABC_juan.py` supports the following arguments:

```bash
python runABC_juan.py \
  --min_choices 4 \
  --max_choices 13 \
  --modes A,B,C \
  --base_output /root/dire/data/Analyser/processed/R/model/output/temp
```

| Argument | Default | Description |
|---|---:|---|
| `--min_choices` | `4` | Lowest number of AST results used as input. |
| `--max_choices` | `13` | Highest number of AST results used as input. |
| `--modes` | `A,B,C` | Comma-separated list of translation modes to run. |
| `--base_output` | `/root/dire/data/Analyser/processed/R/model/output/temp` | Base output directory. Outputs are written to mode-specific subfolders. |

For example, to run only the six-AST-result analyses for all three modes:

```bash
python runABC_juan.py --min_choices 6 --max_choices 6 --modes A,B,C
```

To run only Mode A for six AST results used as input:

```bash
python runABC_juan.py --min_choices 6 --max_choices 6 --modes A
```

## Output files

For each mode and each `number_of_choices`, output files are written to:

```text
<base_output>/Mode-<MODE>/preds/
```

For example:

```text
/root/dire/data/Analyser/processed/R/model/output/temp/Mode-A/preds/
```

The main output files are:

| File | Description |
|---|---|
| `modelOutput_preds_<N>.csv` | Standard model predictions. |
| `modelOutput_answer_<N>.csv` | Reference AST categories. |
| `modelOutput_confpreds_<N>_01.csv` | Conformal prediction labels at significance level 0.10, corresponding to 90% confidence. |
| `modelOutput_confpreds_<N>_005.csv` | Conformal prediction labels at significance level 0.05, corresponding to 95% confidence. |
| `modelOutput_confpreds_<N>_0025.csv` | Conformal prediction labels at significance level 0.025, corresponding to 97.5% confidence. |

Here, `<N>` is the number of AST results used as input.

## Conformal prediction labels

The conformal prediction output uses the following labels:

| Label | Meaning |
|---|---|
| `S` | Susceptible prediction. |
| `R` | Resistant prediction. |
| `SR` | Ambiguous prediction set containing both susceptible and resistant. |
| empty string | No class passed the conformal prediction threshold. |

## Notes on manuscript terminology

In the manuscript, the preferred terminology is:

```text
AST results for six antibiotics were used as input.
```

rather than:

```text
six input antibiotics
```

This is because the model input consists of interpreted AST results, not the antibiotics themselves.

## Minimal reproducible run

A minimal run for the primary six-AST-result setting in Mode A is:

```bash
Rscript flatten_calibration.R
python test_env.py
python runABC_juan.py --min_choices 6 --max_choices 6 --modes A
```

Expected output directory:

```text
/root/dire/data/Analyser/processed/R/model/output/temp/Mode-A/preds/
```

Expected output files:

```text
modelOutput_preds_6.csv
modelOutput_answer_6.csv
modelOutput_confpreds_6_01.csv
modelOutput_confpreds_6_005.csv
modelOutput_confpreds_6_0025.csv
```

## Troubleshooting

### `CICP/calibration_flat.csv` is missing

Run:

```bash
Rscript flatten_calibration.R
```

and verify that the file was created in:

```text
/root/dire/program/model/Confidence-based-Prediction-of-Antibiotic-Resistance/CICP/
```

### CUDA/GPU is not available

Run:

```bash
python test_env.py
```

If CUDA is not available, check that the correct PyTorch build, GPU drivers, and CUDA runtime are installed.

### Input CSV file not found

Check that all required `combinations<N>.csv`, `modelOutput<N>.csv`, and `sirAntibioticsModelWordsJuan2020_12_Mode-<MODE>.csv` files exist in:

```text
/root/dire/data/Analyser/processed/R/model/input/
```

### Hard-coded paths do not match the local environment

Update the hard-coded path variables in the Python scripts and the working directory logic in `flatten_calibration.R`.

## License and citation

Add repository license information here.

If this code is archived in Zenodo or another repository, cite the archived release DOI in the manuscript data availability statement.
