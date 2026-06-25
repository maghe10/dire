# Preprocessing README

This folder contains the R preprocessing scripts used before running the AI susceptibility prediction model.

The top-level script is:

```r
preprocess.R
```

## Purpose

The preprocessing has two main purposes:

A. **Create model input files in the appropriate format**  
It converts antibiotic susceptibility testing (AST) results and demographic information into the input templates required by the downstream Python model scripts.

B. **Prepare patient and sample information for downstream analyses**  
It imports and formats patient/sample metadata so that the same sample information can be used consistently in later analysis and reporting steps.

## Main scripts

| Script | Purpose |
|---|---|
| `preprocess.R` | Top-level script that runs the preprocessing workflow. |
| `model_input.R` | Creates model input files from AST results and demographic data. |
| `patientsAndSamples.R` | Imports and formats patient and sample information for downstream use. |

## How to run

From the project root, run:

```r
source("preprocess.R")
```

The script depends on shared project setup files and path definitions. Input and output locations are therefore controlled by the common project configuration rather than by this README.
