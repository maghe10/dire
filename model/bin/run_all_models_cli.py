#!/usr/bin/env python3
import os
import argparse
import subprocess
import sys

# -------------------------------------------------------------------
# Global defaults (easy to maintain)
# -------------------------------------------------------------------
DEFAULT_MODEL_VERSION = "251204"
DEFAULT_MIN_CHOICES   = 4
DEFAULT_MAX_CHOICES   = 13
DEFAULT_MODES         = ["Mode-A", "Mode-B", "Mode-C"]

DIRE_PATH = "/root/dire/program/model/Confidence-based-Prediction-of-Antibiotic-Resistance"

class cd:
    def __init__(self, newPath):
        self.newPath = newPath
    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)
    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

# -------------------------------------------------------------------
# Main program
# -------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Run run_model_cpu.py for multiple modes and index ranges."
    )

    parser.add_argument(
        "--modes",
        nargs="+",
        choices=["Mode-A", "Mode-B", "Mode-C"],
        default=DEFAULT_MODES,
        help=f"Lägen som ska köras. Default: {' '.join(DEFAULT_MODES)}"
    )

    parser.add_argument(
        "--min_choices",
        type=int,
        default=DEFAULT_MIN_CHOICES,
        help=f"Minsta number_of_choices (inklusive). Default: {DEFAULT_MIN_CHOICES}"
    )

    parser.add_argument(
        "--max_choices",
        type=int,
        default=DEFAULT_MAX_CHOICES,
        help=f"Största number_of_choices (inklusive). Default: {DEFAULT_MAX_CHOICES}"
    )

    parser.add_argument(
        "--model-version",
        default=DEFAULT_MODEL_VERSION,
        help=f"Modelversion (används i output-path). Default: {DEFAULT_MODEL_VERSION}"
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Skriv bara ut kommandon utan att köra dem."
    )

    args = parser.parse_args()

    if args.min_choices > args.max_choices:
        parser.error("--min_choices får inte vara större än --max_choices")

    modes        = args.modes
    indices      = range(args.min_choices, args.max_choices + 1)
    modelVersion = args.model_version

    # -------------------------------------------------------------------
    # Paths
    # -------------------------------------------------------------------
    base_root   = "/root/dire/data/Analyser/processed/R/model"
    base_input  = f"{base_root}/input"
    base_output = f"{base_root}/output/{modelVersion}"

    model_name  = "model_with_patient_data"
    config_file = "config.yaml"
    mode_type   = "both"

    # -------------------------------------------------------------------
    # Ensure output/mode/predsraw folders exist
    # -------------------------------------------------------------------
    for mode in modes:
        preds_dir = os.path.join(base_output, mode, "predsraw")
        os.makedirs(preds_dir, exist_ok=True)

    # -------------------------------------------------------------------
    # Perform runs
    # -------------------------------------------------------------------
    for mode in modes:
        preds_dir = os.path.join(base_output, mode, "predsraw")

        for idx in indices:
            for sample in range(1, 100):   # 1..99
                input_file  = f"{base_input}/{mode}/sirAntibioticsModelWordsJuan_{idx}_{sample}.csv"
                output_file = f"{preds_dir}/sirAntibioticsModelWordsJuan_{idx}_{sample}.csv"

                cmd = [
                    "python",
                    "run_model_cpu.py",
                    model_name,
                    config_file,
                    mode_type,
                    input_file,
                    output_file
                ]

                print("Running:", " ".join(cmd))

                # Run inside the model directory
                with cd(DIRE_PATH):
                    if not args.dry_run:
                        subprocess.run(cmd, check=True)

# -------------------------------------------------------------------
# Entry point
# -------------------------------------------------------------------
if __name__ == "__main__":
    main()
