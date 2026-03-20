#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path

# -------------------------------------------------------------------
# Global defaults
# -------------------------------------------------------------------
DEFAULT_MODEL_VERSION = "251204"
DEFAULT_MIN_CHOICES   = 4
DEFAULT_MAX_CHOICES   = 13
DEFAULT_MODES         = ["Mode-A", "Mode-B", "Mode-C"]

RSCRIPT_PATH = r"C:\Users\xhessm\AppData\Local\Programs\R\R-4.5.2\bin\Rscript.exe"

# Folder where THIS Python script lives
SCRIPT_DIR = Path(__file__).resolve().parent

# Project directory used as working dir for Rscript
PROJECT_DIR = (SCRIPT_DIR / ".." / "Confidence-based-Prediction-of-Antibiotic-Resistance").resolve()

# Root where mode/modelVersion folders live
BASE_OUTPUT_ROOT = Path(
    r"C:\Users\xhessm\OneDrive - Västra Götalandsregionen\DIRE\Analyser\processed\R\model\output"
)

# Path to return_prediction_sets.R
R_SCRIPT = PROJECT_DIR / "CICP" / "return_prediction_sets.R"

RUN_MODE = "full"


# -------------------------------------------------------------------
# Run one (mode, idx, sample)
# -------------------------------------------------------------------
def run_single(mode: str, idx: int, sample: int, model_version: str):
    """Run return_prediction_sets.R for a specific mode, index, and sample."""

    mode_dir = BASE_OUTPUT_ROOT / model_version / mode

    # Input: predsraw/
    input_csv = mode_dir / "predsraw" / f"sirAntibioticsModelWordsJuan_{idx}_{sample}.csv"

    # Output: predscp/
    predscp_dir = mode_dir / "predscp"
    predscp_dir.mkdir(parents=True, exist_ok=True)
    output_csv = predscp_dir / f"sirAntibioticsModelWordsJuan_{idx}_{sample}_cp.csv"

    if not input_csv.exists():
        print(f"Skipping missing input file: {input_csv}")
        return

    cmd = [
        RSCRIPT_PATH,
        "--verbose",
        str(R_SCRIPT),
        str(input_csv),
        str(output_csv),
        RUN_MODE,
    ]

    print("Running:", " ".join(cmd))
    result = subprocess.run(cmd, cwd=PROJECT_DIR)
    print(f" → Exit code: {result.returncode}")

    if result.returncode != 0:
        print(f"!! ERROR for {mode}, idx={idx}, sample={sample}\n")


# -------------------------------------------------------------------
# CLI / main
# -------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Create conformal prediction output for multiple modes, indices, and samples."
    )

    parser.add_argument(
        "--modes",
        nargs="+",
        choices=["Mode-A", "Mode-B", "Mode-C"],
        default=DEFAULT_MODES,
        help=f"Default: {' '.join(DEFAULT_MODES)}",
    )

    parser.add_argument(
        "--min_choices",
        type=int,
        default=DEFAULT_MIN_CHOICES,
        help=f"Minimum index (inclusive). Default: {DEFAULT_MIN_CHOICES}",
    )

    parser.add_argument(
        "--max_choices",
        type=int,
        default=DEFAULT_MAX_CHOICES,
        help=f"Maximum index (inclusive). Default: {DEFAULT_MAX_CHOICES}",
    )

    parser.add_argument(
        "--model-version",
        default=DEFAULT_MODEL_VERSION,
        help="Model version folder name. Default hidden.",
    )

    args = parser.parse_args()

    if args.min_choices > args.max_choices:
        parser.error("--min_choices cannot be greater than --max_choices")

    modes   = args.modes
    indices = range(args.min_choices, args.max_choices + 1)
    samples = range(1, 100)     # sample 1 → 99
    model_version = args.model_version

    print(f"PROJECT_DIR   = {PROJECT_DIR}")
    print(f"OUTPUT_ROOT   = {BASE_OUTPUT_ROOT}")
    print(f"Model version = {model_version}")
    print(f"Modes         = {modes}")
    print(f"Indices       = {args.min_choices}–{args.max_choices}")
    print(f"Samples       = 1–99\n")

    for mode in modes:
        print("="*60)
        print(f"Processing mode: {mode}")
        print("="*60)

        for idx in indices:
            for sample in samples:
                run_single(mode, idx, sample, model_version)


if __name__ == "__main__":
    main()
