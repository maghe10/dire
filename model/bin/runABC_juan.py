#!/usr/bin/env python3
import os
import argparse
import subprocess


def ensure_dir(path: str) -> None:
    """Create folder if missing."""
    os.makedirs(path, exist_ok=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Kör runmodel_juan.py för flera modes och number_of_choices."
    )

    parser.add_argument(
        "--min_choices",
        type=int,
        default=4,
        help="Minsta number_of_choices (inklusive). Default: 4"
    )

    parser.add_argument(
        "--max_choices",
        type=int,
        default=13,
        help="Största number_of_choices (inklusive). Default: 13"
    )

    parser.add_argument(
        "--modes",
        type=str,
        default="A,B,C",
        help="Kommaseparerad lista över modes, t.ex. 'A,B,C'. Default: A,B,C"
    )

    parser.add_argument(
        "--base_output",
        type=str,
        default="/root/dire/data/Analyser/processed/R/model/output/temp",
        help="Basoutputkatalog där Mode-X/preds skapas."
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    # Konvertera 'A,B,C' → ['A','B','C']
    modes = [m.strip() for m in args.modes.split(",") if m.strip()]

    choices_range = range(args.min_choices, args.max_choices + 1)

    # Hårdkodade standardvärden
    csv_template = "sirAntibioticsModelWordsJuan2020_12_Mode-{}.csv"
    run_script = "./runmodel_juan_conformal.py"
    significant_level = "0.10"

    for mode in modes:
        for choice in choices_range:

            print(f"Running model {mode} with number_of_choices = {choice}")

            # Outputdir för mode
            output_folder = os.path.join(args.base_output, f"Mode-{mode}", "preds")
            ensure_dir(output_folder)

            # CSV för just mode
            csv_file = csv_template.format(mode)

            cmd = [
                "python", run_script,
                "--number_of_choices", str(choice),
                "--significant_level", significant_level,
                "--output_folder", output_folder,
                "--load_sir_csvfile", csv_file
            ]

            subprocess.run(cmd, check=True)

            print(f"Completed run for model {mode}, number_of_choices = {choice}")
            print("-" * 40)


if __name__ == "__main__":
    main()
