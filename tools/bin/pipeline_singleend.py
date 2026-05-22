#!/usr/bin/env python3

import argparse
import subprocess
import sys


def run_step(cmd):
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def main():
    parser = argparse.ArgumentParser(
        description="Run full single-end pipeline for one sample"
    )
    parser.add_argument(
        "-s",
        "--sample",
        required=True,
        help="Sample number, e.g. 134",
    )
    args = parser.parse_args()

    sample = str(args.sample)

    steps = [
        ["python", "./trimgalore_singleend.py", "-s", sample],
        ["python", "./uploadtrimmed.py", "-s", sample],
        ["python", "./spades_singleend_standard.py", "-s", sample],
        ["python", "./uploadspades_standard.py", "-s", sample],
        ["python", "./filter_assembly_dire_spades_standard.py", "-s", sample],
        ["python", "./amrfinder_dire_spades_standard.py", "-s", sample],
        ["python", "./kma_amrfinder.py", "-s", sample],
    ]

    for step in steps:
        try:
            run_step(step)
        except subprocess.CalledProcessError as e:
            print(f"Step failed with exit code {e.returncode}: {' '.join(step)}", file=sys.stderr)
            sys.exit(e.returncode)

    print(f"Pipeline completed for sample{sample}")


if __name__ == "__main__":
    main()