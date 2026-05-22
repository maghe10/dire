#!/usr/bin/env python3

import argparse
import os
import subprocess
import sys

# --- defaults ---
DB_DIR = "/root/sequencing/kma_db"
AMR_DB = "/root/miniconda3/envs/dire_p31012/share/amrfinderplus/data/2026-03-24.1/AMR_CDS.fa"
TRIM_DIR = "/root/sequencing/storage/trimmed"
OUT_DIR = "/root/sequencing/storage/kma"


def run_cmd(cmd):
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def main():
    parser = argparse.ArgumentParser(description="Run KMA for a sample")
    parser.add_argument("-s", "--sample", required=True, help="Sample number (e.g. 134)")
    args = parser.parse_args()

    sample_num = args.sample
    sample = f"sample{sample_num}"

    reads = os.path.join(TRIM_DIR, f"{sample}_trimmed.fq.gz")
    out_prefix = os.path.join(OUT_DIR, sample)
    db_prefix = os.path.join(DB_DIR, "amrfinder_cds")

    # --- checks ---
    if not os.path.exists(reads):
        print(f"Error: reads not found: {reads}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(DB_DIR, exist_ok=True)
    os.makedirs(OUT_DIR, exist_ok=True)

    # --- build DB if needed ---
    if not os.path.exists(db_prefix + ".length.b"):
        print("Indexing AMR database...")
        run_cmd([
            "kma", "index",
            "-i", AMR_DB,
            "-o", db_prefix
        ])

    # --- run KMA ---
    print(f"Running KMA for {sample}...")
    run_cmd([
        "kma",
        "-i", reads,
        "-o", out_prefix,
        "-t_db", db_prefix,
        "-mem_mode",
        "-ef"
    ])

    # --- filter high-confidence hits ---
    res_file = out_prefix + ".res"
    hc_file = out_prefix + "_hc.res"

    print("Filtering high-confidence hits...")
    with open(res_file) as fin, open(hc_file, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            try:
                identity = float(cols[4])
                coverage = float(cols[5])
                if identity > 90 and coverage > 90:
                    fout.write(line)
            except (IndexError, ValueError):
                continue

    print(f"Done: {sample}")


if __name__ == "__main__":
    main()
