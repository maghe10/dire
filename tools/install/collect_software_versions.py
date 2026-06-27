#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import datetime as dt
import glob
import platform
import re
import subprocess
from pathlib import Path
from typing import Optional


DEFAULT_OUTPUT_DIR = "/root/dire/data/Analyser/processed/python"


ENVS = {
    "base": "/root/miniconda3",
    "ariba": "/root/miniconda3/envs/ariba",
    "checkm": "/root/miniconda3/envs/checkm",
    "confindr": "/root/miniconda3/envs/confindr",
    "dire_all": "/root/miniconda3/envs/dire_all",
    "juan": "/root/miniconda3/envs/juan",
    "mlst": "/root/miniconda3/envs/mlst",
}


TOOLS = [
    # base
    {"env": "base", "software": "python", "command": ["python", "--version"]},
    {"env": "base", "software": "conda", "command": ["conda", "--version"]},

    # dire_all
    {"env": "dire_all", "software": "python", "command": ["python", "--version"]},
    {"env": "dire_all", "software": "Trim Galore", "command": ["trim_galore", "--version"]},
    {"env": "dire_all", "software": "Cutadapt", "command": ["cutadapt", "--version"]},
    {"env": "dire_all", "software": "FastQC", "command": ["fastqc", "--version"]},
    {"env": "dire_all", "software": "SPAdes", "command": ["spades.py", "--version"]},
    {"env": "dire_all", "software": "QUAST", "command": ["quast.py", "--version"]},
    {"env": "dire_all", "software": "MultiQC", "command": ["multiqc", "--version"]},
    {"env": "dire_all", "software": "SeqKit", "command": ["seqkit", "version"]},
    {"env": "dire_all", "software": "AMRFinderPlus", "command": ["amrfinder", "--version"]},

    # confindr
    {"env": "confindr", "software": "python", "command": ["python", "--version"]},
    {"env": "confindr", "software": "ConFindr", "command": ["confindr.py", "--version"]},
    {"env": "confindr", "software": "BBDuk", "command": ["bbduk.sh", "--version"]},
    {"env": "confindr", "software": "KMA", "command": "kma 2>&1 | head -n 5"},
    {"env": "confindr", "software": "samtools", "command": "samtools --version | head -n 3"},
    {"env": "confindr", "software": "SeqKit", "command": ["seqkit", "version"]},

    # ariba
    {"env": "ariba", "software": "python", "command": ["python", "--version"]},
    {"env": "ariba", "software": "ARIBA", "command": ["ariba", "version"]},
    {"env": "ariba", "software": "Bowtie2", "command": "bowtie2 --version | head -n 1"},
    {"env": "ariba", "software": "samtools", "command": "samtools --version | head -n 3"},

    # checkm
    {"env": "checkm", "software": "python", "command": ["python", "--version"]},
    {"env": "checkm", "software": "CheckM", "command": "checkm 2>&1 | head -n 20"},
    {"env": "checkm", "software": "pplacer", "command": ["pplacer", "--version"]},
    {"env": "checkm", "software": "Prodigal", "command": ["prodigal", "-v"]},
    {"env": "checkm", "software": "HMMER hmmsearch", "command": "hmmsearch -h | head -n 3"},
    # mlst
    {"env": "mlst", "software": "python", "command": ["python", "--version"]},
    {"env": "mlst", "software": "mlst", "command": ["mlst", "--version"]},

    # juan
    {"env": "juan", "software": "python", "command": ["python", "--version"]},
]


DATABASE_FILES = [
    {
        "software": "AMRFinderPlus database",
        "paths": [
            "/root/miniconda3/envs/dire_all/share/amrfinderplus/data/*/version.txt",
            "/root/miniconda3/envs/dire_all/bin/data/*/version.txt",
        ],
    },
    {
        "software": "ResFinder ARIBA database",
        "paths": [
            "/root/resfinder_ariba_db/RESFINDER_DB_VERSION.txt",
        ],
    },
]


def now_timestamp() -> str:
    return dt.datetime.now().strftime("%Y%m%d_%H%M%S")


def command_to_string(command: list[str] | str) -> str:
    if isinstance(command, str):
        return command
    return " ".join(command)


def run_command(
    env_name: str,
    command: list[str] | str,
    timeout: int = 60,
) -> tuple[str, int]:
    if env_name == "base":
        if isinstance(command, str):
            cmd = ["bash", "-lc", command]
        else:
            cmd = command
    else:
        if isinstance(command, str):
            cmd = ["conda", "run", "-n", env_name, "bash", "-lc", command]
        else:
            cmd = ["conda", "run", "-n", env_name] + command

    try:
        result = subprocess.run(
            cmd,
            text=True,
            capture_output=True,
            timeout=timeout,
            check=False,
        )

        output = "\n".join(
            part for part in [result.stdout.strip(), result.stderr.strip()] if part
        )

        return output.strip(), result.returncode

    except subprocess.TimeoutExpired as exc:
        output_parts = []

        if exc.stdout:
            output_parts.append(str(exc.stdout))

        if exc.stderr:
            output_parts.append(str(exc.stderr))

        return "TIMEOUT: " + "\n".join(output_parts).strip(), 124

    except FileNotFoundError as exc:
        return f"ERROR: {exc}", 127


def first_nonempty_line(text: str) -> str:
    for line in str(text).splitlines():
        line = line.strip()
        if line:
            return line
    return ""


def extract_version(raw_output: str) -> str:
    text = raw_output.strip()

    if not text:
        return ""

    lines = [line.strip() for line in text.splitlines() if line.strip()]
    first_line = lines[0] if lines else ""

    joined = " ".join(lines[:5])

    patterns = [
        r"Python\s+([0-9][A-Za-z0-9._+\-]*)",
        r"conda\s+([0-9][A-Za-z0-9._+\-]*)",
        r"cutadapt\s+([0-9][A-Za-z0-9._+\-]*)",
        r"multiqc[, ]+version\s+([0-9][A-Za-z0-9._+\-]*)",
        r"SPAdes\s+genome\s+assembler\s+v?([0-9][A-Za-z0-9._+\-]*)",
        r"QUAST\s+v?([0-9][A-Za-z0-9._+\-]*)",
        r"FastQC\s+v?([0-9][A-Za-z0-9._+\-]*)",
        r"Trim\s+Galore.*?version\s+([0-9][A-Za-z0-9._+\-]*)",
        r"ARIBA\s+version[: ]+([0-9][A-Za-z0-9._+\-]*)",
        r"bowtie2-align-s\s+version\s+([0-9][A-Za-z0-9._+\-]*)",
        r"samtools\s+([0-9][A-Za-z0-9._+\-]*)",
        r"CheckM\s+v?([0-9][A-Za-z0-9._+\-]*)",
        r"pplacer\s+v?([0-9][A-Za-z0-9._+\-]*)",
        r"Prodigal\s+V?([0-9][A-Za-z0-9._+\-]*)",
        r"HMMER\s+([0-9][A-Za-z0-9._+\-]*)",
        r"mlst\s+([0-9][A-Za-z0-9._+\-]*)",
        r"AMRFinderPlus\s+version\s+([0-9][A-Za-z0-9._+\-]*)",
        r"version\s+([0-9][A-Za-z0-9._+\-]*)",
        r"\bv([0-9][A-Za-z0-9._+\-]*)",
        r"\b([0-9]+(?:\.[0-9A-Za-z_+\-]+)+)\b",
        r"\b([0-9]{4}-[0-9]{2}-[0-9]{2}(?:\.[0-9]+)?)\b",
    ]

    for pattern in patterns:
        match = re.search(pattern, joined, flags=re.IGNORECASE)
        if match:
            return match.group(1)

    return first_line


def find_latest_file(patterns: list[str]) -> Optional[Path]:
    candidates: list[Path] = []

    for pattern in patterns:
        candidates.extend(Path(p) for p in glob.glob(pattern))

    candidates = [p for p in candidates if p.is_file()]

    if not candidates:
        return None

    return max(candidates, key=lambda p: p.stat().st_mtime)


def read_text_file(path: Path) -> tuple[str, str]:
    text = path.read_text(errors="replace").strip()
    version = first_nonempty_line(text)
    return version, text


def collect_tool_versions() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    collected_at = dt.datetime.now().isoformat(timespec="seconds")
    host = platform.node()

    for item in TOOLS:
        env_name = item.get("env", item.get("environment"))
        if env_name is None:
            raise KeyError(f"Missing 'env' key in TOOLS entry: {item}")

        software = item["software"]
        command = item["command"]
        env_path = ENVS.get(env_name, "")


        if env_path and not Path(env_path).exists():
            rows.append(
                {
                    "collected_at": collected_at,
                    "host": host,
                    "environment": env_name,
                    "environment_path": env_path,
                    "software": software,
                    "version": "",
                    "command": command_to_string(command),
                    "return_code": "NA",
                    "status": "environment not found",
                    "raw_output": "",
                }
            )
            continue

        raw_output, return_code = run_command(env_name, command)
        version = extract_version(raw_output)

        rows.append(
            {
                "collected_at": collected_at,
                "host": host,
                "environment": env_name,
                "environment_path": env_path,
                "software": software,
                "version": version,
                "command": command_to_string(command),
                "return_code": str(return_code),
                "status": "ok" if return_code == 0 else "error",
                "raw_output": raw_output,
            }
        )

    return rows


def collect_database_versions() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    collected_at = dt.datetime.now().isoformat(timespec="seconds")
    host = platform.node()

    for db in DATABASE_FILES:
        db_file = find_latest_file(db["paths"])

        if db_file is None:
            rows.append(
                {
                    "collected_at": collected_at,
                    "host": host,
                    "environment": "database",
                    "environment_path": "",
                    "software": db["software"],
                    "version": "",
                    "command": "find " + " ".join(db["paths"]),
                    "return_code": "1",
                    "status": "database version file not found",
                    "raw_output": "",
                }
            )
            continue

        version, raw_output = read_text_file(db_file)

        rows.append(
            {
                "collected_at": collected_at,
                "host": host,
                "environment": "database",
                "environment_path": str(db_file.parent),
                "software": db["software"],
                "version": version,
                "command": f"cat {db_file}",
                "return_code": "0" if version else "1",
                "status": "ok" if version else "error",
                "raw_output": raw_output,
            }
        )

    return rows


def collect_versions() -> list[dict[str, str]]:
    rows = []
    rows.extend(collect_tool_versions())
    rows.extend(collect_database_versions())
    return rows


def write_versions_csv(
    rows: list[dict[str, str]],
    output_dir: str | Path,
) -> Path:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    host = platform.node() or "unknown_host"
    out_file = output_dir / f"software_versions_{host}_{now_timestamp()}.csv"

    fieldnames = [
        "collected_at",
        "host",
        "environment",
        "environment_path",
        "software",
        "version",
        "command",
        "return_code",
        "status",
        "raw_output",
    ]

    with out_file.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter=";",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)

    return out_file


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Collect software and local database versions from DIRE conda "
            "environments and write a semicolon-separated CSV file."
        )
    )

    parser.add_argument(
        "--output-dir",
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory for the CSV file. Default: {DEFAULT_OUTPUT_DIR}",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    rows = collect_versions()
    out_file = write_versions_csv(rows, output_dir=args.output_dir)

    n_ok = sum(row["status"] == "ok" for row in rows)
    n_problem = len(rows) - n_ok

    print(f"Wrote: {out_file}")
    print(f"Rows: {len(rows)}")
    print(f"OK: {n_ok}")
    print(f"Problems: {n_problem}")

    if n_problem:
        print("")
        print("Rows with problems:")
        for row in rows:
            if row["status"] != "ok":
                print(
                    f"- {row['environment']} / {row['software']}: "
                    f"{row['status']} "
                    f"(return_code={row['return_code']})"
                )


if __name__ == "__main__":
    main()