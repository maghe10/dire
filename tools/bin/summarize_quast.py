#!/usr/bin/env python3

from pathlib import Path
import argparse
import csv


KEY_METRICS = [
    "# contigs",
    "Largest contig",
    "Total length",
    "GC (%)",
    "N50",
    "N90",
    "L50",
    "L90",
    "# N's per 100 kbp",
]


def find_quast_reports(quast_dir: Path) -> list[Path]:
    return sorted(quast_dir.glob("**/report.tsv"))


def sample_from_path(path: Path) -> str:
    for part in path.parts:
        if part.startswith("sample"):
            return part.replace("sample", "")
    return path.parent.name.replace("sample", "")


def read_quast_report(path: Path) -> dict[str, str]:
    metrics = {
        "sample": sample_from_path(path),
        "source_file": str(path),
    }

    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) < 2:
                continue
            metric = row[0].strip()
            value = row[1].strip()
            metrics[metric] = value

    return metrics


def main() -> None:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--quast-dir",
        default="/root/dire/data/Analyser/processed/python/quast",
    )
    parser.add_argument(
        "--output",
        default="/root/dire/data/Analyser/processed/python/quality/quast_summary.tsv",
    )
    args = parser.parse_args()

    quast_dir = Path(args.quast_dir)
    output = Path(args.output)

    if not quast_dir.exists():
        raise FileNotFoundError(f"Missing QUAST directory: {quast_dir}")

    reports = find_quast_reports(quast_dir)

    if not reports:
        raise FileNotFoundError(f"No QUAST report.tsv files found under: {quast_dir}")

    rows = [read_quast_report(report) for report in reports]

    all_columns = []
    preferred_columns = ["sample", *KEY_METRICS, "source_file"]

    for col in preferred_columns:
        if any(col in row for row in rows) and col not in all_columns:
            all_columns.append(col)

    for row in rows:
        for col in row:
            if col not in all_columns:
                all_columns.append(col)

    output.parent.mkdir(parents=True, exist_ok=True)

    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=all_columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote: {output}")
    print(f"Samples summarized: {len(rows)}")


if __name__ == "__main__":
    main()