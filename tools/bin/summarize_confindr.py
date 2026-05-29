#!/usr/bin/env python3

from pathlib import Path
import argparse
import csv


def find_confindr_reports(confindr_dir: Path) -> list[Path]:
    patterns = [
        "**/confindr_report.csv",
        "**/confindr_report.tsv",
        "**/*confindr*report*.csv",
        "**/*confindr*report*.tsv",
    ]

    reports = []
    for pattern in patterns:
        reports.extend(confindr_dir.glob(pattern))

    return sorted(set(reports))


def sniff_delimiter(path: Path) -> str:
    if path.suffix.lower() == ".tsv":
        return "\t"
    return ","


def sample_from_path(path: Path) -> str:
    for part in path.parts:
        if part.startswith("sample"):
            return part.replace("sample", "")
    return path.parent.name.replace("sample", "")


def classify_confindr(row: dict[str, str]) -> str:
    text = " ".join(str(x).lower() for x in row.values())

    if "contaminated" in text:
        return "contaminated"

    if "true" in text and "contam" in text:
        return "review"

    if "false" in text and "contam" in text:
        return "pass"

    for col, value in row.items():
        col_l = col.lower()
        if "contam" in col_l and "snv" in col_l:
            try:
                if float(value) > 0:
                    return "review"
            except Exception:
                pass

    return "pass"


def read_report(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    delimiter = sniff_delimiter(path)

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        rows = list(reader)
        fields = list(reader.fieldnames or [])

    for row in rows:
        row["source_file"] = str(path)
        if "sample" not in row or not row["sample"]:
            row["sample"] = sample_from_path(path)
        row["confindr_interpretation"] = classify_confindr(row)

    return rows, fields


def main() -> None:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--confindr-dir",
        default="/root/dire/data/Analyser/processed/python/confindr",
    )
    parser.add_argument(
        "--output",
        default="/root/dire/data/Analyser/processed/python/quality/confindr_summary.tsv",
    )
    args = parser.parse_args()

    confindr_dir = Path(args.confindr_dir)
    output = Path(args.output)

    if not confindr_dir.exists():
        raise FileNotFoundError(f"Missing ConFindr directory: {confindr_dir}")

    reports = find_confindr_reports(confindr_dir)

    if not reports:
        raise FileNotFoundError(f"No ConFindr reports found under: {confindr_dir}")

    all_rows = []
    all_fields = []
    failed = []

    for report in reports:
        try:
            rows, fields = read_report(report)
            all_rows.extend(rows)
            for field in fields:
                if field not in all_fields:
                    all_fields.append(field)
        except Exception as exc:
            failed.append((report, exc))

    if not all_rows:
        raise RuntimeError("No readable ConFindr reports found.")

    first_columns = [
        "sample",
        "confindr_interpretation",
        "genus",
        "num_contaminated_snvs",
        "percent_contam",
        "contam_status",
        "source_file",
    ]

    columns = []
    for col in first_columns + all_fields:
        if col not in columns:
            columns.append(col)

    output.parent.mkdir(parents=True, exist_ok=True)

    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in all_rows:
            writer.writerow(row)

    counts = {}
    for row in all_rows:
        key = row.get("confindr_interpretation", "unknown")
        counts[key] = counts.get(key, 0) + 1

    print(f"Wrote: {output}")
    print()
    print("ConFindr interpretation counts:")
    for key, value in sorted(counts.items()):
        print(f"{key}: {value}")

    if failed:
        print()
        print("Unreadable reports:")
        for report, exc in failed:
            print(f"{report}: {exc}")


if __name__ == "__main__":
    main()