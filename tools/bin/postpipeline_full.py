#!/usr/bin/env python3

from pathlib import Path
import argparse
import shutil
import subprocess
import os
import csv


PIPELINE_STEPS = [
    "checkm_lineage",
    "checkm_qa",
    "copy_checkm",
    "copy_checkm_summary",
    "summarize_confindr",
    "summarize_quast",
    "mlst",
    "multiqc",
]

DEFAULT_IGNORE_PATTERNS = [
    "*/sample014*",
    "*/sample038*",
]

DEFAULT_EXCLUDED_MLST_SAMPLES = [
    "sample014",
    "sample038",
]


def should_run(step: str, start_at: str, stop_after: str | None = None) -> bool:
    step_idx = PIPELINE_STEPS.index(step)
    start_idx = PIPELINE_STEPS.index(start_at)

    if step_idx < start_idx:
        return False

    if stop_after is not None:
        stop_idx = PIPELINE_STEPS.index(stop_after)
        if step_idx > stop_idx:
            return False

    return True


def conda_cmd(env: str, cmd: list[str | Path]) -> list[str]:
    return ["conda", "run", "-n", env, "--no-capture-output"] + [str(x) for x in cmd]


def run(cmd: list[str], log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)

    print("\nRunning:")
    print(" ".join(map(str, cmd)))

    with open(log_file, "a") as log:
        log.write("\n\nRunning: " + " ".join(map(str, cmd)) + "\n")
        subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)


def is_nonempty_file(path: Path) -> bool:
    return path.exists() and path.is_file() and path.stat().st_size > 0


def is_nonempty_dir(path: Path) -> bool:
    return path.exists() and path.is_dir() and any(path.iterdir())


def require_dir(path: Path) -> None:
    if not path.exists() or not path.is_dir():
        raise FileNotFoundError(f"Missing directory: {path}")


def require_file(path: Path) -> None:
    if not path.exists() or not path.is_file():
        raise FileNotFoundError(f"Missing file: {path}")


def copy_or_replace_dir(src: Path, dst: Path) -> None:
    require_dir(src)
    dst.parent.mkdir(parents=True, exist_ok=True)

    if dst.exists():
        shutil.rmtree(dst)

    shutil.copytree(src, dst)


def clear_checkm_dir_if_forced(args: argparse.Namespace) -> None:
    checkm_dir = Path(args.checkm_dir)

    if args.force and checkm_dir.exists():
        print(f"Removing existing CheckM directory because --force was used: {checkm_dir}")
        shutil.rmtree(checkm_dir)

    checkm_dir.mkdir(parents=True, exist_ok=True)


def run_checkm_lineage(args: argparse.Namespace) -> None:
    assembly_dir = Path(args.assembly_dir)
    checkm_dir = Path(args.checkm_dir)
    log_dir = Path(args.log_dir)

    require_dir(assembly_dir)

    expected = checkm_dir / "lineage.ms"
    if not args.force and is_nonempty_file(expected):
        print(f"Skipping CheckM lineage_wf: output already exists: {expected}")
        return

    clear_checkm_dir_if_forced(args)

    command = conda_cmd(args.checkm_env, [
        "checkm",
        "lineage_wf",
        "-r",
        "-t", str(args.checkm_threads),
        "-x", args.assembly_extension,
        assembly_dir,
        checkm_dir,
    ])

    run(command, log_dir / "checkm.lineage_wf.log")


def run_checkm_qa(args: argparse.Namespace) -> None:
    """Create both a TSV file for analysis and an -o 1 file for MultiQC."""
    checkm_dir = Path(args.checkm_dir)
    log_dir = Path(args.log_dir)

    lineage_file = checkm_dir / "lineage.ms"
    qa_file = checkm_dir / "qa.tsv"
    multiqc_file = checkm_dir / "checkm_mqc.txt"

    require_file(lineage_file)

    if not args.force and is_nonempty_file(qa_file) and is_nonempty_file(multiqc_file):
        print(f"Skipping CheckM QA: outputs already exist: {qa_file}, {multiqc_file}")
        return

    if args.force or not is_nonempty_file(qa_file):
        command_tsv = conda_cmd(args.checkm_env, [
            "checkm",
            "qa",
            lineage_file,
            checkm_dir,
            "--tab_table",
            "-f", qa_file,
        ])
        run(command_tsv, log_dir / "checkm.qa.tsv.log")

    if args.force or not is_nonempty_file(multiqc_file):
        command_multiqc = conda_cmd(args.checkm_env, [
            "checkm",
            "qa",
            lineage_file,
            checkm_dir,
            "-o", "1",
            "-f", multiqc_file,
        ])
        run(command_multiqc, log_dir / "checkm.qa.multiqc.log")


def copy_checkm_to_storage(args: argparse.Namespace) -> None:
    checkm_dir = Path(args.checkm_dir)
    storage_checkm_dir = Path(args.storage_checkm_dir)

    require_dir(checkm_dir)

    if not args.force and is_nonempty_dir(storage_checkm_dir):
        print(f"Skipping CheckM copy: destination already exists: {storage_checkm_dir}")
        return

    print(f"Copying CheckM directory: {checkm_dir} -> {storage_checkm_dir}")
    copy_or_replace_dir(checkm_dir, storage_checkm_dir)


def copy_checkm_summary_to_quality(args: argparse.Namespace) -> None:
    storage_checkm_dir = Path(args.storage_checkm_dir)
    quality_dir = Path(args.quality_dir)
    source = storage_checkm_dir / "qa.tsv"
    destination = quality_dir / "checkm_summary.tsv"

    require_file(source)
    quality_dir.mkdir(parents=True, exist_ok=True)

    if not args.force and is_nonempty_file(destination):
        print(f"Skipping CheckM summary copy: output already exists: {destination}")
        return

    print(f"Copying CheckM QA summary: {source} -> {destination}")
    shutil.copy2(source, destination)


def run_summary_script(
    args: argparse.Namespace,
    script_path: Path,
    output_path: Path,
    script_args: list[str | Path],
    log_name: str,
    description: str,
) -> None:
    log_dir = Path(args.log_dir)

    require_file(script_path)

    if not args.force and is_nonempty_file(output_path):
        print(f"Skipping {description}: output already exists: {output_path}")
        return

    output_path.parent.mkdir(parents=True, exist_ok=True)

    command = conda_cmd(args.main_env, [
        "python",
        script_path,
        *script_args,
    ])

    run(command, log_dir / log_name)


def run_summarize_confindr(args: argparse.Namespace) -> None:
    tools_dir = Path(args.tools_dir)
    script = tools_dir / args.summarize_confindr_script

    run_summary_script(
        args=args,
        script_path=Path(script),
        output_path=Path(args.confindr_summary_output),
        script_args=[
            "--confindr-dir", Path(args.confindr_dir),
            "--output", Path(args.confindr_summary_output),
        ],
        log_name="summarize.confindr.log",
        description="ConFindr summary",
    )


def run_summarize_quast(args: argparse.Namespace) -> None:
    tools_dir = Path(args.tools_dir)
    script = tools_dir / args.summarize_quast_script

    run_summary_script(
        args=args,
        script_path=Path(script),
        output_path=Path(args.quast_summary_output),
        script_args=[
            "--quast-dir", Path(args.quast_dir),
            "--output", Path(args.quast_summary_output),
        ],
        log_name="summarize.quast.log",
        description="QUAST summary",
    )


def sample_name_from_assembly(path: Path) -> str:
    """Return a stable sample name from common assembly layouts."""
    stem = path.stem

    # SPAdes-style per-sample directories often contain contigs.fasta or scaffolds.fasta.
    if stem in {"contigs", "scaffolds", "assembly"}:
        return path.parent.name

    for suffix in (
        ".assembly",
        "_assembly",
        ".contigs",
        "_contigs",
        ".scaffolds",
        "_scaffolds",
        ".spades",
        "_spades",
    ):
        if stem.endswith(suffix):
            return stem[: -len(suffix)]

    return stem


def normalize_sample_identifier(value: str) -> str:
    """Normalize sample identifiers for exclusion checks."""
    return value.strip().lower().replace("-", "_")


def is_excluded_sample(sample: str, excluded_samples: list[str]) -> bool:
    """Return True when the sample should be excluded from MLST."""
    sample_norm = normalize_sample_identifier(sample)
    excluded_norm = {normalize_sample_identifier(x) for x in excluded_samples}

    if sample_norm in excluded_norm:
        return True

    # Also treat sample14, sample_14, DIRE_EC_014, and DIRE_EC_14 as equivalent.
    digits = "".join(ch for ch in sample_norm if ch.isdigit())
    if digits:
        digits_int = str(int(digits))
        for excluded in excluded_norm:
            excluded_digits = "".join(ch for ch in excluded if ch.isdigit())
            if excluded_digits and str(int(excluded_digits)) == digits_int:
                return True

    return False



def find_assembly_files(
    assembly_dir: Path,
    extension: str,
    excluded_samples: list[str] | None = None,
) -> tuple[list[Path], list[Path]]:
    require_dir(assembly_dir)
    clean_extension = extension.lstrip(".")
    all_files = sorted(assembly_dir.rglob(f"*.{clean_extension}"))

    if not all_files:
        raise FileNotFoundError(
            f"No assembly files with extension .{clean_extension} found under: {assembly_dir}"
        )

    excluded_samples = excluded_samples or []
    included_files: list[Path] = []
    skipped_files: list[Path] = []

    for path in all_files:
        sample = sample_name_from_assembly(path)
        if is_excluded_sample(sample, excluded_samples):
            skipped_files.append(path)
        else:
            included_files.append(path)

    if not included_files:
        raise FileNotFoundError(
            "No assembly files remain after applying MLST sample exclusions."
        )

    return included_files, skipped_files


def parse_mlst_gene_call(call: str) -> tuple[str | None, str]:
    """Parse calls such as adk(53), gyrB(19?), mdh(~9)."""
    if "(" in call and call.endswith(")"):
        gene, value = call.split("(", 1)
        return gene, value[:-1]
    return None, call


def normalize_mlst_output(raw_output: Path, normalized_output: Path) -> None:
    rows: list[dict[str, str]] = []
    gene_order: list[str] = []

    with raw_output.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="	")
        for parts in reader:
            if not parts or len(parts) < 3:
                continue

            assembly_path = Path(parts[0])
            row = {
                "sample": sample_name_from_assembly(assembly_path),
                "file": str(assembly_path),
                "scheme": parts[1],
                "sequence_type": parts[2],
            }

            for call in parts[3:]:
                gene, allele = parse_mlst_gene_call(call)
                if gene is None:
                    gene = f"allele_{len(gene_order) + 1}"
                if gene not in gene_order:
                    gene_order.append(gene)
                row[gene] = allele

            rows.append(row)

    preferred_genes = ["adk", "fumC", "gyrB", "icd", "mdh", "purA", "recA"]
    genes = [g for g in preferred_genes if g in gene_order]
    genes.extend(g for g in gene_order if g not in genes)

    fieldnames = ["sample", "file", "scheme", "sequence_type", *genes]
    normalized_output.parent.mkdir(parents=True, exist_ok=True)

    with normalized_output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="	", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in sorted(rows, key=lambda item: item["sample"]):
            writer.writerow(row)


def run_mlst_all_assemblies(args: argparse.Namespace) -> None:
    assembly_dir = Path(args.assembly_dir)
    mlst_dir = Path(args.mlst_dir)
    output = Path(args.mlst_output)
    raw_output = mlst_dir / "mlst_raw.tsv"
    log_file = Path(args.log_dir) / "mlst.log"

    if not args.force and is_nonempty_file(output):
        print(f"Skipping MLST: output already exists: {output}")
        return

    assembly_files, skipped_files = find_assembly_files(
        assembly_dir,
        args.assembly_extension,
        args.exclude_mlst_sample,
    )
    mlst_dir.mkdir(parents=True, exist_ok=True)
    log_file.parent.mkdir(parents=True, exist_ok=True)

    command = conda_cmd(args.mlst_env, [
        "mlst",
        "--scheme", args.mlst_scheme,
        *assembly_files,
    ])

    print("\nRunning:")
    print(" ".join(map(str, command)))
    print(f"MLST assemblies: {len(assembly_files)}")
    print(f"MLST skipped assemblies: {len(skipped_files)}")
    for skipped in skipped_files:
        print(f"Skipping MLST sample: {sample_name_from_assembly(skipped)} ({skipped})")
    print(f"MLST output: {output}")

    with raw_output.open("w") as out, log_file.open("a") as log:
        log.write("\n\nRunning: " + " ".join(map(str, command)) + "\n")
        log.write(f"MLST assemblies: {len(assembly_files)}\n")
        log.write(f"MLST skipped assemblies: {len(skipped_files)}\n")
        for skipped in skipped_files:
            log.write(f"Skipping MLST sample: {sample_name_from_assembly(skipped)} ({skipped})\n")
        subprocess.run(command, check=True, stdout=out, stderr=log)

    normalize_mlst_output(raw_output, output)


def run_multiqc_all_samples(args: argparse.Namespace) -> None:
    trimmed_dir = Path(args.trimmed_dir)
    storage_checkm_dir = Path(args.storage_checkm_dir)
    multiqc_dir = Path(args.multiqc_dir)
    log_dir = Path(args.log_dir)

    require_dir(trimmed_dir)
    require_dir(storage_checkm_dir)
    multiqc_dir.mkdir(parents=True, exist_ok=True)

    expected = multiqc_dir / "multiqc_report.html"
    if not args.force and is_nonempty_file(expected):
        print(f"Skipping MultiQC: output already exists: {expected}")
        return

    command = conda_cmd(args.main_env, [
        "multiqc",
        "--outdir", multiqc_dir,
        "--force",
    ])

    for pattern in args.ignore:
        command.extend(["--ignore", pattern])

    command.extend([
        str(trimmed_dir),
        str(storage_checkm_dir),
    ])

    run(command, log_dir / "multiqc.log")


def main() -> None:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--force", action="store_true")

    parser.add_argument(
        "--start-at",
        choices=PIPELINE_STEPS,
        default="checkm_lineage",
    )
    parser.add_argument(
        "--stop-after",
        choices=PIPELINE_STEPS,
        default=None,
    )

    parser.add_argument("--main-env", default="dire_all")
    parser.add_argument("--checkm-env", default="checkm")
    parser.add_argument("--mlst-env", default="mlst")

    parser.add_argument("--trimmed-dir", default="/root/dire/data/Analyser/processed/python/trimmed")
    parser.add_argument("--multiqc-dir", default="/root/dire/data/Analyser/processed/python/multiqc")
    parser.add_argument("--assembly-dir", default="/root/sequencing/intermediate/assembly/spades_standard")
    parser.add_argument("--checkm-dir", default="/root/sequencing/intermediate/checkm")
    parser.add_argument("--storage-checkm-dir", default="/root/dire/data/Analyser/processed/python/checkm")
    parser.add_argument("--quality-dir", default="/root/dire/data/Analyser/processed/python/quality")
    parser.add_argument("--confindr-dir", default="/root/dire/data/Analyser/processed/python/confindr")
    parser.add_argument("--quast-dir", default="/root/dire/data/Analyser/processed/python/quast")
    parser.add_argument("--mlst-dir", default="/root/dire/data/Analyser/processed/python/mlst")
    parser.add_argument("--confindr-summary-output", default="/root/dire/data/Analyser/processed/python/quality/confindr_summary.tsv")
    parser.add_argument("--quast-summary-output", default="/root/dire/data/Analyser/processed/python/quality/quast_summary.tsv")
    parser.add_argument("--mlst-output", default="/root/dire/data/Analyser/processed/python/mlst/mlst_ecoli_achtman.tsv")
    parser.add_argument("--summarize-confindr-script", default=str("summarize_confindr.py"))
    parser.add_argument("--summarize-quast-script", default=str("summarize_quast.py"))
    parser.add_argument("--log-dir", default="/root/sequencing/intermediate/logs")

    parser.add_argument("--checkm-threads", type=int, default=4)
    parser.add_argument("--assembly-extension", default="fasta")
    parser.add_argument("--mlst-scheme", default="ecoli_achtman_4")

    parser.add_argument(
        "--ignore",
        action="append",
        default=DEFAULT_IGNORE_PATTERNS.copy(),
        help="MultiQC ignore pattern. Can be repeated.",
    )

    parser.add_argument(
        "--exclude-mlst-sample",
        action="append",
        default=DEFAULT_EXCLUDED_MLST_SAMPLES.copy(),
        help="Sample identifier to exclude from the MLST step. Can be repeated.",
    )
    
    parser.add_argument(
    "--tools-dir",
    default="/root/dire/program/tools/bin",
    help="Directory containing summarize_confindr.py and summarize_quast.py",
)
    args = parser.parse_args()

    if args.stop_after is not None:
        if PIPELINE_STEPS.index(args.stop_after) < PIPELINE_STEPS.index(args.start_at):
            parser.error("--stop-after cannot be before --start-at")

    print("\nPostpipeline steps:")
    print(PIPELINE_STEPS)
    print(f"Start at: {args.start_at}")
    print(f"Stop after: {args.stop_after}")
    print(f"Main env: {args.main_env}")
    print(f"CheckM env: {args.checkm_env}")
    print(f"MLST env: {args.mlst_env}")
    print(f"CheckM threads: {args.checkm_threads}")
    print(f"Quality dir: {args.quality_dir}")
    print(f"MLST dir: {args.mlst_dir}")
    print(f"MLST excluded samples: {args.exclude_mlst_sample}")

    if should_run("checkm_lineage", args.start_at, args.stop_after):
        run_checkm_lineage(args)

    if should_run("checkm_qa", args.start_at, args.stop_after):
        run_checkm_qa(args)

    if should_run("copy_checkm", args.start_at, args.stop_after):
        copy_checkm_to_storage(args)

    if should_run("copy_checkm_summary", args.start_at, args.stop_after):
        copy_checkm_summary_to_quality(args)

    if should_run("summarize_confindr", args.start_at, args.stop_after):
        run_summarize_confindr(args)

    if should_run("summarize_quast", args.start_at, args.stop_after):
        run_summarize_quast(args)

    if should_run("mlst", args.start_at, args.stop_after):
        run_mlst_all_assemblies(args)

    if should_run("multiqc", args.start_at, args.stop_after):
        run_multiqc_all_samples(args)

    print("\nPostpipeline finished successfully.")


if __name__ == "__main__":
    main()
