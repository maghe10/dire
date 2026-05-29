#!/usr/bin/env python3

from pathlib import Path
import argparse
import shutil
import subprocess
import os


PIPELINE_STEPS = [
    "checkm_lineage",
    "checkm_qa",
    "copy_checkm",
    "copy_checkm_summary",
    "summarize_confindr",
    "summarize_quast",
    "multiqc",
]

DEFAULT_IGNORE_PATTERNS = [
    "*/sample14*",
    "*/sample38*",
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

    parser.add_argument("--trimmed-dir", default="/root/dire/data/Analyser/processed/python/trimmed")
    parser.add_argument("--multiqc-dir", default="/root/dire/data/Analyser/processed/python/multiqc")
    parser.add_argument("--assembly-dir", default="/root/sequencing/intermediate/assembly/spades_standard")
    parser.add_argument("--checkm-dir", default="/root/sequencing/intermediate/checkm")
    parser.add_argument("--storage-checkm-dir", default="/root/dire/data/Analyser/processed/python/checkm")
    parser.add_argument("--quality-dir", default="/root/dire/data/Analyser/processed/python/quality")
    parser.add_argument("--confindr-dir", default="/root/dire/data/Analyser/processed/python/confindr")
    parser.add_argument("--quast-dir", default="/root/dire/data/Analyser/processed/python/quast")
    parser.add_argument("--confindr-summary-output", default="/root/dire/data/Analyser/processed/python/quality/confindr_summary.tsv")
    parser.add_argument("--quast-summary-output", default="/root/dire/data/Analyser/processed/python/quality/quast_summary.tsv")
    parser.add_argument("--summarize-confindr-script", default=str("summarize_confindr.py"))
    parser.add_argument("--summarize-quast-script", default=str("summarize_quast.py"))
    parser.add_argument("--log-dir", default="/root/sequencing/intermediate/logs")

    parser.add_argument("--checkm-threads", type=int, default=4)
    parser.add_argument("--assembly-extension", default="fasta")

    parser.add_argument(
        "--ignore",
        action="append",
        default=DEFAULT_IGNORE_PATTERNS.copy(),
        help="MultiQC ignore pattern. Can be repeated.",
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
    print(f"CheckM threads: {args.checkm_threads}")
    print(f"Quality dir: {args.quality_dir}")

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

    if should_run("multiqc", args.start_at, args.stop_after):
        run_multiqc_all_samples(args)

    print("\nPostpipeline finished successfully.")


if __name__ == "__main__":
    main()
