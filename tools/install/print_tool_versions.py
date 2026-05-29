#!/usr/bin/env python3

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import os
import platform
import shutil
import subprocess
from pathlib import Path
from typing import Sequence


DEFAULT_OUTPUT_DIR = "/root/dire/data/Analyser/processed/python"
DEFAULT_ENV_YMLS = [
    "dire_all.yml",
    "confindr.yml",
    "ariba.yml",
    "checkm.yml",
]
DEFAULT_DB_PATHS = {
    "resfinder_ariba_db": "/root/resfinder_ariba_db",
    "confindr_db": "/root/.confindr_db",
}

ENV_TOOLS: dict[str, list[list[str]]] = {
    "dire_all": [
        ["python", "--version"],
        ["trim_galore", "--version"],
        ["cutadapt", "--version"],
        ["fastqc", "--version"],
        ["spades.py", "--version"],
        ["seqkit", "version"],
        ["quast.py", "--version"],
        ["amrfinder", "--version"],
        ["amrfinder_update", "--version"],
        ["multiqc", "--version"],
    ],
    "confindr": [
        ["python", "--version"],
        ["confindr.py", "--version"],
        ["bbduk.sh", "--version"],
        ["bash", "-lc", "kma 2>&1 | head -n 3"],
        ["bash", "-lc", "samtools --version | head -n 3"],
        ["seqkit", "version"],
        ["java", "-version"],
    ],
    "ariba": [
        ["python", "--version"],
        ["ariba", "version"],
        ["bash", "-lc", "bowtie2 --version | head -n 1"],
        ["bash", "-lc", "samtools --version | head -n 3"],
    ],
    "checkm": [
        ["python", "--version"],
        ["bash", "-lc", "checkm 2>&1 | grep -m 1 'CheckM v'"],
        ["pplacer", "--version"],
        ["prodigal", "-v"],
        ["bash", "-lc", "hmmsearch -h | head -n 3"],
    ],
}

CONDA_ENVS = ["dire_all", "confindr", "ariba", "checkm"]


def safe_filename_component(value: str) -> str:
    cleaned = []
    for char in value:
        if char.isalnum() or char in "._-":
            cleaned.append(char)
        else:
            cleaned.append("_")
    result = "".join(cleaned).strip("_")
    return result or "unknown"


def default_output_path(run_tag: str) -> Path:
    return Path(DEFAULT_OUTPUT_DIR) / f"tool_versions_{run_tag}.txt"


def run_command(cmd: Sequence[str], timeout: int = 60) -> tuple[int, str]:
    try:
        result = subprocess.run(
            list(map(str, cmd)),
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=timeout,
        )
        return result.returncode, result.stdout.strip()
    except FileNotFoundError as exc:
        return 127, f"COMMAND NOT FOUND: {exc}"
    except subprocess.TimeoutExpired as exc:
        out = exc.stdout or ""
        return 124, f"TIMEOUT after {timeout} s\n{out}".strip()


def conda_cmd(env: str, inner_cmd: Sequence[str]) -> list[str]:
    return ["conda", "run", "-n", env, "--no-capture-output", *map(str, inner_cmd)]


def write_section(lines: list[str], title: str) -> None:
    lines.append("")
    lines.append("=" * 80)
    lines.append(title)
    lines.append("=" * 80)


def write_subsection(lines: list[str], title: str) -> None:
    lines.append("")
    lines.append("-" * 80)
    lines.append(title)
    lines.append("-" * 80)


def add_command_result(
    lines: list[str],
    cmd: Sequence[str],
    timeout: int,
    label: str | None = None,
) -> None:
    returncode, output = run_command(cmd, timeout=timeout)
    lines.append("")
    if label:
        lines.append(label)
    lines.append(f"$ {' '.join(map(str, cmd))}")
    lines.append(f"exit_code: {returncode}")
    lines.append(output if output else "<no output>")


def add_tool_versions(lines: list[str], timeout: int) -> None:
    write_section(lines, "Tool versions")
    for env, commands in ENV_TOOLS.items():
        write_subsection(lines, f"Environment: {env}")
        for command in commands:
            add_command_result(
                lines,
                conda_cmd(env, command),
                timeout=timeout,
                label=f"[{env}] {' '.join(command)}",
            )


def add_base_info(lines: list[str], timeout: int) -> None:
    write_section(lines, "Base system and conda info")
    for cmd in [
        ["python", "--version"],
        ["conda", "--version"],
        ["conda", "info"],
        ["conda", "env", "list"],
    ]:
        add_command_result(lines, cmd, timeout=timeout)




def add_system_info(lines: list[str], timeout: int) -> None:
    write_section(lines, "Computer and operating system information")

    lines.append("")
    lines.append("Python/platform module:")
    lines.append(f"platform.node: {platform.node()}")
    lines.append(f"platform.platform: {platform.platform()}")
    lines.append(f"platform.system: {platform.system()}")
    lines.append(f"platform.release: {platform.release()}")
    lines.append(f"platform.version: {platform.version()}")
    lines.append(f"platform.machine: {platform.machine()}")
    lines.append(f"platform.processor: {platform.processor()}")
    lines.append(f"python_version: {platform.python_version()}")

    for cmd in [
        ["uname", "-a"],
        ["bash", "-lc", "cat /etc/os-release 2>/dev/null || true"],
        ["hostnamectl"],
        ["nproc"],
        ["lscpu"],
        ["bash", "-lc", "grep -E 'MemTotal|MemFree|MemAvailable|SwapTotal|SwapFree' /proc/meminfo 2>/dev/null || true"],
        ["free", "-h"],
        ["df", "-h"],
        ["lsblk"],
    ]:
        add_command_result(lines, cmd, timeout=timeout)


def add_conda_list(lines: list[str], env: str, timeout: int) -> None:
    write_section(lines, f"Conda package list: {env}")
    add_command_result(lines, ["conda", "list", "-n", env], timeout=timeout)


def add_conda_export(lines: list[str], env: str, timeout: int) -> None:
    write_section(lines, f"Conda environment export: {env}")
    add_command_result(
        lines,
        ["conda", "env", "export", "-n", env, "--no-builds"],
        timeout=timeout,
    )


def sha256_file(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def add_file_snapshot(lines: list[str], path: Path) -> None:
    lines.append("")
    lines.append(f"File: {path}")
    if not path.exists():
        lines.append("status: missing")
        return
    lines.append(f"status: present")
    lines.append(f"size_bytes: {path.stat().st_size}")
    lines.append(f"sha256: {sha256_file(path)}")
    try:
        lines.append("content:")
        lines.append(path.read_text(errors="replace").rstrip())
    except Exception as exc:
        lines.append(f"content: <could not read text: {exc}>")


def add_input_file_snapshots(lines: list[str], paths: list[str]) -> None:
    write_section(lines, "Installation and environment definition files")
    for path_string in paths:
        add_file_snapshot(lines, Path(path_string))


def add_database_info(lines: list[str], timeout: int, db_paths: dict[str, str]) -> None:
    write_section(lines, "Database and reference data information")

    add_command_result(
        lines,
        conda_cmd("dire_all", ["amrfinder_update", "--version"]),
        timeout=timeout,
        label="AMRFinderPlus database/update tool",
    )
    amrfinder_db_version = Path(
        "/root/miniconda3/envs/dire_all/share/amrfinderplus/data/latest/version.txt"
    )

    if amrfinder_db_version.exists():
        lines.append("AMRFinderPlus local database version:")
        lines.append(amrfinder_db_version.read_text(errors="replace").strip())
    else:
        lines.append("AMRFinderPlus local database version: not found")


    for name, path_string in db_paths.items():
        path = Path(path_string)
        write_subsection(lines, f"Database path: {name}")
        lines.append(f"path: {path}")
        lines.append(f"exists: {path.exists()}")
        if not path.exists():
            continue

        if path.is_dir():
            lines.append("top-level files:")
            for child in sorted(path.iterdir())[:50]:
                try:
                    size = child.stat().st_size
                except OSError:
                    size = -1
                kind = "dir" if child.is_dir() else "file"
                lines.append(f"  {kind}\t{size}\t{child.name}")

            for candidate_name in [
                "RESFINDER_DB_VERSION.txt",
                "VERSION",
                "version.txt",
                "README",
                "README.md",
            ]:
                candidate = path / candidate_name
                if candidate.exists() and candidate.is_file():
                    add_file_snapshot(lines, candidate)
        else:
            add_file_snapshot(lines, path)


def is_warning_line(line: str) -> bool:
    """Return True for warning lines that should not be used as summary values."""
    lowered = line.strip().lower()
    return (
        lowered.startswith("warning:")
        or lowered.startswith("[warning]")
        or "python locale settings can't be changed" in lowered
    )


def first_informative_line(text: str) -> str:
    """Return the first useful non-empty line, ignoring warning-only lines.

    The detailed report still contains the complete raw output. This filter only
    affects the short summary, so tools such as QUAST and CheckM are summarized
    by their version line rather than by Python locale warnings.
    """
    first_nonempty = None
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if first_nonempty is None:
            first_nonempty = line
        if is_warning_line(line):
            continue
        return line
    return first_nonempty or "<no output>"


def add_summary_command_result(
    lines: list[str],
    cmd: Sequence[str],
    timeout: int,
    label: str,
) -> None:
    returncode, output = run_command(cmd, timeout=timeout)
    summary = first_informative_line(output)
    lines.append(f"{label}: {summary} (exit_code={returncode})")


def add_summary_tool_versions(lines: list[str], timeout: int) -> None:
    write_section(lines, "Main tool versions")
    for env, commands in ENV_TOOLS.items():
        write_subsection(lines, f"Environment: {env}")
        for command in commands:
            add_summary_command_result(
                lines,
                conda_cmd(env, command),
                timeout=timeout,
                label=" ".join(command),
            )


def add_summary_database_info(lines: list[str], timeout: int, db_paths: dict[str, str]) -> None:
    write_section(lines, "Database and reference data summary")

    add_summary_command_result(
        lines,
        conda_cmd("dire_all", ["amrfinder_update", "--version"]),
        timeout=timeout,
        label="AMRFinderPlus database/update tool",
    )

    amrfinder_db_version = Path(
        "/root/miniconda3/envs/dire_all/share/amrfinderplus/data/latest/version.txt"
    )

    if amrfinder_db_version.exists():
        lines.append("AMRFinderPlus local database version:")
        lines.append(amrfinder_db_version.read_text(errors="replace").strip())
    else:
        lines.append("AMRFinderPlus local database version: not found")
    
    for name, path_string in db_paths.items():
        path = Path(path_string)
        write_subsection(lines, name)
        lines.append(f"path: {path}")
        lines.append(f"exists: {path.exists()}")
        if not path.exists():
            continue
        try:
            stat = path.stat()
            lines.append(f"modified: {dt.datetime.fromtimestamp(stat.st_mtime).isoformat(timespec='seconds')}")
        except OSError as exc:
            lines.append(f"modified: <could not read: {exc}>")

        if path.is_file():
            lines.append(f"size_bytes: {path.stat().st_size}")
            lines.append(f"sha256: {sha256_file(path)}")
            continue

        children = sorted(path.iterdir())
        files = [child for child in children if child.is_file()]
        dirs = [child for child in children if child.is_dir()]
        lines.append(f"top_level_directories: {len(dirs)}")
        lines.append(f"top_level_files: {len(files)}")

        version_files = [
            "RESFINDER_DB_VERSION.txt",
            "VERSION",
            "version.txt",
            "README",
            "README.md",
        ]
        for candidate_name in version_files:
            candidate = path / candidate_name
            if candidate.exists() and candidate.is_file():
                lines.append("")
                lines.append(f"{candidate_name}:")
                try:
                    content = candidate.read_text(errors="replace").strip()
                    lines.append(content if content else "<empty>")
                except Exception as exc:
                    lines.append(f"<could not read: {exc}>")

        lines.append("")
        lines.append("top-level entries, first 20:")
        for child in children[:20]:
            try:
                size = child.stat().st_size
            except OSError:
                size = -1
            kind = "dir" if child.is_dir() else "file"
            lines.append(f"  {kind}\t{size}\t{child.name}")


def default_summary_output_path(run_tag: str) -> Path:
    return Path(DEFAULT_OUTPUT_DIR) / f"tool_versions_summary_{run_tag}.txt"


def write_summary_report(
    summary_output_path: Path,
    created: str,
    run_tag: str,
    timeout: int,
    db_paths: dict[str, str],
) -> None:
    summary_output_path.parent.mkdir(parents=True, exist_ok=True)

    lines: list[str] = []
    lines.append("DIRE pipeline short tool/database summary")
    lines.append(f"Created: {created}")
    lines.append(f"Run tag: {run_tag}")
    lines.append(f"Computer/host: {platform.node()}")
    lines.append(f"Platform: {platform.platform()}")
    lines.append(f"Python running this script: {platform.python_version()}")
    lines.append(f"Working directory: {Path.cwd()}")
    lines.append(f"Conda executable: {shutil.which('conda') or '<not found>'}")

    write_section(lines, "Computer and operating system summary")
    lines.append(f"platform.node: {platform.node()}")
    lines.append(f"platform.platform: {platform.platform()}")
    lines.append(f"platform.system: {platform.system()}")
    lines.append(f"platform.release: {platform.release()}")
    lines.append(f"platform.version: {platform.version()}")
    lines.append(f"platform.machine: {platform.machine()}")
    lines.append(f"python_version: {platform.python_version()}")
    for cmd, label in [
        (["uname", "-a"], "uname"),
        (["bash", "-lc", "cat /etc/os-release 2>/dev/null | head -n 6 || true"], "os-release"),
        (["nproc"], "nproc"),
        (["bash", "-lc", "grep -E 'MemTotal|MemAvailable' /proc/meminfo 2>/dev/null || true"], "memory"),
        (["bash", "-lc", "df -h / /root /mnt/c 2>/dev/null || df -h"], "disk"),
    ]:
        add_command_result(lines, cmd, timeout=timeout, label=label)

    write_section(lines, "Conda summary")
    for cmd in [["conda", "--version"], ["conda", "env", "list"]]:
        add_command_result(lines, cmd, timeout=timeout)

    add_summary_tool_versions(lines, timeout=timeout)
    add_summary_database_info(lines, timeout=timeout, db_paths=db_paths)

    summary_output_path.write_text("\n".join(lines) + "\n")


def parse_db_paths(values: list[str]) -> dict[str, str]:
    db_paths = dict(DEFAULT_DB_PATHS)
    for value in values:
        if "=" not in value:
            raise SystemExit(f"Invalid --db-path value, expected NAME=PATH: {value}")
        name, path = value.split("=", 1)
        db_paths[name] = path
    return db_paths


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Write DIRE pipeline tool versions, computer/OS information, conda environments, install files, database snapshots, and a short summary report.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output report file. Default: tool_versions_<computer>_<timestamp>.txt in /root/dire/data/Analyser/processed/python.",
    )
    parser.add_argument(
        "--summary-output",
        default=None,
        help="Short summary report file. Default: tool_versions_summary_<computer>_<timestamp>.txt in /root/dire/data/Analyser/processed/python.",
    )
    parser.add_argument("--timeout", type=int, default=60)
    parser.add_argument("--skip-conda-list", action="store_true")
    parser.add_argument("--skip-conda-export", action="store_true")
    parser.add_argument("--skip-database-info", action="store_true")
    parser.add_argument(
        "--env-yml",
        action="append",
        default=[],
        help="Environment YAML file to include in the report. Can be repeated.",
    )
    parser.add_argument(
        "--install-script",
        default="install_all.sh",
        help="Installation script to include in the report if present.",
    )
    parser.add_argument(
        "--db-path",
        action="append",
        default=[],
        help="Database path to summarize, format NAME=PATH. Can be repeated.",
    )
    args = parser.parse_args()

    created_dt = dt.datetime.now()
    created = created_dt.isoformat(timespec="seconds")
    computer_name = safe_filename_component(platform.node() or "unknown_host")
    run_tag = f"{computer_name}_{created_dt.strftime('%Y%m%d_%H%M%S')}"

    output_path = Path(args.output) if args.output else default_output_path(run_tag)
    summary_output_path = Path(args.summary_output) if args.summary_output else default_summary_output_path(run_tag)
    db_paths = parse_db_paths(args.db_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    lines: list[str] = []
    lines.append("DIRE pipeline tool versions and system information")
    lines.append(f"Created: {created}")
    lines.append(f"Run tag: {run_tag}")
    lines.append(f"Host: {platform.node()}")
    lines.append(f"Platform: {platform.platform()}")
    lines.append(f"Python running this script: {platform.python_version()}")
    lines.append(f"Working directory: {Path.cwd()}")
    lines.append(f"Conda executable: {shutil.which('conda') or '<not found>'}")
    lines.append(f"PATH: {os.environ.get('PATH', '<not set>')}")

    add_system_info(lines, timeout=args.timeout)
    add_base_info(lines, timeout=args.timeout)
    add_tool_versions(lines, timeout=args.timeout)

    input_files = [args.install_script, *DEFAULT_ENV_YMLS, *args.env_yml]
    add_input_file_snapshots(lines, input_files)

    if not args.skip_database_info:
        add_database_info(lines, timeout=args.timeout, db_paths=db_paths)

    if not args.skip_conda_list:
        for env in CONDA_ENVS:
            add_conda_list(lines, env, timeout=args.timeout)

    if not args.skip_conda_export:
        for env in CONDA_ENVS:
            add_conda_export(lines, env, timeout=args.timeout)

    output_path.write_text("\n".join(lines) + "\n")

    write_summary_report(
        summary_output_path=summary_output_path,
        created=created,
        run_tag=run_tag,
        timeout=args.timeout,
        db_paths=db_paths,
    )

    print(f"Wrote detailed tool versions to: {output_path}")
    print(f"Wrote short tool/database summary to: {summary_output_path}")


if __name__ == "__main__":
    main()
