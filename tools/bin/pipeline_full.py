#!/usr/bin/env python3

from pathlib import Path
import argparse
import shutil
import subprocess
import sys


BATCHES = {
    "A": ["002", "003", "004", "005", "006", "008"],
    "B": ["009", "011", "012", "013", "018", "019", "022"],
    "C": ["024", "025", "026", "028", "029", "034", "035", "036"],
    "D": ["039", "040", "043", "044", "045", "047", "048", "050"],
    "E": ["051", "052", "056", "057", "059", "061", "062", "064"],
    "F": ["066", "067", "069", "070", "072", "076", "080", "082"],
    "G": ["083", "084", "086", "087", "088", "089", "091", "093"],
    "H": ["094", "096", "097", "098", "101", "102", "103", "105"],
    "I": ["109", "111", "112", "113", "120"],
    "J": ["017", "020", "023", "027", "032", "037", "041", "046"],
    "K": ["049", "053", "054", "058", "075", "095", "100", "115", "116"],
    "L": ["021", "063", "077", "078", "114"],
    "M": ["055", "060", "065", "079", "085", "110", "121", "122", "123", "124"],
    "N": ["125"],
}

ALL_SAMPLES = [s for batch in BATCHES.values() for s in batch]

PIPELINE_STEPS = [
    "trimgalore",
    "spades",
    "filter",
    "quast",
    "amrfinder",
    "confindr",
    "ariba",
]


def should_run(step, start_at, stop_after=None):
    step_idx = PIPELINE_STEPS.index(step)
    start_idx = PIPELINE_STEPS.index(start_at)

    if step_idx < start_idx:
        return False

    if stop_after is not None:
        stop_idx = PIPELINE_STEPS.index(stop_after)
        if step_idx > stop_idx:
            return False

    return True


def conda_cmd(env, cmd):
    return ["conda", "run", "-n", env, "--no-capture-output"] + [str(x) for x in cmd]


def run(cmd, log_file):
    log_file.parent.mkdir(parents=True, exist_ok=True)
    print("Running:", " ".join(map(str, cmd)))

    with open(log_file, "a") as log:
        log.write("\n\nRunning: " + " ".join(map(str, cmd)) + "\n")
        subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)


def is_nonempty_file(path: Path):
    return path.exists() and path.is_file() and path.stat().st_size > 0


def is_nonempty_dir(path: Path):
    return path.exists() and path.is_dir() and any(path.iterdir())


def require_file(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")


def copy_or_replace(src: Path, dst: Path):
    dst.parent.mkdir(parents=True, exist_ok=True)

    if src.is_dir():
        if dst.exists():
            shutil.rmtree(dst)
        shutil.copytree(src, dst)
    else:
        shutil.copy2(src, dst)


def trimgalore(sample, raw_dir, work_dir, storage, threads, force, env):
    work_out = work_dir / "trimmed" / f"sample{sample}"
    storage_out = storage / "trimmed" / f"sample{sample}"

    expected_r1 = work_out / f"sample{sample}_val_1.fq.gz"
    expected_r2 = work_out / f"sample{sample}_val_2.fq.gz"

    if not force and is_nonempty_file(expected_r1) and is_nonempty_file(expected_r2):
        print(f"Skipping Trim Galore: intermediate output already exists for sample{sample}")
        return work_out

    work_out.mkdir(parents=True, exist_ok=True)

    r1 = raw_dir / f"DIRE_EC_{sample}_R1.fastq.gz"
    r2 = raw_dir / f"DIRE_EC_{sample}_R2.fastq.gz"

    require_file(r1)
    require_file(r2)

    run(
        conda_cmd(env, [
            "trim_galore",
            "--fastqc",
            "--paired",
            "--cores", str(threads),
            "--basename", f"sample{sample}",
            "-o", work_out,
            r1, r2,
        ]),
        work_dir / "logs" / f"sample{sample}.trimgalore.log",
    )

    copy_or_replace(work_out, storage_out)
    return work_out


def spades(sample, trimmed_dir, work_dir, storage, threads, memory, kmers, force, env):
    work_out = work_dir / "spades_standard" / f"sample{sample}"
    storage_out = storage / "spades" / "standard" / f"sample{sample}"

    expected = work_out / "scaffolds.fasta"

    if not force and is_nonempty_file(expected):
        print(f"Skipping SPAdes: intermediate scaffolds already exist for sample{sample}")
        return work_out

    work_out.mkdir(parents=True, exist_ok=True)

    r1 = trimmed_dir / f"sample{sample}_val_1.fq.gz"
    r2 = trimmed_dir / f"sample{sample}_val_2.fq.gz"

    require_file(r1)
    require_file(r2)

    run(
        conda_cmd(env, [
            "spades.py",
            "--isolate",
            "--cov-cutoff", "auto",
            "-k", kmers,
            "-1", r1,
            "-2", r2,
            "-o", work_out,
            "-t", str(threads),
            "-m", str(memory),
        ]),
        work_dir / "logs" / f"sample{sample}.spades.log",
    )

    copy_or_replace(work_out, storage_out)
    return work_out


def filter_assembly(sample, spades_dir, work_dir, storage, min_len, force, env):
    work_outdir = work_dir / "assembly" / "spades_standard"
    storage_outdir = storage / "assembly" / "spades_standard"

    work_outdir.mkdir(parents=True, exist_ok=True)
    storage_outdir.mkdir(parents=True, exist_ok=True)

    scaffolds = spades_dir / "scaffolds.fasta"
    filtered_work = work_outdir / f"sample{sample}.fasta"
    filtered_storage = storage_outdir / f"sample{sample}.fasta"

    if not force and is_nonempty_file(filtered_work):
        print(f"Skipping assembly filtering: intermediate filtered assembly already exists for sample{sample}")
        return filtered_work

    require_file(scaffolds)

    print(f"Filtering assembly: sample{sample}, min length {min_len}")

    with open(filtered_work, "w") as outfile:
        subprocess.run(
            conda_cmd(env, ["seqkit", "seq", "-m", str(min_len), scaffolds]),
            check=True,
            stdout=outfile,
        )

    copy_or_replace(filtered_work, filtered_storage)
    return filtered_work


def quast(sample, assembly, work_dir, storage, threads, force, env):
    work_out = work_dir / "quast" / f"sample{sample}"
    storage_out = storage / "quast" / f"sample{sample}"

    expected = work_out / "report.tsv"

    if not force and is_nonempty_file(expected):
        print(f"Skipping QUAST: intermediate report already exists for sample{sample}")
        return work_out

    work_out.mkdir(parents=True, exist_ok=True)
    require_file(assembly)

    run(
        conda_cmd(env, [
            "quast.py",
            assembly,
            "-o", work_out,
            "-t", str(threads),
        ]),
        work_dir / "logs" / f"sample{sample}.quast.log",
    )

    copy_or_replace(work_out, storage_out)
    return work_out


def amrfinder(sample, assembly, work_dir, storage, organism, force, env):
    work_outdir = work_dir / "amrfinder" / "spades_standard"
    storage_outdir = storage / "amrfinder" / "spades_standard"

    work_outdir.mkdir(parents=True, exist_ok=True)
    storage_outdir.mkdir(parents=True, exist_ok=True)

    work_out = work_outdir / f"sample{sample}.tsv"
    storage_out = storage_outdir / f"sample{sample}.tsv"

    if not force and is_nonempty_file(work_out):
        print(f"Skipping AMRFinder: intermediate output already exists for sample{sample}")
        return work_out

    require_file(assembly)

    print(f"Running AMRFinder: sample{sample}")

    with open(work_out, "w") as outfile:
        subprocess.run(
            conda_cmd(env, [
                "amrfinder",
                "--plus",
                "-n", assembly,
                "-O", organism,
            ]),
            check=True,
            stdout=outfile,
        )

    copy_or_replace(work_out, storage_out)
    return work_out


def confindr(sample, trimmed_dir, work_dir, storage, db, threads, base_cutoff, force, env):
    work_out = work_dir / "confindr" / f"sample{sample}"
    storage_out = storage / "confindr" / f"sample{sample}"

    if not force and is_nonempty_dir(work_out):
        print(f"Skipping ConFindr: intermediate output already exists for sample{sample}")
        return work_out

    if force and work_out.exists():
        shutil.rmtree(work_out)

    work_out.mkdir(parents=True, exist_ok=True)

    r1 = trimmed_dir / f"sample{sample}_val_1.fq.gz"
    r2 = trimmed_dir / f"sample{sample}_val_2.fq.gz"

    require_file(r1)
    require_file(r2)

    code = (
        "from confindr_src import confindr; "
        "import sys; "
        "pair=[sys.argv[1], sys.argv[2]]; "
        "confindr.find_contamination("
        "pair=pair, "
        "forward_id='_1', "
        "threads=int(sys.argv[3]), "
        "base_cutoff=int(sys.argv[4]), "
        "output_folder=sys.argv[5], "
        "databases_folder=sys.argv[6]"
        ")"
    )

    run(
        conda_cmd(env, [
            "python",
            "-c", code,
            r1,
            r2,
            str(threads),
            str(base_cutoff),
            work_out,
            db,
        ]),
        work_dir / "logs" / f"sample{sample}.confindr.log",
    )

    copy_or_replace(work_out, storage_out)
    return work_out


def ariba(sample, trimmed_dir, work_dir, storage, db, force, env):
    work_out = work_dir / "ariba" / f"sample{sample}"
    storage_out = storage / "ariba" / f"sample{sample}"
    expected = work_out / "report.tsv"

    if not force and is_nonempty_file(expected):
        print(f"Skipping ARIBA: intermediate report already exists for sample{sample}")
        return work_out

    if work_out.exists():
        if force or not is_nonempty_file(expected):
            print(f"Removing incomplete ARIBA output: {work_out}")
            shutil.rmtree(work_out)

    work_out.parent.mkdir(parents=True, exist_ok=True)

    r1 = trimmed_dir / f"sample{sample}_val_1.fq.gz"
    r2 = trimmed_dir / f"sample{sample}_val_2.fq.gz"

    require_file(r1)
    require_file(r2)

    run(
        conda_cmd(env, [
            "ariba",
            "run",
            db,
            r1,
            r2,
            work_out,
        ]),
        work_dir / "logs" / f"sample{sample}.ariba.log",
    )

    copy_or_replace(work_out, storage_out)
    return work_out


def run_sample(sample, args):
    raw_dir = Path(args.raw_dir)
    work_dir = Path(args.work_dir)
    storage = Path(args.storage_dir)

    print(f"\n=== Processing sample{sample} ===")
    print(f"Start at: {args.start_at}")
    print(f"Stop after: {args.stop_after}")

    trimmed_work = work_dir / "trimmed" / f"sample{sample}"
    spades_work = work_dir / "spades_standard" / f"sample{sample}"
    assembly_work = work_dir / "assembly" / "spades_standard" / f"sample{sample}.fasta"

    if should_run("trimgalore", args.start_at, args.stop_after):
        trimmed_work = trimgalore(
            sample=sample,
            raw_dir=raw_dir,
            work_dir=work_dir,
            storage=storage,
            threads=args.threads,
            force=args.force,
            env=args.main_env,
        )

    if should_run("spades", args.start_at, args.stop_after):
        require_file(trimmed_work / f"sample{sample}_val_1.fq.gz")
        require_file(trimmed_work / f"sample{sample}_val_2.fq.gz")

        spades_work = spades(
            sample=sample,
            trimmed_dir=trimmed_work,
            work_dir=work_dir,
            storage=storage,
            threads=args.threads,
            memory=args.memory,
            kmers=args.kmers,
            force=args.force,
            env=args.main_env,
        )

    if should_run("filter", args.start_at, args.stop_after):
        require_file(spades_work / "scaffolds.fasta")

        assembly_work = filter_assembly(
            sample=sample,
            spades_dir=spades_work,
            work_dir=work_dir,
            storage=storage,
            min_len=args.min_contig_len,
            force=args.force,
            env=args.main_env,
        )

    if should_run("quast", args.start_at, args.stop_after):
        require_file(assembly_work)

        quast(
            sample=sample,
            assembly=assembly_work,
            work_dir=work_dir,
            storage=storage,
            threads=args.threads,
            force=args.force,
            env=args.main_env,
        )

    if should_run("amrfinder", args.start_at, args.stop_after):
        require_file(assembly_work)

        amrfinder(
            sample=sample,
            assembly=assembly_work,
            work_dir=work_dir,
            storage=storage,
            organism=args.organism,
            force=args.force,
            env=args.main_env,
        )

    if should_run("confindr", args.start_at, args.stop_after) and not args.skip_confindr:
        require_file(trimmed_work / f"sample{sample}_val_1.fq.gz")
        require_file(trimmed_work / f"sample{sample}_val_2.fq.gz")

        confindr(
            sample=sample,
            trimmed_dir=trimmed_work,
            work_dir=work_dir,
            storage=storage,
            db=args.confindr_db,
            threads=args.confindr_threads,
            base_cutoff=args.confindr_base_cutoff,
            force=args.force,
            env=args.confindr_env,
        )

    if should_run("ariba", args.start_at, args.stop_after) and not args.skip_ariba:
        require_file(trimmed_work / f"sample{sample}_val_1.fq.gz")
        require_file(trimmed_work / f"sample{sample}_val_2.fq.gz")

        ariba(
            sample=sample,
            trimmed_dir=trimmed_work,
            work_dir=work_dir,
            storage=storage,
            db=args.ariba_db,
            force=args.force,
            env=args.ariba_env,
        )

    print(f"Finished sample{sample}")


def select_samples(args, parser):
    selected = sum([
        args.sample is not None,
        args.batch is not None,
        args.all,
    ])

    if selected == 0:
        parser.error("No samples selected. Use one of: -s SAMPLE, -b BATCH, or --all")

    if selected > 1:
        parser.error("Use only one of: --sample, --batch, or --all")

    if args.sample:
        return [args.sample]

    if args.batch:
        batch = args.batch.upper()
        if batch not in BATCHES:
            parser.error(
                f"Unknown batch '{batch}'. Available batches: {', '.join(BATCHES.keys())}"
            )
        return BATCHES[batch]

    confirm = input(f"Run ALL {len(ALL_SAMPLES)} samples? [y/N]: ")
    if confirm.lower() != "y":
        print("Aborted.")
        raise SystemExit

    return ALL_SAMPLES


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("-s", "--sample", default=None)
    parser.add_argument("-b", "--batch", default=None)
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--force", action="store_true")

    parser.add_argument(
        "--start-at",
        choices=PIPELINE_STEPS,
        default="trimgalore",
    )

    parser.add_argument(
        "--stop-after",
        choices=PIPELINE_STEPS,
        default=None,
    )

    parser.add_argument("--raw-dir", default="/root/dire/data/Illumina/sequences")
    parser.add_argument("--work-dir", default="/root/sequencing/intermediate/")
    parser.add_argument("--storage-dir", default="/root/dire/data/Analyser/processed/python")

    parser.add_argument("--main-env", default="dire_all")
    parser.add_argument("--confindr-env", default="confindr")
    parser.add_argument("--ariba-env", default="ariba")

    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--memory", type=int, default=32)
    parser.add_argument("--kmers", default="21,33,55,77,99,127")
    parser.add_argument("--min-contig-len", type=int, default=500)
    parser.add_argument("--organism", default="Escherichia")

    parser.add_argument("--skip-confindr", action="store_true")
    parser.add_argument("--skip-ariba", action="store_true")
    parser.add_argument("--confindr-db", default="/root/.confindr_db")
    parser.add_argument("--confindr-threads", type=int, default=4)
    parser.add_argument("--confindr-base-cutoff", type=int, default=3)
    parser.add_argument("--ariba-db", default="/root/resfinder_ariba_db")

    args = parser.parse_args()

    if args.stop_after is not None:
        if PIPELINE_STEPS.index(args.stop_after) < PIPELINE_STEPS.index(args.start_at):
            parser.error("--stop-after cannot be before --start-at")

    samples = select_samples(args, parser)

    print("\nSamples to process:")
    print(samples)
    print(f"Start at: {args.start_at}")
    print(f"Stop after: {args.stop_after}")
    print(f"Main env: {args.main_env}")
    print(f"ConFindr env: {args.confindr_env}")
    print(f"ARIBA env: {args.ariba_env}")
    print(f"Threads per sample: {args.threads}")
    print(f"SPAdes memory limit: {args.memory} GB")
    print(f"SPAdes k-mers: {args.kmers}")
    print(f"Run ConFindr: {not args.skip_confindr}")
    print(f"Run ARIBA: {not args.skip_ariba}")

    failed = []

    for sample in samples:
        try:
            run_sample(sample, args)
        except Exception as exc:
            print(f"ERROR in sample{sample}: {exc}", file=sys.stderr)
            failed.append(sample)

    if failed:
        print("\nFailed samples:")
        print(failed)
        raise SystemExit(1)

    print("\nAll selected samples processed successfully.")


if __name__ == "__main__":
    main()