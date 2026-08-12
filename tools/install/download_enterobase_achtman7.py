#!/usr/bin/env python3

from pathlib import Path
from urllib.request import urlopen, Request
from urllib.error import URLError, HTTPError
import gzip
import shutil
import argparse
import sys


BASE_URL = "https://enterobase.warwick.ac.uk/schemes/Escherichia.Achtman7GeneMLST"

FILES = [
    "MLST_Achtman_ref.fasta",
    "adk.fasta.gz",
    "fumC.fasta.gz",
    "gyrB.fasta.gz",
    "icd.fasta.gz",
    "mdh.fasta.gz",
    "profiles.list.gz",
    "purA.fasta.gz",
    "recA.fasta.gz",
]


def download_file(url: str, destination: Path, force: bool = False) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)

    if destination.exists() and destination.stat().st_size > 0 and not force:
        print(f"Skipping existing file: {destination}")
        return

    print(f"Downloading: {url}")
    request = Request(url, headers={"User-Agent": "Mozilla/5.0"})

    try:
        with urlopen(request, timeout=60) as response, destination.open("wb") as out:
            shutil.copyfileobj(response, out)
    except HTTPError as e:
        raise RuntimeError(f"HTTP error while downloading {url}: {e.code} {e.reason}") from e
    except URLError as e:
        raise RuntimeError(f"URL error while downloading {url}: {e.reason}") from e

    if destination.stat().st_size == 0:
        raise RuntimeError(f"Downloaded file is empty: {destination}")


def read_gzip_text(path: Path):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            yield line.rstrip("\n")


def check_profile(profiles_file: Path, st: str, expected_profile: list[str]) -> bool:
    for line in read_gzip_text(profiles_file):
        if not line or line.startswith("#"):
            continue

        parts = line.split()
        if parts[0] == st:
            found_profile = parts[1:8]
            if found_profile == expected_profile:
                print(f"OK: ST{st} profile found: {' '.join(parts)}")
                return True

            print(f"ERROR: ST{st} found, but profile differs.")
            print(f"Expected: {' '.join([st, *expected_profile])}")
            print(f"Found:    {' '.join(parts)}")
            return False

    print(f"ERROR: ST{st} not found in {profiles_file}")
    return False


def fasta_contains_allele(fasta_gz: Path, allele_number: str) -> bool:
    possible_headers = [
        f">{allele_number}",
        f">gyrB_{allele_number}",
        f">mdh_{allele_number}",
        f">{allele_number} ",
        f">{allele_number}\t",
    ]

    for line in read_gzip_text(fasta_gz):
        if not line.startswith(">"):
            continue
        if allele_number in line:
            print(f"OK: allele {allele_number} found in {fasta_gz.name}: {line}")
            return True
        if any(line.startswith(header) for header in possible_headers):
            print(f"OK: allele {allele_number} found in {fasta_gz.name}: {line}")
            return True

    print(f"WARNING: allele {allele_number} not found by simple header search in {fasta_gz.name}")
    return False


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Download current EnteroBase Escherichia Achtman 7-gene MLST scheme files."
    )
    parser.add_argument(
        "--outdir",
        default="enterobase_achtman7",
        help="Output directory for downloaded scheme files.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-download files even if they already exist.",
    )
    parser.add_argument(
        "--skip-checks",
        action="store_true",
        help="Skip verification of ST19543/ST19544 and key alleles.",
    )

    args = parser.parse_args()
    outdir = Path(args.outdir)

    for filename in FILES:
        url = f"{BASE_URL}/{filename}"
        destination = outdir / filename
        download_file(url, destination, force=args.force)

    if args.skip_checks:
        print("\nDownloaded files without validation checks.")
        return 0

    print("\nValidating key profiles and alleles...")

    profiles_file = outdir / "profiles.list.gz"
    gyrb_file = outdir / "gyrB.fasta.gz"
    mdh_file = outdir / "mdh.fasta.gz"

    checks = []

    checks.append(
        check_profile(
            profiles_file,
            "19543",
            ["6", "65", "32", "26", "1932", "8", "2"],
        )
    )
    checks.append(
        check_profile(
            profiles_file,
            "19544",
            ["37", "38", "1888", "37", "17", "11", "26"],
        )
    )

    checks.append(fasta_contains_allele(mdh_file, "1932"))
    checks.append(fasta_contains_allele(gyrb_file, "1888"))

    if all(checks):
        print("\nAll checks passed.")
        return 0

    print("\nDownload completed, but one or more validation checks failed or gave warnings.")
    return 1


if __name__ == "__main__":
    sys.exit(main())