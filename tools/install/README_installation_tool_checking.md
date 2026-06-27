# DIRE installation and version reporting

Condensed instructions for installing the DIRE conda environments and documenting software, database, computer, operating system, and R package versions.

## Files

- `install_all.sh` installs or updates Miniconda, creates or updates the required conda environments, updates selected databases, validates the installation, and writes installation logs and system information.
- `collect_software_versions.py` is a read-only reporting script for the Linux/WSL pipeline environment. It records software and local database versions from the DIRE conda environments into a timestamped semicolon-separated CSV file.
- `collect_r_package_versions.py` is a read-only reporting script for the R analysis environment. It scans the R project for used packages, records installed R package versions, writes R session information, and optionally creates an Excel workbook.

## Required environment files

Run the installation script from a directory containing the conda environment YAML files used by the project, typically:

```text
dire_all.yml
confindr.yml
ariba.yml
checkm.yml
mlst.yml
```

The installation script should check that the required files exist before updating the environments.

## Conda environments

The pipeline/version reporting setup uses these conda environments:

| Environment | Main purpose | Main tools |
|---|---|---|
| `dire_all` | Main sample pipeline and summaries | Trim Galore, Cutadapt, FastQC, SPAdes, seqkit, QUAST, AMRFinderPlus, MultiQC |
| `confindr` | Contamination screening | ConFindr, BBDuk/BBMap, KMA, Samtools, Java |
| `ariba` | Resistance-gene read mapping | ARIBA, Bowtie2, Samtools |
| `checkm` | Assembly completeness/contamination QC | CheckM, pplacer, Prodigal, HMMER |
| `mlst` | Sequence typing | mlst |
| `juan` | Project environment for running AI model | Python |

`collect_software_versions.py` records all of the above when the corresponding environment directory exists. If an environment is missing, the row is kept in the output CSV with status `environment not found`.

## Installation

Run:

```bash
bash install_all.sh
```

The installation script is expected to:

1. Install Miniconda if missing.
2. Accept required Anaconda terms of service.
3. Update conda.
4. Create or update the required conda environments from the YAML files, typically using `conda env update --prune`.
5. Update the AMRFinderPlus database.
6. Download and prepare the ARIBA ResFinder database.
7. Write system information before and after installation.
8. Export conda environment definitions and package lists.
9. Validate key tools by running version commands.
10. Write a summary and log file.

Use the installation script only when you intentionally want to create/update environments and databases.

## Installation output

Installation records are written to:

```text
/root/dire/data/Analyser/processed/python/installation
```

Generated files typically include the computer name and timestamp in the filename, for example:

```text
install_all_<computer>_<YYYYMMDD_HHMMSS>.log
system_info_before_install_<computer>_<YYYYMMDD_HHMMSS>.txt
system_info_after_install_<computer>_<YYYYMMDD_HHMMSS>.txt
install_summary_<computer>_<YYYYMMDD_HHMMSS>.txt
input_files_<computer>_<YYYYMMDD_HHMMSS>.sha256.txt
dire_all_export_<computer>_<YYYYMMDD_HHMMSS>.yml
dire_all_conda_list_<computer>_<YYYYMMDD_HHMMSS>.txt
```

Equivalent export and package-list files may also be written for `confindr`, `ariba`, `checkm`, and `mlst`.

## Important database behavior

`install_all.sh` updates databases.

It runs the AMRFinderPlus database update:

```bash
conda run -n dire_all --no-capture-output amrfinder_update --force_update
```

It also rebuilds the ARIBA ResFinder database at:

```text
/root/resfinder_ariba_db
```

The prepared ResFinder database records a local preparation file:

```text
/root/resfinder_ariba_db/RESFINDER_DB_VERSION.txt
```

## Linux/WSL software and database version report

Run:

```bash
python collect_software_versions.py
```

The script is read-only. It should not update databases or modify environments.

By default it writes a timestamped semicolon-separated CSV to:

```text
/root/dire/data/Analyser/processed/python
```

Default filename pattern:

```text
software_versions_<computer>_<YYYYMMDD_HHMMSS>.csv
```

Use a custom output directory if needed:

```bash
python collect_software_versions.py   --output-dir /root/dire/data/Analyser/processed/python
```

## What `collect_software_versions.py` records

The output CSV contains one row per software/database item with these columns:

```text
collected_at
host
environment
environment_path
software
version
command
return_code
status
raw_output
```

The script currently checks:

| Environment | Recorded tools |
|---|---|
| `base` | Python, conda |
| `dire_all` | Python, Trim Galore, Cutadapt, FastQC, SPAdes, QUAST, MultiQC, SeqKit, AMRFinderPlus |
| `confindr` | Python, ConFindr, BBDuk, KMA, Samtools, SeqKit |
| `ariba` | Python, ARIBA, Bowtie2, Samtools |
| `checkm` | Python, CheckM, pplacer, Prodigal, HMMER |
| `mlst` | Python, mlst |
| `juan` | Python |

The script also records local database version files for:

| Database | Checked path |
|---|---|
| AMRFinderPlus database | `/root/miniconda3/envs/dire_all/share/amrfinderplus/data/*/version.txt` and `/root/miniconda3/envs/dire_all/bin/data/*/version.txt` |
| ResFinder ARIBA database | `/root/resfinder_ariba_db/RESFINDER_DB_VERSION.txt` |

Rows with problems are printed at the end of the run. A failing row does not necessarily mean the whole report failed; it can indicate a missing optional environment, a version command that returned a non-zero exit code, or a missing local database version file.

## Manual database checks

AMRFinderPlus database:

```bash
cat /root/miniconda3/envs/dire_all/share/amrfinderplus/data/latest/version.txt
```

ResFinder ARIBA database:

```bash
cat /root/resfinder_ariba_db/RESFINDER_DB_VERSION.txt
ls -lah /root/resfinder_ariba_db
```

ConFindr database:

```bash
ls -lah /root/.confindr_db
```

## R package version report

Run this on the system where the R analysis environment is installed:

```bash
python collect_r_package_versions.py
```

The script is read-only. It scans R files for package use and records the installed package versions from R.

By default it scans:

```text
~/OneDrive - Västra Götalandsregionen/git/dire/R
```

and writes output to:

```text
~/OneDrive - Västra Götalandsregionen/DIRE/Analyser/processed/R
```

Default output files:

```text
r_used_packages_<computer>_<YYYYMMDD_HHMMSS>.csv
r_installed_packages_<computer>_<YYYYMMDD_HHMMSS>.csv
r_session_info_<computer>_<YYYYMMDD_HHMMSS>.txt
r_package_versions_<computer>_<YYYYMMDD_HHMMSS>.xlsx
```

The CSV files are semicolon-separated.

## What `collect_r_package_versions.py` records

`collect_r_package_versions.py` writes:

- a table of packages used in the R project scripts, based on `library()`, `require()`, and `pkg::function` usage
- installed R package versions from `installed.packages()`
- R session information from `sessionInfo()`
- an optional Excel workbook with sheets for used and installed packages, if `openpyxl` is available

The used-package table includes:

```text
package
version
status
lib_path
built
source_files
```

The installed-package table includes package metadata such as version, library path, dependencies/imports, license, and R build version.

## Useful R version-report options

Use a custom R project directory:

```bash
python collect_r_package_versions.py   --project-dir "C:/Users/xhessm/OneDrive - Västra Götalandsregionen/git/dire/R"
```

Use a custom output directory:

```bash
python collect_r_package_versions.py   --output-dir "C:/Users/xhessm/OneDrive - Västra Götalandsregionen/DIRE/Analyser/processed/R"
```

Use a specific `Rscript.exe`:

```bash
python collect_r_package_versions.py   --rscript "C:/Program Files/R/R-4.5.2/bin/Rscript.exe"
```

If `Rscript` is not on `PATH`, the script also checks common Windows R installation paths.

## Recommended workflow

For a new installation or environment/database update:

```bash
bash install_all.sh
python collect_software_versions.py
```

For R package documentation:

```bash
python collect_r_package_versions.py
```

For documentation before manuscript submission, without updating environments or databases:

```bash
python collect_software_versions.py
python collect_r_package_versions.py
```

Archive the software-version CSV, R package CSV files, R session information file, optional R Excel workbook, and installation logs together with the pipeline output.
