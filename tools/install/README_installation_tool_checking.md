# DIRE installation and tool checking

Condensed instructions for installing the DIRE conda environments and documenting tool, database, computer, and operating system versions.

## Files

- `install_all.sh` installs or updates Miniconda, creates or updates the required conda environments, updates selected databases, validates the installation, and writes installation logs and system information.
- `print_tool_versions.py` is a read-only reporting script. It prints tool versions, conda environment information, database/reference information, and computer/OS information to timestamped report files.

## Required environment files

Run the installation script from a directory containing:

```text
dire_all.yml
confindr.yml
ariba.yml
checkm.yml
```

The installation script checks that these files exist before updating the environments.

## Conda environments

The pipeline uses four conda environments:

| Environment | Main purpose | Main tools |
|---|---|---|
| `dire_all` | Main pipeline tools | Trim Galore, Cutadapt, FastQC, SPAdes, seqkit, QUAST, AMRFinderPlus, MultiQC |
| `confindr` | Contamination screening | ConFindr, BBDuk/BBMap, KMA, Samtools, Java |
| `ariba` | Resistance-gene read mapping | ARIBA, Bowtie2, Samtools |
| `checkm` | Assembly completeness/contamination QC | CheckM, pplacer, Prodigal, HMMER |

## Installation

Run:

```bash
bash install_all.sh
```

The script:

1. Installs Miniconda if missing.
2. Accepts required Anaconda terms of service.
3. Updates conda.
4. Creates or updates the four conda environments from the YAML files using `conda env update --prune`.
5. Updates the AMRFinderPlus database.
6. Downloads and prepares the ARIBA ResFinder database.
7. Writes system information before and after installation.
8. Exports conda environment definitions and package lists.
9. Validates key tools by running version commands.
10. Writes a summary and log file.

## Installation output

Installation records are written to:

```text
/root/dire/data/Analyser/processed/python/installation
```

Generated files include the computer name and timestamp in the filename, for example:

```text
install_all_<computer>_<YYYYMMDD_HHMMSS>.log
system_info_before_install_<computer>_<YYYYMMDD_HHMMSS>.txt
system_info_after_install_<computer>_<YYYYMMDD_HHMMSS>.txt
install_summary_<computer>_<YYYYMMDD_HHMMSS>.txt
input_files_<computer>_<YYYYMMDD_HHMMSS>.sha256.txt
dire_all_export_<computer>_<YYYYMMDD_HHMMSS>.yml
dire_all_conda_list_<computer>_<YYYYMMDD_HHMMSS>.txt
```

Equivalent export and package-list files are also written for `confindr`, `ariba`, and `checkm`.

## Important database behavior

`install_all.sh` updates databases.

It runs AMRFinderPlus database update:

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

Use the installation script only when you intentionally want to create/update environments and databases.

## Tool and version report

Run:

```bash
python print_tool_versions.py
```

The script is read-only. It should not update databases or modify environments.

By default it writes reports to:

```text
/root/dire/data/Analyser/processed/python
```

Default report filenames include computer name and timestamp:

```text
tool_versions_<computer>_<YYYYMMDD_HHMMSS>.txt
tool_versions_summary_<computer>_<YYYYMMDD_HHMMSS>.txt
```

Use a fixed output filename if needed:

```bash
python print_tool_versions.py \
  -o /root/dire/data/Analyser/processed/python/tool_versions.txt
```

Use a fixed summary filename if needed:

```bash
python print_tool_versions.py \
  --summary-output /root/dire/data/Analyser/processed/python/tool_versions_summary.txt
```

## What the versioning script records

`print_tool_versions.py` records:

- computer name, OS, kernel, CPU, memory, and disk information
- base Python and conda information
- conda environment list
- tool versions from all four conda environments
- installation script and YAML file snapshots, with SHA256 checksums
- conda package lists
- conda environment exports
- database/reference information
- a short summary report

## AMRFinderPlus database version

The script reports the local AMRFinderPlus database version from:

```text
/root/miniconda3/envs/dire_all/share/amrfinderplus/data/latest/version.txt
```

Manual check:

```bash
cat /root/miniconda3/envs/dire_all/share/amrfinderplus/data/latest/version.txt
```

This only reads the installed local database version. It does not update the database.

## ResFinder ARIBA database version

The script checks:

```text
/root/resfinder_ariba_db/RESFINDER_DB_VERSION.txt
```

and lists the top-level files in:

```text
/root/resfinder_ariba_db
```

## ConFindr database information

The script also summarizes:

```text
/root/.confindr_db
```

## Useful options

Skip conda package lists:

```bash
python print_tool_versions.py --skip-conda-list
```

Skip conda environment exports:

```bash
python print_tool_versions.py --skip-conda-export
```

Skip database/reference summaries:

```bash
python print_tool_versions.py --skip-database-info
```

Add another database/reference path to summarize:

```bash
python print_tool_versions.py --db-path my_database=/path/to/database
```

Include additional environment YAML files in the report:

```bash
python print_tool_versions.py --env-yml extra_environment.yml
```

Use a longer timeout for slow commands:

```bash
python print_tool_versions.py --timeout 120
```

## Recommended workflow

For a new installation or environment update:

```bash
bash install_all.sh
python print_tool_versions.py
```

For documentation before manuscript submission, without updating anything:

```bash
python print_tool_versions.py
```

Archive both the detailed and summary reports together with the pipeline output.
