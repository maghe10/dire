# DIRE sequencing pipeline

Condensed overview of the DIRE sample-level sequencing pipeline and project-level postpipeline.

## Scripts

- `pipeline_full.py`: runs sample-level processing from raw reads to trimmed reads, assemblies, QC outputs, AMR gene calls, contamination screening, and ARIBA/ResFinder results.
- `postpipeline_full.py`: runs project-level processing after sample-level analysis, including CheckM, CheckM result copying, ConFindr and QUAST summary tables, MLST, and MultiQC.

## Conda environments

| Environment | Used for |
|---|---|
| `dire_all` | Trim Galore, FastQC, SPAdes, seqkit, QUAST, AMRFinderPlus, MultiQC, summary scripts |
| `confindr` | ConFindr contamination screening |
| `ariba` | ARIBA / ResFinder analysis |
| `checkm` | CheckM lineage workflow and QA export |
| `mlst` | MLST typing of filtered assemblies |

## Sample naming convention

The current pipeline uses three-digit sample identifiers throughout processed outputs, for example `sample002`, `sample021`, and `sample065`.

Raw read files are expected to use the DIRE input naming convention:

```text
DIRE_EC_<SAMPLE>_R1.fastq.gz
DIRE_EC_<SAMPLE>_R2.fastq.gz
```

where `<SAMPLE>` is normally zero-padded, for example:

```text
DIRE_EC_002_R1.fastq.gz
DIRE_EC_002_R2.fastq.gz
```

When running a single sample, prefer passing the zero-padded ID:

```bash
python pipeline_full.py -s 002
```

The built-in batches already use zero-padded IDs.

Trim Galore is run with `--basename sample<SAMPLE>`. Therefore, the current expected trimmed-read names are:

```text
sample<SAMPLE>_val_1.fq.gz
sample<SAMPLE>_val_2.fq.gz
```

For example:

```text
sample002_val_1.fq.gz
sample002_val_2.fq.gz
```

## Main pipeline

The main pipeline processes one sample, one predefined batch, or all samples.

Run one sample:

```bash
python pipeline_full.py --sample 002
```

Short option:

```bash
python pipeline_full.py -s 002
```

Run one batch:

```bash
python pipeline_full.py --batch A
```

Short option:

```bash
python pipeline_full.py -b A
```

Run all samples:

```bash
python pipeline_full.py --all
```

When `--all` is used, the script asks for confirmation before starting all samples.

Resume from a specific step:

```bash
python pipeline_full.py --sample 002 --start-at spades
```

Stop after a specific step:

```bash
python pipeline_full.py --sample 002 --stop-after quast
```

Run a step range:

```bash
python pipeline_full.py --sample 002 --start-at filter --stop-after amrfinder
```

Force rerun:

```bash
python pipeline_full.py --sample 002 --force
```

Skip optional steps:

```bash
python pipeline_full.py --sample 002 --skip-confindr
python pipeline_full.py --sample 002 --skip-ariba
```

Custom resources and databases:

```bash
python pipeline_full.py --sample 002 --threads 16 --memory 64
python pipeline_full.py --sample 002 --confindr-db /path/to/confindr_db --ariba-db /path/to/resfinder_ariba_db
```

## Main pipeline sample batches

| Batch | Samples |
|---|---|
| A | 002, 003, 004, 005, 006, 008 |
| B | 009, 011, 012, 013, 018, 019, 022 |
| C | 024, 025, 026, 028, 029, 034, 035, 036 |
| D | 039, 040, 043, 044, 045, 047, 048, 050 |
| E | 051, 052, 056, 057, 059, 061, 062, 064 |
| F | 066, 067, 069, 070, 072, 076, 080, 082 |
| G | 083, 084, 086, 087, 088, 089, 091, 093 |
| H | 094, 096, 097, 098, 101, 102, 103, 105 |
| I | 109, 111, 112, 113, 120 |
| J | 017, 020, 023, 027, 032, 037, 041, 046 |
| K | 049, 053, 054, 058, 075, 095, 100, 115, 116 |
| L | 021, 063, 077, 078, 114 |
| M | 055, 060, 065, 079, 085, 110, 121, 122, 123, 124 |
| N | 125 |

## Main pipeline steps

Steps run in this order:

1. `trimgalore` - trims paired-end reads and runs FastQC.
2. `spades` - assembles reads using SPAdes isolate mode.
3. `filter` - filters assembled scaffolds by minimum contig length using seqkit.
4. `quast` - runs assembly quality assessment.
5. `amrfinder` - detects antimicrobial resistance genes using AMRFinderPlus.
6. `confindr` - screens reads for contamination unless `--skip-confindr` is used.
7. `ariba` - runs ARIBA against the ResFinder database unless `--skip-ariba` is used.

## Main pipeline default options

| Option | Default |
|---|---|
| `--raw-dir` | `/root/dire/data/Illumina/sequences` |
| `--work-dir` | `/root/sequencing/intermediate/` |
| `--storage-dir` | `/root/dire/data/Analyser/processed/python` |
| `--main-env` | `dire_all` |
| `--confindr-env` | `confindr` |
| `--ariba-env` | `ariba` |
| `--threads` | `8` |
| `--memory` | `32` |
| `--kmers` | `21,33,55,77,99,127` |
| `--min-contig-len` | `500` |
| `--organism` | `Escherichia` |
| `--confindr-db` | `/root/.confindr_db` |
| `--confindr-threads` | `4` |
| `--confindr-base-cutoff` | `3` |
| `--ariba-db` | `/root/resfinder_ariba_db` |

## Main pipeline important outputs

| Output | Default location |
|---|---|
| Trimmed reads | `/root/sequencing/intermediate/trimmed/sample<SAMPLE>/` and `/root/dire/data/Analyser/processed/python/trimmed/sample<SAMPLE>/` |
| Expected trimmed FASTQ files | `sample<SAMPLE>_val_1.fq.gz` and `sample<SAMPLE>_val_2.fq.gz` |
| SPAdes output | `/root/sequencing/intermediate/spades_standard/sample<SAMPLE>/` and `/root/dire/data/Analyser/processed/python/spades/standard/sample<SAMPLE>/` |
| Filtered assembly | `/root/sequencing/intermediate/assembly/spades_standard/sample<SAMPLE>.fasta` and `/root/dire/data/Analyser/processed/python/assembly/spades_standard/sample<SAMPLE>.fasta` |
| QUAST report | `/root/sequencing/intermediate/quast/sample<SAMPLE>/` and `/root/dire/data/Analyser/processed/python/quast/sample<SAMPLE>/` |
| AMRFinderPlus result | `/root/sequencing/intermediate/amrfinder/spades_standard/sample<SAMPLE>.tsv` and `/root/dire/data/Analyser/processed/python/amrfinder/spades_standard/sample<SAMPLE>.tsv` |
| ConFindr result | `/root/sequencing/intermediate/confindr/sample<SAMPLE>/` and `/root/dire/data/Analyser/processed/python/confindr/sample<SAMPLE>/` |
| ARIBA result | `/root/sequencing/intermediate/ariba/sample<SAMPLE>/` and `/root/dire/data/Analyser/processed/python/ariba/sample<SAMPLE>/` |
| Logs | `/root/sequencing/intermediate/logs/` |

## Postpipeline

The postpipeline runs project-level steps after the main pipeline has produced filtered assemblies and sample-level output directories.

Run the full postpipeline:

```bash
python postpipeline_full.py
```

Run from a specific step:

```bash
python postpipeline_full.py --start-at checkm_qa
```

Stop after a specific step:

```bash
python postpipeline_full.py --stop-after summarize_quast
```

Run a step range:

```bash
python postpipeline_full.py --start-at summarize_confindr --stop-after multiqc
```

Run only MLST:

```bash
python postpipeline_full.py --start-at mlst --stop-after mlst
```

Run only MultiQC:

```bash
python postpipeline_full.py --start-at multiqc --stop-after multiqc
```

Run only CheckM QA if lineage has already completed:

```bash
python postpipeline_full.py --start-at checkm_qa --stop-after checkm_qa
```

Run only CheckM result copy:

```bash
python postpipeline_full.py --start-at copy_checkm --stop-after copy_checkm
```

Force rerun:

```bash
python postpipeline_full.py --force
```

When `--force` is used and `checkm_lineage` is run, the CheckM intermediate directory is removed before CheckM starts.

## Postpipeline steps

Steps run in this order:

1. `checkm_lineage` - runs CheckM lineage workflow on all filtered assemblies.
2. `checkm_qa` - exports CheckM QA results both as a tab-separated table and as a CheckM `-o 1` file for MultiQC.
3. `copy_checkm` - copies the CheckM result directory to processed storage.
4. `copy_checkm_summary` - copies `qa.tsv` from stored CheckM output to `quality/checkm_summary.tsv`.
5. `summarize_confindr` - runs `summarize_confindr.py` and writes a project-level ConFindr summary table.
6. `summarize_quast` - runs `summarize_quast.py` and writes a project-level QUAST summary table.
7. `mlst` - runs MLST on all filtered assemblies and writes normalized tab-separated output.
8. `multiqc` - runs MultiQC on the processed trimmed-read directory and stored CheckM directory.

## Postpipeline default options

| Option | Default |
|---|---|
| `--main-env` | `dire_all` |
| `--checkm-env` | `checkm` |
| `--mlst-env` | `mlst` |
| `--trimmed-dir` | `/root/dire/data/Analyser/processed/python/trimmed` |
| `--multiqc-dir` | `/root/dire/data/Analyser/processed/python/multiqc` |
| `--assembly-dir` | `/root/sequencing/intermediate/assembly/spades_standard` |
| `--checkm-dir` | `/root/sequencing/intermediate/checkm` |
| `--storage-checkm-dir` | `/root/dire/data/Analyser/processed/python/checkm` |
| `--quality-dir` | `/root/dire/data/Analyser/processed/python/quality` |
| `--confindr-dir` | `/root/dire/data/Analyser/processed/python/confindr` |
| `--quast-dir` | `/root/dire/data/Analyser/processed/python/quast` |
| `--mlst-dir` | `/root/dire/data/Analyser/processed/python/mlst` |
| `--confindr-summary-output` | `/root/dire/data/Analyser/processed/python/quality/confindr_summary.tsv` |
| `--quast-summary-output` | `/root/dire/data/Analyser/processed/python/quality/quast_summary.tsv` |
| `--mlst-output` | `/root/dire/data/Analyser/processed/python/mlst/mlst_ecoli_achtman.tsv` |
| `--summarize-confindr-script` | `summarize_confindr.py` |
| `--summarize-quast-script` | `summarize_quast.py` |
| `--tools-dir` | `/root/dire/program/tools/bin` |
| `--log-dir` | `/root/sequencing/intermediate/logs` |
| `--checkm-threads` | `4` |
| `--assembly-extension` | `fasta` |
| `--mlst-scheme` | `ecoli_achtman_4` |
| `--exclude-mlst-sample` | `sample014`, `sample038` |

## Postpipeline important outputs

| Output | Default location |
|---|---|
| CheckM lineage and QA intermediate files | `/root/sequencing/intermediate/checkm/` |
| Stored CheckM output | `/root/dire/data/Analyser/processed/python/checkm/` |
| CheckM summary table | `/root/dire/data/Analyser/processed/python/quality/checkm_summary.tsv` |
| ConFindr summary table | `/root/dire/data/Analyser/processed/python/quality/confindr_summary.tsv` |
| QUAST summary table | `/root/dire/data/Analyser/processed/python/quality/quast_summary.tsv` |
| MLST raw output | `/root/dire/data/Analyser/processed/python/mlst/mlst_raw.tsv` |
| MLST normalized output | `/root/dire/data/Analyser/processed/python/mlst/mlst_ecoli_achtman.tsv` |
| MultiQC report | `/root/dire/data/Analyser/processed/python/multiqc/multiqc_report.html` |
| Logs | `/root/sequencing/intermediate/logs/` |

## MLST exclusions

The postpipeline excludes these samples from MLST by default:

```text
sample014
sample038
```

Add additional excluded samples with repeated `--exclude-mlst-sample` arguments:

```bash
python postpipeline_full.py --exclude-mlst-sample sample014 --exclude-mlst-sample sample038 --exclude-mlst-sample sample999
```

Because `--exclude-mlst-sample` uses `action="append"`, custom values are added to the defaults unless the script is edited.

## MultiQC default exclusions

MultiQC ignores these patterns by default:

```text
*/sample014*
*/sample038*
```

Add additional ignore patterns with repeated `--ignore` arguments:

```bash
python postpipeline_full.py --ignore "*/sample014*" --ignore "*/sample038*" --ignore "*/other_pattern*"
```

Because `--ignore` uses `action="append"`, custom patterns are added to the defaults unless the script is edited.

## Resume and rerun behavior

Both scripts skip steps when their expected output already exists, unless `--force` is used.

For the main pipeline, outputs are usually written first to the intermediate work directory and then copied to the processed storage directory.

For the postpipeline, CheckM output is generated in the intermediate directory, copied to processed storage, and then summarized into the processed `quality` directory. MLST output is written directly to the processed `mlst` directory.


## Typical full run

```bash
python pipeline_full.py --all
python postpipeline_full.py
```
