# DIRE sequencing pipeline

Condensed overview of the DIRE sample-level sequencing pipeline and project-level postpipeline.

## Scripts

- `pipeline_full.py`: runs sample-level processing from raw reads to trimmed reads, assemblies, QC outputs, AMR gene calls, contamination screening, and ARIBA/ResFinder results.
- `postpipeline_full.py`: runs project-level processing after sample-level analysis, including CheckM, CheckM result copying, ConFindr and QUAST summary tables, and MultiQC.

## Conda environments

| Environment | Used for |
|---|---|
| `dire_all` | Trim Galore, FastQC, SPAdes, seqkit, QUAST, AMRFinderPlus, MultiQC, summary scripts |
| `confindr` | ConFindr contamination screening |
| `ariba` | ARIBA / ResFinder analysis |
| `checkm` | CheckM lineage workflow and QA export |

## Main pipeline

The main pipeline processes one sample, one predefined batch, or all samples.

Run one sample:

```bash
python pipeline_full.py --sample 2
```

Short option:

```bash
python pipeline_full.py -s 2
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
python pipeline_full.py --sample 2 --start-at spades
```

Stop after a specific step:

```bash
python pipeline_full.py --sample 2 --stop-after quast
```

Run a step range:

```bash
python pipeline_full.py --sample 2 --start-at filter --stop-after amrfinder
```

Force rerun:

```bash
python pipeline_full.py --sample 2 --force
```

Skip optional steps:

```bash
python pipeline_full.py --sample 2 --skip-confindr
python pipeline_full.py --sample 2 --skip-ariba
```

Custom resources and databases:

```bash
python pipeline_full.py --sample 2 --threads 16 --memory 64
python pipeline_full.py --sample 2 --confindr-db /path/to/confindr_db --ariba-db /path/to/resfinder_ariba_db
```

## Main pipeline sample batches

| Batch | Samples |
|---|---|
| A | 2, 3, 4, 5, 6, 8 |
| B | 9, 11, 12, 13, 14, 18, 19, 22 |
| C | 24, 25, 26, 28, 29, 34, 35, 36 |
| D | 39, 40, 43, 44, 45, 47, 48, 50 |
| E | 51, 52, 56, 57, 59, 61, 62, 64 |
| F | 66, 67, 69, 70, 72, 76, 80, 82 |
| G | 83, 84, 86, 87, 88, 89, 91, 93 |
| H | 94, 96, 97, 98, 101, 102, 103, 105 |
| I | 109, 111, 112, 113, 120 |
| J | 17, 20, 23, 27, 32, 37, 38, 41, 46 |
| K | 49, 53, 54, 58, 75, 95, 100, 115, 116 |
| L | 21, 63, 77, 78, 114 |
| M | 121, 122, 123, 124, 55B, 60B, 65B, 79B, 85B, 110B |
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
| `--raw-dir` | `/root/sequencing/in/reads` |
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
7. `multiqc` - runs MultiQC on the processed trimmed-read directory and stored CheckM directory.

## Postpipeline default options

| Option | Default |
|---|---|
| `--main-env` | `dire_all` |
| `--checkm-env` | `checkm` |
| `--trimmed-dir` | `/root/dire/data/Analyser/processed/python/trimmed` |
| `--multiqc-dir` | `/root/dire/data/Analyser/processed/python/multiqc` |
| `--assembly-dir` | `/root/sequencing/intermediate/assembly/spades_standard` |
| `--checkm-dir` | `/root/sequencing/intermediate/checkm` |
| `--storage-checkm-dir` | `/root/dire/data/Analyser/processed/python/checkm` |
| `--quality-dir` | `/root/dire/data/Analyser/processed/python/quality` |
| `--confindr-dir` | `/root/dire/data/Analyser/processed/python/confindr` |
| `--quast-dir` | `/root/dire/data/Analyser/processed/python/quast` |
| `--confindr-summary-output` | `/root/dire/data/Analyser/processed/python/quality/confindr_summary.tsv` |
| `--quast-summary-output` | `/root/dire/data/Analyser/processed/python/quality/quast_summary.tsv` |
| `--summarize-confindr-script` | `summarize_confindr.py` |
| `--summarize-quast-script` | `summarize_quast.py` |
| `--tools-dir` | `/root/dire/program/tools/bin` |
| `--log-dir` | `/root/sequencing/intermediate/logs` |
| `--checkm-threads` | `4` |
| `--assembly-extension` | `fasta` |

## Postpipeline important outputs

| Output | Default location |
|---|---|
| CheckM lineage and QA intermediate files | `/root/sequencing/intermediate/checkm/` |
| Stored CheckM output | `/root/dire/data/Analyser/processed/python/checkm/` |
| CheckM summary table | `/root/dire/data/Analyser/processed/python/quality/checkm_summary.tsv` |
| ConFindr summary table | `/root/dire/data/Analyser/processed/python/quality/confindr_summary.tsv` |
| QUAST summary table | `/root/dire/data/Analyser/processed/python/quality/quast_summary.tsv` |
| MultiQC report | `/root/dire/data/Analyser/processed/python/multiqc/multiqc_report.html` |
| Logs | `/root/sequencing/intermediate/logs/` |

## MultiQC default exclusions

MultiQC ignores these patterns by default:

```text
*/sample14*
*/sample38*
```

Add additional ignore patterns with repeated `--ignore` arguments:

```bash
python postpipeline_full.py --ignore "*/sample14*" --ignore "*/sample38*" --ignore "*/other_pattern*"
```

Because `--ignore` uses `action="append"`, custom patterns are added to the defaults unless the script is edited.

## Resume and rerun behavior

Both scripts skip steps when their expected output already exists, unless `--force` is used.

For the main pipeline, outputs are usually written first to the intermediate work directory and then copied to the processed storage directory.

For the postpipeline, CheckM output is generated in the intermediate directory, copied to processed storage, and then summarized into the processed `quality` directory.

## Typical full run

```bash
python pipeline_full.py --all
python postpipeline_full.py
```

