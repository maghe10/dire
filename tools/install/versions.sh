#!/usr/bin/env bash

# use with:
#   bash versions.sh &> versions.log

set -euo pipefail
set -x

echo "========================================"
echo "base environment"
echo "========================================"

python --version
conda --version

echo ""
echo "========================================"
echo "dire_all environment"
echo "========================================"

conda run -n dire_all python --version

echo ""
echo "raw read processing"

conda run -n dire_all trim_galore --version
conda run -n dire_all cutadapt --version
conda run -n dire_all fastqc --version

echo ""
echo "assembly"

conda run -n dire_all spades.py --version
conda run -n dire_all seqkit version

echo ""
echo "post assembly tools"

conda run -n dire_all quast.py --version
conda run -n dire_all multiqc --version
conda run -n dire_all amrfinder --version

echo ""
echo "========================================"
echo "confindr environment"
echo "========================================"

conda run -n confindr python --version
conda run -n confindr confindr.py --version
conda run -n confindr bbduk.sh --version

echo ""
echo "KMA:"
conda run -n confindr bash -c "kma 2>&1 | head -n 3"

echo ""
echo "samtools:"
conda run -n confindr bash -c "samtools --version | head -n 3"

echo ""
echo "========================================"
echo "ariba environment"
echo "========================================"

conda run -n ariba python --version
conda run -n ariba ariba version

echo ""
echo "bowtie2:"
conda run -n ariba bash -c "bowtie2 --version | head -n 1"

echo ""
echo "samtools:"
conda run -n ariba bash -c "samtools --version | head -n 3"

echo ""
echo "========================================"
echo "database versions"
echo "========================================"

echo ""
echo "AMRFinder database:"
conda run -n dire_all amrfinder_update --version

echo ""
echo "ResFinder ARIBA database:"
if [ -f /root/resfinder_ariba_db/RESFINDER_DB_VERSION.txt ]; then
    cat /root/resfinder_ariba_db/RESFINDER_DB_VERSION.txt
else
    echo "No RESFINDER_DB_VERSION.txt found"
fi

echo ""
echo "ResFinder ARIBA database files:"
ls -lh /root/resfinder_ariba_db | head