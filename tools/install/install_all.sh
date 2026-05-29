#!/usr/bin/env bash

set -euo pipefail

CONDA_DIR="$HOME/miniconda3"
INSTALL_OUT_DIR="/root/dire/data/Analyser/processed/python/installation"
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
COMPUTER_NAME="$(hostname -s 2>/dev/null || hostname 2>/dev/null || echo unknown_host)"
COMPUTER_NAME="$(printf '%s' "$COMPUTER_NAME" | tr -c 'A-Za-z0-9_.-' '_')"
RUN_TAG="${COMPUTER_NAME}_${TIMESTAMP}"
LOG_FILE="$INSTALL_OUT_DIR/install_all_${RUN_TAG}.log"
SYSTEM_INFO_BEFORE="$INSTALL_OUT_DIR/system_info_before_install_${RUN_TAG}.txt"
SYSTEM_INFO_AFTER="$INSTALL_OUT_DIR/system_info_after_install_${RUN_TAG}.txt"
INSTALL_SUMMARY="$INSTALL_OUT_DIR/install_summary_${RUN_TAG}.txt"
INPUT_CHECKSUMS="$INSTALL_OUT_DIR/input_files_${RUN_TAG}.sha256.txt"
INPUT_SNAPSHOT_DIR="$INSTALL_OUT_DIR/input_files_${RUN_TAG}"
ENV_FILES=(dire_all.yml confindr.yml ariba.yml checkm.yml)
ENV_NAMES=(dire_all confindr ariba checkm)

mkdir -p "$INSTALL_OUT_DIR"
exec > >(tee -a "$LOG_FILE") 2>&1
set -x

write_system_info() {
    local outfile="$1"

    {
        echo "DIRE installation system information"
        echo "Created: $(date --iso-8601=seconds 2>/dev/null || date)"
        echo "Run tag: $RUN_TAG"
        echo "Computer name: $COMPUTER_NAME"
        echo "Working directory: $(pwd)"
        echo "User: $(whoami)"
        echo "Host: $(hostname)"
        echo "Shell: ${SHELL:-unknown}"
        echo "PATH: ${PATH:-unknown}"
        echo ""

        echo "========================================"
        echo "Operating system"
        echo "========================================"
        if [ -f /etc/os-release ]; then
            cat /etc/os-release
        else
            echo "/etc/os-release not found"
        fi
        echo ""
        uname -a || true
        hostnamectl 2>/dev/null || true
        echo ""

        echo "========================================"
        echo "CPU"
        echo "========================================"
        nproc || true
        lscpu 2>/dev/null || true
        echo ""

        echo "========================================"
        echo "Memory"
        echo "========================================"
        free -h 2>/dev/null || true
        grep -E 'MemTotal|MemFree|MemAvailable|SwapTotal|SwapFree' /proc/meminfo 2>/dev/null || true
        echo ""

        echo "========================================"
        echo "Disk"
        echo "========================================"
        df -h || true
        echo ""
        lsblk 2>/dev/null || true
        echo ""

        echo "========================================"
        echo "Conda before/after availability"
        echo "========================================"
        command -v conda || true
        conda --version 2>/dev/null || true
        conda info 2>/dev/null || true
    } > "$outfile"
}

require_files() {
    for f in "${ENV_FILES[@]}"; do
        if [ ! -f "$f" ]; then
            echo "Missing required environment file: $f" >&2
            exit 1
        fi
    done
}

snapshot_install_inputs() {
    mkdir -p "$INPUT_SNAPSHOT_DIR"

    if [ -f "$0" ]; then
        cp -f "$0" "$INPUT_SNAPSHOT_DIR/$(basename "$0")"
    fi

    for f in "${ENV_FILES[@]}"; do
        cp -f "$f" "$INPUT_SNAPSHOT_DIR/$f"
    done

    sha256sum "$INPUT_SNAPSHOT_DIR"/* > "$INPUT_CHECKSUMS" 2>/dev/null || true
}

export_conda_envs() {
    for env in "${ENV_NAMES[@]}"; do
        conda env export -n "$env" --no-builds > "$INSTALL_OUT_DIR/${env}_export_${RUN_TAG}.yml"
        conda list -n "$env" > "$INSTALL_OUT_DIR/${env}_conda_list_${RUN_TAG}.txt"
    done
}

validate_installation() {
    echo "Validating installation..."

    conda run -n dire_all --no-capture-output python --version
    conda run -n dire_all --no-capture-output trim_galore --version
    conda run -n dire_all --no-capture-output cutadapt --version
    conda run -n dire_all --no-capture-output fastqc --version
    conda run -n dire_all --no-capture-output spades.py --version
    conda run -n dire_all --no-capture-output seqkit version
    conda run -n dire_all --no-capture-output quast.py --version
    conda run -n dire_all --no-capture-output amrfinder --version
    conda run -n dire_all --no-capture-output multiqc --version

    conda run -n confindr --no-capture-output python --version
    conda run -n confindr --no-capture-output confindr.py --version
    conda run -n confindr --no-capture-output bbduk.sh --version

    conda run -n ariba --no-capture-output python --version
    conda run -n ariba --no-capture-output ariba version

    conda run -n checkm --no-capture-output python --version
    conda run -n checkm --no-capture-output checkm version || true
}

write_install_summary() {
    {
        echo "DIRE installation summary"
        echo "Completed: $(date --iso-8601=seconds 2>/dev/null || date)"
        echo "Run tag: $RUN_TAG"
        echo "Computer name: $COMPUTER_NAME"
        echo "Host: $(hostname)"
        echo "User: $(whoami)"
        echo "Working directory: $(pwd)"
        echo "Conda base: $(conda info --base)"
        echo "Log file: $LOG_FILE"
        echo "System info before install: $SYSTEM_INFO_BEFORE"
        echo "System info after install: $SYSTEM_INFO_AFTER"
        echo "Environment exports: $INSTALL_OUT_DIR/*_export_${RUN_TAG}.yml"
        echo "Conda package lists: $INSTALL_OUT_DIR/*_conda_list_${RUN_TAG}.txt"
        echo "Input file snapshot directory: $INPUT_SNAPSHOT_DIR"
        echo "Input file checksums: $INPUT_CHECKSUMS"
    } > "$INSTALL_SUMMARY"
}

write_system_info "$SYSTEM_INFO_BEFORE"
require_files
snapshot_install_inputs

# Install Miniconda only if missing
if ! command -v conda >/dev/null 2>&1; then
    echo "Conda not found. Installing Miniconda..."

    mkdir -p "$CONDA_DIR"

    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
        -O "$CONDA_DIR/miniconda.sh"

    bash "$CONDA_DIR/miniconda.sh" -b -u -p "$CONDA_DIR"

    rm -f "$CONDA_DIR/miniconda.sh"

    source "$CONDA_DIR/etc/profile.d/conda.sh"
    conda init bash
else
    echo "Conda already installed."
    source "$(conda info --base)/etc/profile.d/conda.sh"
fi

conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Update conda
conda update -n base -c defaults conda -y

# Create/update conda envs
for f in "${ENV_FILES[@]}"; do
    conda env update -f "$f" --prune
done

echo "Updating AMRFinder database..."
conda run -n dire_all --no-capture-output amrfinder_update --force_update

echo "Installing/updating ARIBA ResFinder database..."
echo "Removing existing ARIBA ResFinder database at /root/resfinder_ariba_db"
( 
cd ~
conda run -n ariba --no-capture-output ariba getref resfinder /root/resfinder_ariba_db_download
rm -rf /root/resfinder_ariba_db
conda run -n ariba --no-capture-output ariba prepareref \
  -f /root/resfinder_ariba_db_download.fa \
  -m /root/resfinder_ariba_db_download.tsv \
  /root/resfinder_ariba_db
)
{
    echo "ResFinder database prepared: $(date --iso-8601=seconds 2>/dev/null || date)"
    echo "Run tag: $RUN_TAG"
    echo "Computer name: $COMPUTER_NAME"
} > /root/resfinder_ariba_db/RESFINDER_DB_VERSION.txt

write_system_info "$SYSTEM_INFO_AFTER"
export_conda_envs
validate_installation
write_install_summary

echo "Done."
echo "Run tag: $RUN_TAG"
echo "Installation log: $LOG_FILE"
echo "Installation summary: $INSTALL_SUMMARY"
