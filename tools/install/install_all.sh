#!/usr/bin/env bash

set -euo pipefail

CONDA_DIR="$HOME/miniconda3"

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
conda env update -f dire_all.yml --prune
conda env update -f confindr.yml --prune
conda env update -f ariba.yml --prune

echo "Updating AMRFinder database..."
conda run -n dire_all amrfinder_update --force_update

echo "Installing/updating ARIBA ResFinder database..."
rm -rf /root/resfinder_ariba_db
mkdir -p /root/resfinder_ariba_db
conda run -n ariba ariba getref resfinder /root/resfinder_ariba_db
conda run -n ariba --no-capture-output ariba prepareref \
  -f /root/resfinder_ariba_db.fa \
  -m /root/resfinder_ariba_db.tsv \
  /root/resfinder_ariba_db

echo "ResFinder database prepared: $(date)" \
    > /root/resfinder_ariba_db/RESFINDER_DB_VERSION.txt

echo "Done."