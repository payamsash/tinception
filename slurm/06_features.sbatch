#!/bin/bash -l
# Extract all FreeSurfer feature tables + fsaverage5 vertex matrices -> derivatives/features/
# Submit from repo root (after slurm/06_annot_stats.sbatch): sbatch slurm/06_features.sbatch
#SBATCH --job-name=features
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --exclude=u24-cva0ls0-509,u24-cva0000-129
#SBATCH --output=logs/features_%j.out
set -euo pipefail
export TINCEPTION_ENV=cluster
cd "$HOME/Tinception"
.venv/bin/python -W ignore scripts/06_features.py
