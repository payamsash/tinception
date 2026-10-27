#!/bin/bash -l
# Extract all FreeSurfer feature tables + fsaverage5 vertex matrices -> derivatives/features/
# Submit from repo root (after slurm/06_annot_stats.sbatch): sbatch slurm/06_features.sbatch
# FS_ROOT=<SUBJECTS_DIR> (with --export=ALL) reads another FreeSurfer tree than config paths.freesurfer;
# FEATURES_OUT=<dir> writes somewhere else than <derivatives>/features (set it for any non-primary tree)
#SBATCH --job-name=features
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --exclude=u24-cva0ls0-509,u24-cva0000-129
#SBATCH --output=logs/features_%j.out
set -euo pipefail
export TINCEPTION_ENV=cluster
cd "$HOME/Tinception"
.venv/bin/python -W ignore scripts/06_features.py ${FS_ROOT:+--fs-root "$FS_ROOT"} ${FEATURES_OUT:+--out "$FEATURES_OUT"}
