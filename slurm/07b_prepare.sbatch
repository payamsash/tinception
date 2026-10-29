#!/bin/bash -l
# Phase-2 normative inputs: fsaverage5 vertex matrices (cortex only) and CAT12 grey matter resampled to 3 mm.
# Submit from the repo root: sbatch slurm/07b_prepare.sbatch [set ...]
#SBATCH --job-name=nmprep
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --exclude=u24-cva0ls0-509,u24-cva0000-129,u24-chaiam0-613   # broken apptainer on these nodes (user ns / resolv.conf mount), 2026-09
#SBATCH --output=logs/nmprep_%j.out
set -euo pipefail
export TINCEPTION_ENV=cluster
cd "$HOME/Tinception"
.venv/bin/python -W ignore scripts/07b_normative_maps.py prepare ${1:+--sets "$@"}
