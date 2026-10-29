#!/bin/bash -l
# Assemble and calibrate the phase-2 normative maps (scripts/07b_normative_maps.py collect).
# Submit from the repo root: sbatch slurm/07b_collect.sbatch
#SBATCH --job-name=nmcollect
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --exclude=u24-cva0ls0-509,u24-cva0000-129,u24-chaiam0-613   # broken apptainer on these nodes (user ns / resolv.conf mount), 2026-09
#SBATCH --output=logs/nmcollect_%j.out
set -euo pipefail
export TINCEPTION_ENV=cluster
cd "$HOME/Tinception"
.venv/bin/python -W ignore scripts/07b_normative_maps.py collect
