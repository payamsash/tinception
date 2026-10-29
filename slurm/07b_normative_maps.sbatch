#!/bin/bash -l
# Phase-2 normative models, one task = reference x model x set x chunk (scripts/07b_normative_maps.py tasks).
# Submit from the repo root:
#   T=/shares/tinnitus.orl.med.uzh/tinception/derivatives/normative/maps/tasks.tsv
#   sbatch --array=1-$(($(wc -l < $T)-1))%100 slurm/07b_normative_maps.sbatch $T
#SBATCH --job-name=nmmaps
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --exclude=u24-cva0ls0-509,u24-cva0000-129,u24-chaiam0-613   # broken apptainer on these nodes (user ns / resolv.conf mount), 2026-09
#SBATCH --output=logs/nmmaps_%A_%a.out
set -euo pipefail
export TINCEPTION_ENV=cluster OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
cd "$HOME/Tinception"
.venv/bin/python -W ignore scripts/07b_normative_maps.py fit --task-file "$1" --task "$SLURM_ARRAY_TASK_ID"
