#!/bin/bash -l
# LOSO transfer check of the ROI normative models, one task = reference x model x feature set x site.
# Submit from the repo root:
#   T=/shares/tinnitus.orl.med.uzh/tinception/derivatives/normative/loso/tasks.tsv
#   sbatch --array=1-$(($(wc -l < $T)-1))%100 slurm/07c_normative_loso.sbatch $T
#SBATCH --job-name=nmloso
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --exclude=u24-cva0ls0-509,u24-cva0000-129,u24-chaiam0-613   # broken apptainer on these nodes (user ns / resolv.conf mount), 2026-09
#SBATCH --output=logs/nmloso_%A_%a.out
set -euo pipefail
export TINCEPTION_ENV=cluster OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
cd "$HOME/Tinception"
.venv/bin/python -W ignore scripts/07c_normative_loso.py fit --task-file "$1" --task "$SLURM_ARRAY_TASK_ID"
