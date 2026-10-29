#!/bin/bash -l
# Normative models (BLR), one task = reference x model x feature set (scripts/07_normative.py tasks).
# Submit from the repo root:
#   T=/shares/tinnitus.orl.med.uzh/tinception/derivatives/normative/tasks.tsv
#   sbatch --array=1-$(($(wc -l < $T)-1)) slurm/07_normative.sbatch $T
#SBATCH --job-name=normative
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --exclude=u24-cva0ls0-509,u24-cva0000-129,u24-chaiam0-613   # broken apptainer on these nodes (user ns / resolv.conf mount), 2026-09
#SBATCH --output=logs/normative_%A_%a.out
set -euo pipefail
export TINCEPTION_ENV=cluster OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
cd "$HOME/Tinception"
.venv/bin/python -W ignore scripts/07_normative.py fit --task-file "$1" --task "$SLURM_ARRAY_TASK_ID"
