#!/bin/bash -l
# Gradient contextualization (plan 6.9 C; scripts/11_context.py run). The axes table data/context/dk_gradients.tsv
# is built locally first (scripts/11_context.py prepare; needs internet), committed, and synced to ~/Tinception/data.
# Submit from the repo root:
#   sbatch slurm/11_context.sbatch
#SBATCH --job-name=context
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --exclude=u24-cva0ls0-509,u24-cva0000-129,u24-chaiam0-613   # broken apptainer / share mount on these nodes
#SBATCH --output=logs/context_%j.out
set -euo pipefail
export TINCEPTION_ENV=cluster OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
cd "$HOME/Tinception"
.venv/bin/python -W ignore scripts/11_context.py run "$@"
