#!/bin/bash -l
# Step 08 group models: one subcommand of scripts/08_group.py (h1, h1d-ci, ...). Submit from the repo root:
#   sbatch slurm/08_group.sbatch h1
#SBATCH --job-name=group
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --exclude=u24-cva0ls0-509,u24-cva0000-129,u24-chaiam0-613   # broken apptainer on these nodes (user ns / resolv.conf mount), 2026-09
#SBATCH --output=logs/group_%j.out
set -euo pipefail
export TINCEPTION_ENV=cluster OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
cd "$HOME/Tinception"
.venv/bin/python -W ignore scripts/08_group.py "$@"
