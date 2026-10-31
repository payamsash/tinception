#!/bin/bash -l
# H3: brain-clinical dimensions (PLS) and subtypes in patients (scripts/09_dimensions.py). Submit from the repo root:
#   sbatch slurm/09_dimensions.sbatch
#SBATCH --job-name=h3
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --exclude=u24-cva0ls0-509,u24-cva0000-129,u24-chaiam0-613   # broken apptainer / share mount on these nodes
#SBATCH --output=logs/h3_%j.out
set -euo pipefail
export TINCEPTION_ENV=cluster OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
cd "$HOME/Tinception"
.venv/bin/python -W ignore scripts/09_dimensions.py "$@"
