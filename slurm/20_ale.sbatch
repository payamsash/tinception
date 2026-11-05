#!/bin/bash -l
# Coordinate-based meta-analysis (plan 6.10 B2; scripts/20_ale.py) on data/meta/coordinates_final.tsv:
#   sbatch slurm/20_ale.sbatch --descriptive          (all conditions; < 17 experiments = descriptive only)
#SBATCH --job-name=ale
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --exclude=u24-cva0ls0-509,u24-cva0000-129,u24-chaiam0-613,u24-chaihm0-636   # broken apptainer / share mount on these nodes
#SBATCH --output=logs/ale_%j.out
set -euo pipefail
export TINCEPTION_ENV=cluster OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$HOME/Tinception"
.venv/bin/python -W ignore scripts/20_ale.py run "$@"
