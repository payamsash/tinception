#!/bin/bash -l
# Leave-site-out prediction (plan 6.8; scripts/14_prediction.py). One array task per prediction task:
#   sbatch --array=1-6 slurm/14_prediction.sbatch      then  .venv/bin/python scripts/14_prediction.py collect
#SBATCH --job-name=predict
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --exclude=u24-cva0ls0-509,u24-cva0000-129,u24-chaiam0-613   # broken apptainer / share mount on these nodes
#SBATCH --output=logs/predict_%A_%a.out
set -euo pipefail
export TINCEPTION_ENV=cluster OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$HOME/Tinception"
# FEATURES=pyment-multi-2025 (or pyment-reg-2025) for the frozen-embedding probe (plan 6.11 B; --array=1-6)
.venv/bin/python -W ignore scripts/14_prediction.py --features "${FEATURES:-roi}" run "$@"
