#!/bin/bash -l
# Network context (plan 6.9 B): hub / epicentre models of group maps and individual deviations
# (scripts/10_network.py). Submit from the repo root:  sbatch slurm/10_network.sbatch
#SBATCH --job-name=network
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --exclude=u24-cva0ls0-509,u24-cva0000-129,u24-chaiam0-613   # broken apptainer / share mount on these nodes
#SBATCH --output=logs/network_%j.out
set -euo pipefail
export TINCEPTION_ENV=cluster OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
cd "$HOME/Tinception"
.venv/bin/python -W ignore scripts/10_network.py "$@"
