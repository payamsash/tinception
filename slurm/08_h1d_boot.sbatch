#!/bin/bash -l
# H1d bootstrap of the attenuation (5,000 resamples within site per region, split in chunks).
# Submit from the repo root after `08_group.py h1` and `08_group.py h1d-tasks`:
#   T=/shares/tinnitus.orl.med.uzh/tinception/results/08_group/h1/h1d_tasks.tsv
#   sbatch --array=1-$(($(wc -l < $T)-1))%100 slurm/08_h1d_boot.sbatch
#SBATCH --job-name=h1dboot
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --exclude=u24-cva0ls0-509,u24-cva0000-129,u24-chaiam0-613   # broken apptainer on these nodes (user ns / resolv.conf mount), 2026-09
#SBATCH --output=logs/h1dboot_%A_%a.out
set -euo pipefail
export TINCEPTION_ENV=cluster OMP_NUM_THREADS=1
cd "$HOME/Tinception"
.venv/bin/python -W ignore scripts/08_group.py h1d-boot --task "$SLURM_ARRAY_TASK_ID"
