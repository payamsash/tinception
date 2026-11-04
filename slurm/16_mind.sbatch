#!/bin/bash -l
# MIND networks (plan 6.12; scripts/16_mind.py compute): one subject per task, all parcellations (~35 min).
#   .venv/bin/python scripts/16_mind.py subjects; L=<rawdata>/code/mind_subjects.tsv
#   sbatch --array=1-$(($(wc -l < $L)-1))%200 slurm/16_mind.sbatch   then  .venv/bin/python scripts/16_mind.py collect
#SBATCH --job-name=mind
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --exclude=u24-cva0ls0-509,u24-cva0000-129,u24-chaiam0-613   # broken apptainer / share mount on these nodes
#SBATCH --output=logs/mind_%A_%a.out
set -euo pipefail
export TINCEPTION_ENV=cluster OMP_NUM_THREADS=1
cd "$HOME/Tinception"
.venv/bin/python -W ignore scripts/16_mind.py compute --task "$SLURM_ARRAY_TASK_ID"
