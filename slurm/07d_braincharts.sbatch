#!/bin/bash -l
# External normative reference (plan 6.8): transfer of pretrained braincharts models with 2-fold cross-fitted site
# adaptation (scripts/braincharts_transfer.py, PCNtoolkit 0.35 env). Run `07d_braincharts.py prepare` first.
#   sbatch slurm/07d_braincharts.sbatch          then  .venv/bin/python scripts/07d_braincharts.py analyze
#SBATCH --job-name=braincharts
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --exclude=u24-cva0ls0-509,u24-cva0000-129,u24-chaiam0-613   # broken apptainer / share mount on these nodes
#SBATCH --output=logs/braincharts_%j.out
set -euo pipefail
BASE=/shares/tinnitus.orl.med.uzh
export TMPDIR=/scratch/$USER/braincharts_tmp && mkdir -p "$TMPDIR"
cd "$HOME/Tinception"
$BASE/containers/envs/pcn035/bin/python -W ignore scripts/braincharts_transfer.py \
  --data $BASE/tinception/derivatives/normative_braincharts --braincharts $BASE/containers/braincharts "$@"
