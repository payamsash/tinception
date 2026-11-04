#!/bin/bash -l
# Vertex / VBM group maps with PALM (plan 6.6; scripts/13_maps.py). Prepare first, then one task per analysis:
#   .venv/bin/python scripts/13_maps.py prepare            (run via srun: needs ~40 GB for the VBM stacks)
#   sbatch --array=1-$(($(wc -l < <derivatives>/maps_palm/tasks.tsv)-1)) slurm/13_maps.sbatch
#SBATCH --job-name=palm
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --exclude=u24-cva0ls0-509,u24-cva0000-129,u24-chaiam0-613   # broken apptainer / share mount on these nodes
#SBATCH --output=logs/palm_%A_%a.out
set -euo pipefail
export TINCEPTION_ENV=cluster
cd "$HOME/Tinception"
.venv/bin/python -W ignore scripts/13_maps.py fit --task "$SLURM_ARRAY_TASK_ID" "$@"
