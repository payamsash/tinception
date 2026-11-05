#!/bin/bash -l
# Pretrained brain age (plan 6.11 A): both pyment brain-age models on the FastSurfer-preprocessed images
# (author's container, CPU). Run after slurm/19_brainage_preprocess.sbatch:
#   sbatch --dependency=afterany:<preprocess job> slurm/19_brainage_predict.sbatch
#SBATCH --job-name=pyment_pred
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --exclude=u24-cva0ls0-509,u24-cva0000-129,u24-chaiam0-613   # broken apptainer / share mount on these nodes
#SBATCH --output=logs/pyment_pred_%j.out
set -euo pipefail
module load apptainer/1.5.0
BASE=/shares/tinnitus.orl.med.uzh
SIF=$BASE/containers/pyment-predict-5.1.0.sif
D=$BASE/tinception/derivatives/pyment
run() { apptainer exec --cleanenv -B "$D/fastsurfer:/fastsurfer,$D:/output" "$SIF" pyment-predict /fastsurfer "$@"; }
# multi-task model (default architecture sfcn-multi; all heads written, 'age' is used)
run -w multi-2025 -d /output/predictions_multi-2025.csv
# pure brain-age model needs its own architecture
run -w reg-2025 -m sfcn-reg -t age -d /output/predictions_reg-2025.csv
wc -l "$D"/predictions_*.csv
