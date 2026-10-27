#!/bin/bash -l
#SBATCH --job-name=asegchk
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --exclude=u24-cva0ls0-509,u24-cva0000-129
#SBATCH --output=logs/asegchk_%j.out
cd "$HOME/Tinception"
.venv/bin/python -W ignore scripts/check_aseg_consistency.py /shares/tinnitus.orl.med.uzh/tinception/derivatives/freesurfer-8.0.0 /shares/tinnitus.orl.med.uzh/tinception/derivatives/features/aseg_consistency.tsv
