#!/bin/bash -l
#SBATCH --job-name=asegchk
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --exclude=u24-cva0ls0-509,u24-cva0000-129
#SBATCH --output=logs/asegchk_%j.out
# Usage: sbatch slurm/check_aseg.sbatch [<SUBJECTS_DIR> [<out.tsv>]]
# default output is named after the tree: features/aseg_consistency_<basename of SUBJECTS_DIR>.tsv
set -euo pipefail
D=/shares/tinnitus.orl.med.uzh/tinception/derivatives
FS=${1:-$D/freesurfer-8.0.0}
cd "$HOME/Tinception"
.venv/bin/python -W ignore scripts/check_aseg_consistency.py "$FS" "${2:-$D/features/aseg_consistency_$(basename "$FS").tsv}"
