#!/bin/bash -l
# Frozen pyment embeddings (plan 6.11 B; scripts/19b_pyment_embed.py) inside the pinned predict container, from the
# crops written by slurm/19_brainage_predict.sbatch:
#   sbatch --dependency=afterok:<predict job> slurm/19_brainage_embed.sbatch
#SBATCH --job-name=pyment_emb
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --exclude=u24-cva0ls0-509,u24-cva0000-129,u24-chaiam0-613,u24-chaihm0-636   # broken apptainer / share mount on these nodes
#SBATCH --output=logs/pyment_emb_%j.out
set -euo pipefail
module load apptainer/1.5.0
BASE=/shares/tinnitus.orl.med.uzh
SIF=$BASE/containers/pyment-predict-5.1.0.sif
D=$BASE/tinception/derivatives/pyment
apptainer exec --cleanenv -B "$D/fastsurfer:/fastsurfer:ro,$D:/output,$HOME/Tinception/scripts:/code:ro" "$SIF" \
  python /code/19b_pyment_embed.py /fastsurfer /output
wc -l "$D"/embeddings_*.tsv
