#!/bin/bash -l
# Pilot: N4 bias-field correction of raw T1w (64-ch coil, no prescan normalise) before recon-all.
# Usage: sbatch slurm/pilot_n4.sbatch <list.tsv> <outdir>
#SBATCH --job-name=n4
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --exclude=u24-cva0ls0-509,u24-cva0000-129
#SBATCH --output=logs/n4_%j.out
set -euo pipefail
module load apptainer/1.5.0
LIST=$1; OUT=$2; mkdir -p "$OUT"
RAW=/shares/tinnitus.orl.med.uzh/tinception/rawdata
SIF=/shares/tinnitus.orl.med.uzh/containers/freesurfer-8.0.0-ubuntu20.sif
tail -n +2 "$LIST" | cut -f1 | while read -r SUB; do
  [[ -f $OUT/${SUB}_T1w.nii.gz ]] && continue   # done earlier (outputs are written atomically)
  apptainer exec -B "$RAW,$OUT" "$SIF" bash -c "source \$FREESURFER_HOME/SetUpFreeSurfer.sh >/dev/null; \
    AntsN4BiasFieldCorrectionFs -i $RAW/$SUB/anat/${SUB}_T1w.nii.gz -o $OUT/${SUB}_tmp_T1w.nii.gz" </dev/null
  mv "$OUT/${SUB}_tmp_T1w.nii.gz" "$OUT/${SUB}_T1w.nii.gz"
  echo "n4 done $SUB"
done
