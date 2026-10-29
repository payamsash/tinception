#!/bin/bash -l
# Map fsaverage parcellations (Schaefer-400 7-network, HCP-MMP1) onto each subject and compute
# per-region cortical stats (thickness, area, GM volume) with mris_anatomical_stats.
# Submit: L=<list.tsv, first column participant_id>; sbatch --array=1-N%100 slurm/06_annot_stats.sbatch $L
#SBATCH --job-name=annot
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --exclude=u24-cva0ls0-509,u24-cva0000-129,u24-chaiam0-613   # broken apptainer on these nodes (user ns / resolv.conf mount), 2026-09
#SBATCH --output=logs/annot_%A_%a.out
set -euo pipefail
module load apptainer/1.5.0

LIST=$1
SUBJECTS_DIR=${SUBJECTS_DIR_OVERRIDE:-/shares/tinnitus.orl.med.uzh/tinception/derivatives/freesurfer-8.0.0}
SIF=/shares/tinnitus.orl.med.uzh/containers/freesurfer-8.0.0-ubuntu20.sif
LICENSE=$HOME/.freesurfer/license.txt
SUB=$(awk -F'\t' -v i=$((SLURM_ARRAY_TASK_ID + 1)) 'NR==i {print $1}' "$LIST")

run() { apptainer exec -B "$SUBJECTS_DIR,$LICENSE" --env SUBJECTS_DIR="$SUBJECTS_DIR",FS_LICENSE="$LICENSE" "$SIF" \
          bash -c 'export K=$SUBJECTS_DIR; source $FREESURFER_HOME/SetUpFreeSurfer.sh >/dev/null; export SUBJECTS_DIR=$K; exec "$@"' _ "$@"; }

for ANNOT in Schaefer2018_400Parcels_7Networks_order HCPMMP1; do
  for H in lh rh; do
    A=$SUBJECTS_DIR/$SUB/label/$H.$ANNOT.annot
    S=$SUBJECTS_DIR/$SUB/stats/$H.$ANNOT.stats
    [[ -f $A ]] || run mri_surf2surf --srcsubject fsaverage --trgsubject "$SUB" --hemi $H \
                      --sval-annot "$SUBJECTS_DIR/fsaverage/label/$H.$ANNOT.annot" --tval "$A"
    # same options as recon-all's aparc stats (-th3 volumes, cortex mask) -> comparable to DK/DKT/Destrieux
    [[ -f $S ]] || run mris_anatomical_stats -th3 -mgz -cortex "$SUBJECTS_DIR/$SUB/label/$H.cortex.label" \
                      -a "$A" -f "$S" -b "$SUB" $H white
  done
done
echo "done $SUB"
