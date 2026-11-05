#!/bin/bash -l
# Pretrained brain age (plan 6.11 A): pyment preprocessing (author's container, FastSurfer conform + brain mask)
# for one subject per task. List = rawdata/code/cat12_list.tsv (participant_id, t1: same inputs as FreeSurfer/CAT12).
#   L=/shares/tinnitus.orl.med.uzh/tinception/rawdata/code/cat12_list.tsv
#   sbatch --array=1-$(($(wc -l < $L)-1))%100 slurm/19_brainage_preprocess.sbatch $L
#SBATCH --job-name=pyment_pre
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --exclude=u24-cva0ls0-509,u24-cva0000-129,u24-chaiam0-613   # broken apptainer / share mount on these nodes
#SBATCH --output=logs/pyment_pre_%A_%a.out
set -euo pipefail
module load apptainer/1.5.0
LIST=$1
BASE=/shares/tinnitus.orl.med.uzh
SIF=$BASE/containers/pyment-preprocess-5.0.0.sif
OUT=$BASE/tinception/derivatives/pyment
ROW=$((SLURM_ARRAY_TASK_ID + 1))
SUB=$(awk -F'\t' -v i=$ROW 'NR==i {print $1}' "$LIST")
T1=$(awk -F'\t' -v i=$ROW 'NR==i {print $2}' "$LIST")
MASK="$OUT/fastsurfer/${SUB}_T1w/mri/mask.mgz"   # exact path (a substring match let sub-uiuc1 hit sub-uiuc108)
# set +e before exit: a login shell runs ~/.bash_logout on `exit`, and under set -e a false test there gives exit 1
if [[ -f $MASK ]]; then echo "done: $SUB"; set +e; exit 0; fi
WORK=/scratch/$USER/pyment/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$WORK/input" "$WORK/output" "$WORK/licenses" "$OUT"
trap 'rm -rf "$WORK"' EXIT
cp "$T1" "$WORK/input/${SUB}_T1w.nii.gz"
cp "$HOME/.freesurfer/license.txt" "$WORK/licenses/freesurfer.txt"
apptainer run --cleanenv -B "$WORK/input:/input,$WORK/output:/output,$WORK/licenses:/licenses" "$SIF"
mkdir -p "$OUT/fastsurfer"
[[ -f "$WORK/output/fastsurfer/${SUB}_T1w/mri/mask.mgz" ]] || { echo "no mask.mgz for $SUB"; ls -R "$WORK/output" | head -40; exit 1; }
# copy everything except the mask first, the mask last: a mask on the share then means a complete copy
rsync -a --exclude "fastsurfer/${SUB}_T1w/mri/mask.mgz" "$WORK/output/" "$OUT/"
rsync -a "$WORK/output/fastsurfer/${SUB}_T1w/mri/mask.mgz" "$OUT/fastsurfer/${SUB}_T1w/mri/"
echo "done: $SUB"
