#!/bin/bash -l
# CAT12.9 r2665 volume pipeline: segment (config/cat12_segment.m -> mwp1/mwp2, TIV, IQR, ROI volumes),
# then smooth mwp1 at 6 and 8 mm (FWHM env overrides). One subject per task; finished subjects are skipped.
# List: TSV with columns participant_id, t1 (scripts/05_cat12.py list). Submit from the repo root:
#   L=/shares/tinnitus.orl.med.uzh/tinception/rawdata/code/cat12_list.tsv
#   sbatch --array=1-$(($(wc -l < $L)-1))%100 slurm/05_cat12.sbatch $L
#SBATCH --job-name=cat12
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --exclude=u24-cva0ls0-509,u24-cva0000-129   # broken apptainer on these nodes (user ns / resolv.conf mount), 2026-09
#SBATCH --output=logs/cat12_%A_%a.out
set -euo pipefail
module load apptainer/1.5.0

LIST=$1
BASE=/shares/tinnitus.orl.med.uzh/tinception
OUT=${CAT12_OUT:-$BASE/derivatives/cat12-12.9}   # CAT12_OUT=<dir> for pilots
SIF=/shares/tinnitus.orl.med.uzh/containers/cat12-12.9-r2665.sif
BATCH=$SLURM_SUBMIT_DIR/config/cat12_segment.m
SMOOTH=/opt/spm/standalone/cat_standalone_smooth.m   # inside the image
FWHM=${FWHM:-6 8}                                    # = analysis.yaml cat12.smoothing_fwhm

ROW=$((SLURM_ARRAY_TASK_ID + 1))
SUB=$(awk -F'\t' -v i=$ROW 'NR==i {print $1}' "$LIST")
T1=$(awk -F'\t' -v i=$ROW 'NR==i {print $2}' "$LIST")
D=$OUT/$SUB
MWP=$D/mri/mwp1${SUB}_T1w.nii

is_done() {
  [[ -f $D/report/cat_${SUB}_T1w.xml && -f $MWP ]] || return 1
  for k in $FWHM; do [[ -f $D/mri/s${k}mwp1${SUB}_T1w.nii ]] || return 1; done
}
if is_done; then echo "done: $SUB"; exit 0; fi
[[ -f $T1 ]] || { echo "missing $T1"; exit 1; }

WORK=/scratch/$USER/tinception_work/cat12/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$D" "$WORK"
export MCR_CACHE_ROOT=$WORK/mcr   # per-task MATLAB Runtime cache: a shared one in $HOME breaks with parallel tasks
run() { apptainer exec -B "$BASE,$WORK,$(dirname "$BATCH")" "$SIF" "$@"; }

trap 'rm -f "$D/${SUB}_T1w.nii"; rm -rf "$WORK"' EXIT   # input copy + MCR cache, also on failure
if [[ ! -f $MWP ]]; then
  # CAT12 writes mri/ report/ label/ next to its input -> uncompressed copy inside the subject dir
  gunzip -c "$T1" > "$D/${SUB}_T1w.nii"
  run cat_standalone.sh -b "$BATCH" "$D/${SUB}_T1w.nii" || echo "cat_standalone.sh exited with $?"
fi
[[ -f $MWP ]] || { echo "CAT12 segmentation failed: $SUB (see $D/err or the log above)"; exit 1; }

for k in $FWHM; do
  [[ -f $D/mri/s${k}mwp1${SUB}_T1w.nii ]] || \
    run cat_standalone.sh -b "$SMOOTH" "$MWP" -a1 "[$k $k $k]" -a2 "'s$k'"
done
is_done || { echo "missing outputs: $SUB"; exit 1; }
echo "done: $SUB"
