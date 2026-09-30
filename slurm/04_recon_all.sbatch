#!/bin/bash -l
# FreeSurfer 8.0.0 recon-all (+ qcache) for subjects listed in a TSV (default: rawdata/code/needs_freesurfer.tsv).
# Submit:  N=$(($(wc -l < $LIST) - 1)); sbatch --array=1-${N}%50 slurm/04_recon_all.sbatch $LIST
# Andes T1-only sensitivity run (primary andes outputs keep their -T2pial surfaces):
#   SUBJECTS_DIR_OVERRIDE=/shares/tinnitus.orl.med.uzh/tinception/derivatives/freesurfer-8.0.0-andesT1only SUBSEG=0 \
#   sbatch --export=ALL --array=1-106%50 slurm/04_recon_all.sbatch rawdata/code/andes_t1only.tsv
#SBATCH --job-name=recon
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G   # 24G is OOM-killed during SynthSeg in FS 8.0.0
#SBATCH --exclude=u24-cva0ls0-509,u24-cva0000-129   # broken apptainer on these nodes (user ns / resolv.conf mount), 2026-09
#SBATCH --output=logs/recon_%A_%a.out
set -euo pipefail

LIST=${1:-/shares/tinnitus.orl.med.uzh/tinception/rawdata/code/needs_freesurfer.tsv}
RAW=/shares/tinnitus.orl.med.uzh/tinception/rawdata
SUBJECTS_DIR=${SUBJECTS_DIR_OVERRIDE:-/shares/tinnitus.orl.med.uzh/tinception/derivatives/freesurfer-8.0.0}
SIF=/shares/tinnitus.orl.med.uzh/containers/freesurfer-8.0.0-ubuntu20.sif  # exact legacy build
MCR=/shares/tinnitus.orl.med.uzh/containers/MCR_R2019b/v97                   # for segmentHA_T1.sh etc.
LICENSE=$HOME/.freesurfer/license.txt

SUB=$(awk -F'\t' -v i=$((SLURM_ARRAY_TASK_ID + 1)) 'NR==i {print $1}' "$LIST")
T1=$RAW/$SUB/anat/${SUB}_T1w.nii.gz
[[ -f "$T1" || -f "$SUBJECTS_DIR/$SUB/scripts/recon-all.done" ]] || { echo "missing $T1"; exit 1; }
mkdir -p "$SUBJECTS_DIR"

module load apptainer/1.5.0
# every call sources FreeSurfer's own setup script (FS 8 recon-all needs $FREESURFER etc.)
run() { apptainer exec -B "$RAW,$SUBJECTS_DIR,$LICENSE,$MCR:/usr/local/freesurfer/8.0.0/MCRv97" \
          --env SUBJECTS_DIR="$SUBJECTS_DIR",FS_LICENSE="$LICENSE" "$SIF" \
          bash -c 'export SUBJECTS_DIR_KEEP=$SUBJECTS_DIR; source $FREESURFER_HOME/SetUpFreeSurfer.sh >/dev/null;
                   export SUBJECTS_DIR=$SUBJECTS_DIR_KEEP; exec "$@"' _ "$@"; }
run bash -c 'grep -q "ubuntu20_x86_64-8.0.0-20250204-0f2bd3a" $FREESURFER_HOME/build-stamp.txt' \
  || { echo "wrong FreeSurfer build"; exit 2; }

# New subjects are run T1-only (no -T2pial); andes legacy outputs used -T2pial (see sensitivity run above)
D=$SUBJECTS_DIR/$SUB/scripts
# FS writes recon-all.done even on failure: done == .done present AND no .error
if [[ -f $D/recon-all.done && ! -f $D/recon-all.error ]]; then
  [[ -f $SUBJECTS_DIR/$SUB/surf/lh.thickness.fwhm10.fsaverage.mgh ]] || \
    run recon-all -s "$SUB" -qcache -threads "$SLURM_CPUS_PER_TASK" -no-isrunning   # legacy: qcache only
elif [[ -f $SUBJECTS_DIR/$SUB/mri/orig/001.mgz ]]; then
  rm -f $D/recon-all.error $D/IsRunning.*                                             # resume failed run
  run recon-all -s "$SUB" -all -qcache -threads "$SLURM_CPUS_PER_TASK" -no-isrunning
else
  run recon-all -s "$SUB" -i "$T1" -all -qcache -threads "$SLURM_CPUS_PER_TASK"
fi

[[ "${SUBSEG:-1}" == "1" ]] || exit 0
# T1-only subsegmentations with the SAME tools as the 701 legacy subjects (MATLAB-runtime versions,
# outputs *-T1.v22 / ThalamicNuclei.v13.T1), only where missing. Brainstem/hypothalamus not used.
MRI=$SUBJECTS_DIR/$SUB/mri
if [[ ! -f $MRI/lh.hippoSfVolumes-T1.v22.txt || ! -f $MRI/rh.amygNucVolumes-T1.v22.txt ]]; then
  run segmentHA_T1.sh "$SUB" "$SUBJECTS_DIR"
fi
if [[ ! -f $MRI/ThalamicNuclei.v13.T1.volumes.txt ]]; then
  run segmentThalamicNuclei.sh "$SUB" "$SUBJECTS_DIR"
fi
