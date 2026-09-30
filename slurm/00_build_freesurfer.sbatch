#!/bin/bash -l
# Build the exact-match FreeSurfer 8.0.0 Ubuntu20 image. Submit from repo root.
#SBATCH --job-name=build_fs
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=logs/build_fs_%j.out
set -euo pipefail
module load apptainer/1.5.0
export APPTAINER_CACHEDIR=/scratch/$USER/apptainer_cache
export APPTAINER_TMPDIR=/scratch/$USER/apptainer_tmp
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"
OUT=/scratch/$USER/freesurfer-8.0.0-ubuntu20.sif   # build on scratch (fast), then copy to the share
apptainer build --fakeroot --ignore-fakeroot-command --mksquashfs-args "-processors $SLURM_CPUS_PER_TASK" "$OUT" containers/freesurfer-8.0.0-ubuntu20.def
apptainer exec "$OUT" bash -c 'cat $FREESURFER_HOME/build-stamp.txt; command -v recon-all segment_subregions'
cp "$OUT" /shares/tinnitus.orl.med.uzh/containers/ && sha256sum /shares/tinnitus.orl.med.uzh/containers/freesurfer-8.0.0-ubuntu20.sif
