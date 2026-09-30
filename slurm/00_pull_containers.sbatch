#!/bin/bash -l
# Pull the pinned container images (containers/images.txt) into the project share.
# Submit from the repo root:  sbatch slurm/00_pull_containers.sbatch
#SBATCH --job-name=pull_sif
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/pull_%j.out
set -euo pipefail
module load apptainer/1.5.0
DEST=/shares/tinnitus.orl.med.uzh/containers
export APPTAINER_CACHEDIR=/scratch/$USER/apptainer_cache
export APPTAINER_TMPDIR=/scratch/$USER/apptainer_tmp
mkdir -p "$DEST" "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"
while read -r name uri; do
  [[ -z "$name" || "$name" == \#* ]] && continue
  [[ -f "$DEST/$name" ]] && { echo "exists: $name"; continue; }
  echo "pull $name <- $uri"; apptainer pull "$DEST/$name" "$uri"
done < containers/images.txt
sha256sum "$DEST"/*.sif | tee containers/digests_cluster.txt
apptainer exec "$DEST/freesurfer-8.0.0.sif" cat /usr/local/freesurfer/build-stamp.txt
