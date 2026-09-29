#!/usr/bin/env bash
# Step 00 on the UZH ScienceCluster: check the environment, pull pinned containers, sync the env.
# Usage (on the cluster login node):  bash scripts/00_setup_env.sh [--check-only]
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONTAINERS="${TINCEPTION_CONTAINERS:-/shares/tinnitus.orl.med.uzh/containers}"

echo "== host: $(hostname)  user: $USER"
echo "== storage"
df -h "$HOME" "/scratch/$USER" 2>/dev/null || true
ls -d /shares/* 2>/dev/null | head || true

echo "== scheduler / container runtime"
command -v sbatch && sinfo -s | head
module avail 2>&1 | grep -iE "apptainer|singularity|mamba|anaconda" | head || true
command -v apptainer || command -v singularity || echo "!! load apptainer/singularity module"

echo "== FreeSurfer licence"
[[ -f "$HOME/.freesurfer/license.txt" ]] && echo ok || echo "!! put license.txt at ~/.freesurfer/license.txt"

[[ "${1:-}" == "--check-only" ]] && exit 0

echo "== containers -> $CONTAINERS"
mkdir -p "$CONTAINERS"
RUNTIME=$(command -v apptainer || command -v singularity)
pull() { [[ -f "$CONTAINERS/$1" ]] || "$RUNTIME" pull "$CONTAINERS/$1" "$2"; }
while read -r name uri; do
  [[ -z "$name" || "$name" == \#* ]] && continue
  pull "$name" "$uri"
done < "$REPO/containers/images.txt"
sha256sum "$CONTAINERS"/*.sif > "$REPO/containers/digests_$(hostname -s).txt"

echo "== python env (uv)"
command -v uv || curl -LsSf https://astral.sh/uv/install.sh | sh
cd "$REPO" && uv sync --all-extras
