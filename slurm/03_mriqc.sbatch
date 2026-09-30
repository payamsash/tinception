#!/bin/bash -l
# MRIQC 24.0.2 (T1w, participant level) for subjects listed in a TSV (first column = participant_id).
# Submit from repo root:
#   L=/shares/tinnitus.orl.med.uzh/tinception/rawdata/code/mriqc_list.tsv
#   sbatch --array=1-$(($(wc -l < $L)-1))%60 slurm/03_mriqc.sbatch $L
# The first task pre-fetches TemplateFlow; run `sbatch --array=1 ...` once first to fill the cache.
#SBATCH --job-name=mriqc
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --exclude=u24-cva0ls0-509,u24-cva0000-129   # broken apptainer on these nodes (user ns / resolv.conf mount), 2026-09
#SBATCH --output=logs/mriqc_%A_%a.out
set -euo pipefail
module load apptainer/1.5.0

LIST=$1
BASE=/shares/tinnitus.orl.med.uzh/tinception
RAW=$BASE/rawdata
OUT=$BASE/derivatives/mriqc-24.0.2
SIF=/shares/tinnitus.orl.med.uzh/containers/mriqc-24.0.2.sif
TF=/shares/tinnitus.orl.med.uzh/containers/templateflow
WORK=/scratch/$USER/tinception_work/mriqc/$SLURM_ARRAY_JOB_ID/$SLURM_ARRAY_TASK_ID
mkdir -p "$OUT" "$TF" "$WORK"

SUB=$(awk -F'\t' -v i=$((SLURM_ARRAY_TASK_ID + 1)) 'NR==i {print $1}' "$LIST")
LABEL=${SUB#sub-}
if [[ -f $OUT/${SUB}/anat/${SUB}_T1w.json || -f $OUT/${SUB}_T1w.json ]]; then echo "done: $SUB"; exit 0; fi

# Only this subject's T1w is indexed (bids filter) -> fast, independent tasks.
cat > "$WORK/filter.json" <<EOF
{"t1w": {"datatype": "anat", "suffix": "T1w", "acquisition": null}}
EOF

apptainer run --cleanenv -B "$RAW:/data:ro,$OUT:/out,$WORK:/work,$TF:/templateflow" \
  --env TEMPLATEFLOW_HOME=/templateflow "$SIF" \
  /data /out participant --participant-label "$LABEL" -m T1w \
  --bids-filter-file /work/filter.json --no-sub --nprocs "$SLURM_CPUS_PER_TASK" --mem 14 \
  -w /work --verbose-reports --notrack

rm -rf "$WORK"
