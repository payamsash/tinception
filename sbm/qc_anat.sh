#!/bin/bash

## set Paths
export FREESURFER_HOME=/usr/local/freesurfer/8.0.0
export SUBJECTS_DIR=/home/ubuntu/volume/subjects_fs_dir
export LD_LIBRARY_PATH=$FREESURFER_HOME/MCRv97/runtime/glnxa64:$FREESURFER_HOME/MCRv97/bin/glnxa64:$FREESURFER_HOME/MCRv97/sys/os/glnxa64:$FREESURFER_HOME/MCRv97/extern/bin/glnxa64
source $FREESURFER_HOME/SetUpFreeSurfer.sh
export PATH="/home/ubuntu/ants-2.6.2/ants-2.6.2/bin:$PATH"

# Create output directory if it doesn't exist
OUTPUT_DIR="/home/ubuntu/volume/QC/bids_data"
mkdir -p "$OUTPUT_DIR"

# Loop through all directories in SUBJECTS_DIR
for subj in "$SUBJECTS_DIR"/*; do
    subj_name=$(basename "$subj")

    if [[ -d "$subj" && -d "$subj/hist" ]]; then
        original_file="$subj/mri/orig/001.mgz"

        anat_dir="$OUTPUT_DIR/sub-${subj_name}/anat"
        mkdir -p "$anat_dir"

        output_file="$anat_dir/sub-${subj_name}_T1w.nii.gz"

        echo "Converting $original_file to $output_file"
        mri_convert "$original_file" "$output_file"
        ((counter++))
    fi
done