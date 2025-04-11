#!/bin/bash

## dicom to nifti
neuropren_dir="/home/ubuntu/data/tinception/NEUROPREN"
neuropren_subjects_dir="/home/ubuntu/data/tinception/subjects_raw"

for subject_id in "$neuropren_dir"/*; do
    if [ -d "$subject_id" ]; then
        dcm2niix -o "$neuropren_subjects_dir" "$subject_id"
    fi
done

find $neuropren_subjects_dir -type f ! -name "*_t1_*.nii" -exec rm -f {} \;

for file in "$neuropren_subjects_dir"/*.nii; do
    new_name=$(echo "$file" | sed -E 's/_t1_.*//')
    mv "$file" "$new_name.nii"
done

## running recon-all in parallel
export FREESURFER_HOME=/usr/local/freesurfer/8.0.0
export SUBJECTS_DIR=/home/ubuntu/data/tinception/subjects_fs_dir
export LD_LIBRARY_PATH=$FREESURFER_HOME/MCRv97/runtime/glnxa64:$FREESURFER_HOME/MCRv97/bin/glnxa64:$FREESURFER_HOME/MCRv97/sys/os/glnxa64:$FREESURFER_HOME/MCRv97/extern/bin/glnxa64
source $FREESURFER_HOME/SetUpFreeSurfer.sh
export FS_V8_XOPTS=0

subjects_raw="/home/ubuntu/data/tinception/subjects_raw"
cd $subjects_raw
nohup bash -c 'ls *.nii | parallel --jobs 20 recon-all -s {.} -i {} -all' &

## running extra segmentations