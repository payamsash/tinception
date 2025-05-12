#!/bin/bash

# Written by Payam S. Shabestari, Zurich, 04.2025
# email: payam.sadeghishabestari@uzh.ch
# This script is written mainly for TINCEPTION project. However It could be used for other purposes.

# setup FS
export FREESURFER_HOME=/usr/local/freesurfer/8.0.0
export SUBJECTS_DIR=/home/ubuntu/volume/subjects_fs_dir
export LD_LIBRARY_PATH=$FREESURFER_HOME/MCRv97/runtime/glnxa64:$FREESURFER_HOME/MCRv97/bin/glnxa64:$FREESURFER_HOME/MCRv97/sys/os/glnxa64:$FREESURFER_HOME/MCRv97/extern/bin/glnxa64
source $FREESURFER_HOME/SetUpFreeSurfer.sh
export FS_V8_XOPTS=0

###### NEUROPREN
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

subjects_raw="/home/ubuntu/data/tinception/subjects_raw"
cd $subjects_raw
nohup bash -c 'ls *.mgz | parallel --jobs 20 recon-all -s {.} -i {} -all' &

###### TINSPECT & TINMEG
## .hdr to nifti
tinspect_dir="/home/ubuntu/data/tinception/TINSPECT"
tinspect_subjects_dir="/home/ubuntu/data/tinception/subjects_raw"

## find T1 folder and convert it
for subject_path in "$tinspect_dir"/*; do
    if [ -d "$subject_path" ]; then
        subject_id=$(basename "$subject_path")
        echo -e "\e[33mconverting subject $subject_id ..."
        t1w_dir=$(find "$subject_path" -maxdepth 1 -type d -name '*_T1W_*' | head -n 1)
        img_file=$(find "$t1w_dir" -maxdepth 1 -type f -name '*.img' | head -n 1)
        mri_convert $img_file "$tinspect_subjects_dir/$subject_id.mgz"
    fi
done

subjects_raw="/home/ubuntu/data/tinception/subjects_raw"
cd $subjects_raw
nohup bash -c 'ls *.nii | parallel --jobs 21 "subject={.}; [[ \$subject == MR* ]] && flag=\"-notal-check\" || flag=\"\"; recon-all -s \$subject -i {} -all \$flag"' &

###### COGTAIL
cogtail_dir="/home/ubuntu/data/tinception/cogtail/CogTAiL/ICCAC" # NT-HA  NT-HL  TI-HA  TI-HL ICCAC
cogtail_subjects_dir="/home/ubuntu/data/cogtail/subjects_raw"

mkdir -p "$dest"
for subj_dir in "$cogtail_dir"/*/; do
    subj_id=$(basename "$subj_dir")
    nii_file=$(find "$subj_dir" -maxdepth 1 -name "*.nii" | head -n 1)
    if [ -n "$nii_file" ]; then
        cp "$nii_file" "$cogtail_subjects_dir/${subj_id}.nii"
        echo "Moved $nii_file to $cogtail_subjects_dir/${subj_id}.nii"
    else
        echo "No .nii file found in $subj_dir"
    fi
done

subjects_raw="/home/ubuntu/data/tinception/subjects_raw"
cd $subjects_raw
nohup bash -c 'ls *.nii | parallel --jobs 25 recon-all -s {.} -i {} -all' &


# nohup bash -c "ls *.nii.gz | parallel --jobs 42 'recon-all -s {= s/\.nii\.gz$// =} -i {} -all'" > recon_all.log 2>&1 &
