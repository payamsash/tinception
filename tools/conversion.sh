#!/bin/bash

# Written by Payam S. Shabestari, Zurich, 04.2025
# email: payam.sadeghishabestari@uzh.ch
# This script is written mainly for TINCEPTION project. However It could be used for other purposes.

# setup FS
export FREESURFER_HOME=/usr/local/freesurfer/8.0.0
export SUBJECTS_DIR=/home/ubuntu/volume/Tinception/subjects_fs_dir
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

###### INDITMS
subjects_raw="/Users/payamsadeghishabestari/MRI_indiTMS"
dest_folder="/Users/payamsadeghishabestari/MRI_indiTMS_nifti"
cd "$subjects_raw" || exit
for folder in */; do
    dcm2niix "$folder"
    folder_name="${folder%/}"
    mv $folder/"${folder_name}"_t1_MPRAGE*.nii $dest_folder
done

##### TRIPLE
BASE_DIR="/Users/payamsadeghishabestari/Downloads"
OUT_DIR="$BASE_DIR/TRIPLE"
mkdir -p "$OUT_DIR"

for group in patients controls; do
    GROUP_DIR="$BASE_DIR/$group"

    for subj_dir in "$GROUP_DIR"/*; do
        subj_id=$(basename "$subj_dir")

        # --- controls structure ---
        if [ "$group" == "controls" ]; then
            nifti_dir="$subj_dir/nifti"
            hdr_file=$(find "$nifti_dir" -maxdepth 1 -name "s${subj_id}_vbm.hdr" 2>/dev/null)
        
        # --- patient structure ---
        elif [ "$group" == "patients" ]; then
            nifti_dir="$subj_dir/pre/nifti"
            hdr_file=$(find "$nifti_dir" -maxdepth 1 -name "s${subj_id}_pre-vbm.hdr" 2>/dev/null)
        fi

        # Process if .hdr file found
        if [ -f "$hdr_file" ]; then
            img_file="${hdr_file%.hdr}.img"
            nii_file="$OUT_DIR/${subj_id}_${group}.nii.gz"

            # Skip if already done
            if [ -f "$nii_file" ]; then
                echo "Skipping $nii_file (already exists)"
                continue
            fi

            echo "Processing $subj_id ($group)..."

            # Convert .hdr/.img → .nii (temporary)
            fslchfiletype NIFTI "$hdr_file" /tmp/tmp_${subj_id}.nii

            # Reorient to standard and save to TRIPLE folder
            fslreorient2std /tmp/tmp_${subj_id}.nii "$nii_file"

            # Remove temporary file
            rm /tmp/tmp_${subj_id}.nii
        else
            echo "No .hdr found for $subj_id in expected location"
        fi
    done
done


# nohup bash -c "ls *.nii.gz | parallel --jobs 42 'recon-all -s {= s/\.nii\.gz$// =} -i {} -all'" > recon_all.log 2>&1 &
