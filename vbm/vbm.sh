#!/bin/bash

# Written by Payam S. Shabestari, Zurich, 04.2025
# email: payam.sadeghishabestari@uzh.ch
# This script is written mainly for TINCEPTION project. However It could be used for other purposes.


## create VBM folder
export FREESURFER_HOME=/usr/local/freesurfer/8.0.0
export SUBJECTS_DIR=/home/ubuntu/volume/subjects_fs_dir
source $FREESURFER_HOME/SetUpFreeSurfer.sh

VBM_b_DIR="/home/ubuntu/volume/VBM_b"
VBM_n_DIR="/home/ubuntu/volume/VBM_n"
temp="/home/ubuntu/volume/temp"
ID_FILE="/home/ubuntu/data/src_codes/matched_subjects.txt"  # list of matched subject IDs
mkdir -p "$VBM_b_DIR"
mkdir -p "$VBM_n_DIR"
mkdir -p "$temp"

while IFS= read -r subj; do
    src_file="$SUBJECTS_DIR/$subj/mri/orig/001.mgz"
    temp_nii="/home/ubuntu/volume/temp/${subj}.nii" 

    if [ -f "$src_file" ]; then 

        if [[ "$subj" =~ ^[0-9]{2} ]]; then # tinmeg
            dest_file="$VBM_n_DIR/${subj}.nii.gz"
            if [ ! -f "$dest_file" ]; then
                mri_convert "$src_file" "$temp_nii"
                fslreorient2std "$temp_nii" "$dest_file"
                rm "$temp_nii"
                echo "Saved $dest_file"
            fi

        elif [[ "$subj" == ICC* ]]; then # talaska
            dest_file="$VBM_b_DIR/${subj}.nii.gz"
            if [ ! -f "$dest_file" ]; then
                mri_convert "$src_file" "$temp_nii"
                fslreorient2std "$temp_nii" "$dest_file"
                rm "$temp_nii"
                echo "Saved $dest_file"
            fi

        elif [[ "$subj" == *MRT || "$subj" == *MRT_2 ]]; then # neuropren
            dest_file="$VBM_n_DIR/${subj}.nii.gz"
            if [ ! -f "$dest_file" ]; then
                mri_convert "$src_file" "$temp_nii"
                fslreorient2std "$temp_nii" "$dest_file"
                rm "$temp_nii"
                echo "Saved $dest_file"
            fi

        elif [[ "$subj" == *NT-HL || "$subj" == *NT-HA || "$subj" == *TI-HL || "$subj" == *TI-HA ]]; then
            dest_file="$VBM_n_DIR/${subj}.nii.gz"
            if [ ! -f "$dest_file" ]; then
                mri_convert "$src_file" "$temp_nii"
                fslreorient2std "$temp_nii" "$dest_file"
                rm "$temp_nii"
                echo "Saved $dest_file"
            fi

        elif [[ "$subj" == MR* ]]; then # tinspect
            dest_file="$VBM_b_DIR/${subj}.nii.gz"
            if [ ! -f "$dest_file" ]; then
                mri_convert "$src_file" "$temp_nii"
                fslreorient2std "$temp_nii" "$dest_file"
                python fs_to_fsl.py --input_fname "$dest_file" --output_fname "$dest_file"
                rm "$temp_nii"
                echo "Saved $dest_file"
            fi
        
        fi
    else
        echo "Source file not found for subject $subj"
    
    fi

done < "$ID_FILE"


cd $VBM_DIR
fslvbm_1_bet -b
fslvbm_2_template -n # for some data -R -f 0.6
fslvbm_3_proc

nohup bash -c 'randomise -i GM_mod_merg_s3 -m GM_mask -o fslvbm_s3 -d design.mat -t design.con -T -n 5000; \
                randomise -i GM_mod_merg_s2 -m GM_mask -o fslvbm_s2 -d design.mat -t design.con -T -n 5000; \
                randomise -i GM_mod_merg_s4 -m GM_mask -o fslvbm_s4 -d design.mat -t design.con -T -n 5000' \
                > out.log 2>&1 &

## just using subcortical VBM
nohup bash -c 'randomise -i GM_mod_merg_s2 -m GM_mask -o fslvbm_s2_80 -m subcortical_mask_thr80 -d design.mat -t design.con -T -n 5000; \
                randomise -i GM_mod_merg_s2 -m GM_mask -o fslvbm_s3_80 -m subcortical_mask_thr80 -d design.mat -t design.con -T -n 5000' \
                > out2.log 2>&1 &

fsleyes $FSLDIR/data/standard/MNI152_T1_2mm fslvbm_s3_80_tfce_corrp_tstat1 -cm red-yellow -dr 0.949 1


'''
## creating Hippo mask
fslmaths $FSLDIR/data/standard/MNI152_T1_2mm -mul 0 -add 1 -roi 59 1 54 1 27 1 0 1 hippo_L_point -odt float
fslmaths $FSLDIR/data/standard/MNI152_T1_2mm -mul 0 -add 1 -roi 31 1 54 1 27 1 0 1 hippo_R_point -odt float
fslmaths hippo_L_point -kernel sphere 20 -fmean hippo_L_roi -odt float
fslmaths hippo_R_point -kernel sphere 20 -fmean hippo_R_roi -odt float
fslmaths hippo_L_roi -add hippo_R_roi hippo_bilat_sphere
fslmaths hippo_bilat_sphere -mul 2 -thr `fslstats hippo_bilat_sphere -p 100` -bin sphere_mask
fslmaths hippo_bilat_sphere -mas GM_mask hippo_mask
'''