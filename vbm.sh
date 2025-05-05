#!/bin/bash

# Written by Payam S. Shabestari, Zurich, 04.2025
# email: payam.sadeghishabestari@uzh.ch
# This script is written mainly for TINCEPTION project. However It could be used for other purposes.


## create VBM folder
export FREESURFER_HOME=/usr/local/freesurfer/8.0.0
export SUBJECTS_DIR=/home/ubuntu/volume/subjects_fs_dir
source $FREESURFER_HOME/SetUpFreeSurfer.sh

VBM_DIR="/home/ubuntu/volume/VBM"
ID_FILE="/home/ubuntu/data/src_codes/matched_subjects.txt"  # list of matched subject IDs
mkdir -p "$VBM_DIR"

while IFS= read -r subj; do
    src_file="$SUBJECTS_DIR/$subj/mri/orig/001.mgz"
    dest_file="$VBM_DIR/${subj}.nii.gz"
    temp_nii="/home/ubuntu/volume/temp/${subj}.nii"
    
    if [ -f "$src_file" ]; then
        mri_convert "$src_file" "$temp_nii"
        fslreorient2std "$temp_nii" "$dest_file"
        rm "$temp_nii"
        echo "Saved $dest_file"

    else
        echo "Warning: $src_file not found"
    fi
done < "$ID_FILE"


cd $VBM_DIR
fslvbm_1_bet -b

'''
## checks here slicesdir
fslvbm_2_template -n

## create design.mat with text2vest (demean age, 3 dummy var for site and 2 cols for group)
'''
1 0 -6 1 0 0 0
1 0 -3 0 1 0 0
1 0  4 1 0 1 0
1 0  0 0 0 0 1
0 1 -5 1 0 0 0
0 1  3 1 1 0 0
0 1  6 0 0 1 0
0 1 -2 0 0 0 1
...
'''

## next create design.con with: Text2Vest contrasts.txt design.con

-1  1  0  0  0  0  0  # Tinnitus > Control
 1 -1  0  0  0  0  0  # Control > Tinnitus


-1  1  0  0
1 -1  0  0


### for paper plot
freeview -f \
  $SUBJECTS_DIR/subject/surf/lh.pial:curvature=curv \
  -viewport 3d

'''


