#!/bin/bash

export FREESURFER_HOME=/usr/local/freesurfer/8.0.0
export SUBJECTS_DIR=/home/ubuntu/volume/subjects_fs_dir
export LD_LIBRARY_PATH=$FREESURFER_HOME/MCRv97/runtime/glnxa64:$FREESURFER_HOME/MCRv97/bin/glnxa64:$FREESURFER_HOME/MCRv97/sys/os/glnxa64:$FREESURFER_HOME/MCRv97/extern/bin/glnxa64
source $FREESURFER_HOME/SetUpFreeSurfer.sh

for subj in $(ls $SUBJECTS_DIR); do
    echo "Processing subject: $subj"

    for hemi in lh rh; do
        if [ -f "$SUBJECTS_DIR/$subj/surf/${hemi}.thickness" ]; then
            mri_surf2surf \
                --srcsubject $subj \
                --trgsubject fsaverage5 \
                --hemi $hemi \
                --sval $SUBJECTS_DIR/$subj/surf/${hemi}.thickness \
                --tval $SUBJECTS_DIR/$subj/surf/${hemi}.thickness.fsaverage5.mgh
        else
            echo "  ${hemi}.thickness not found for $subj — skipping."
        fi
    done
done