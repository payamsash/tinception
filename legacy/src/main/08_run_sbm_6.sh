#!/bin/bash

"""
# Surface Resampling: fsaverage → fsaverage5
#
# This script loops over all FreeSurfer subjects in SUBJECTS_DIR and resamples
# their cortical thickness maps from fsaverage to fsaverage5 for both 
# hemispheres (lh and rh). 
#
# Checks are performed to:
#   1. Skip template subjects (fsaverage, fsaverage5, fsaverage6)
#   2. Skip subjects missing the source file
#   3. Skip subjects where the target file already exists
#
# Requires:
#   - FreeSurfer environment set up
#   - Source thickness files: lh.thickness.fwhm10.fsaverage.mgh, rh.thickness.fwhm10.fsaverage.mgh
#
# Output:
#   - lh.thickness.fwhm10.fsaverage5.mgh
#   - rh.thickness.fwhm10.fsaverage5.mgh
#   in the subject's surf folder
"""

export FREESURFER_HOME=/Applications/freesurfer/dev
export SUBJECTS_DIR=/Volumes/Extreme_SSD/payam_data/Tinception/subjects_fs_dir
export LD_LIBRARY_PATH=$FREESURFER_HOME/MCRv97/runtime/glnxa64:$FREESURFER_HOME/MCRv97/bin/glnxa64:$FREESURFER_HOME/MCRv97/sys/os/glnxa64:$FREESURFER_HOME/MCRv97/extern/bin/glnxa64
source $FREESURFER_HOME/SetUpFreeSurfer.sh
export SUBJECTS_DIR

# Loop over all subjects
for subj in "$SUBJECTS_DIR"/*; do
    subj=$(basename "$subj")

    # Skip template subjects
    if [[ "$subj" == "fsaverage" || "$subj" == "fsaverage5" || "$subj" == "fsaverage6" ]]; then
        echo "Skipping template subject: $subj"
        continue
    fi

    echo "Processing subject: $subj"

    lh_src="$SUBJECTS_DIR/$subj/surf/lh.thickness.fwhm10.fsaverage.mgh"
    lh_out="$SUBJECTS_DIR/$subj/surf/lh.thickness.fwhm10.fsaverage5.mgh"

    if [[ ! -f "$lh_src" ]]; then
        echo "  LH source missing — skipping"
    elif [[ -f "$lh_out" ]]; then
        echo "  LH fsaverage5 exists — skipping"
    else
        mri_surf2surf \
            --srcsubject fsaverage \
            --trgsubject fsaverage5 \
            --hemi lh \
            --sval "$lh_src" \
            --tval "$lh_out"
    fi


    rh_src="$SUBJECTS_DIR/$subj/surf/rh.thickness.fwhm10.fsaverage.mgh"
    rh_out="$SUBJECTS_DIR/$subj/surf/rh.thickness.fwhm10.fsaverage5.mgh"

    if [[ ! -f "$rh_src" ]]; then
        echo "  RH source missing — skipping"
    elif [[ -f "$rh_out" ]]; then
        echo "  RH fsaverage5 exists — skipping"
    else
        mri_surf2surf \
            --srcsubject fsaverage \
            --trgsubject fsaverage5 \
            --hemi rh \
            --sval "$rh_src" \
            --tval "$rh_out"
    fi

done