#!/bin/bash

# -------------------------------------------------------
# FreeSurfer setup
# -------------------------------------------------------
export FREESURFER_HOME=/Applications/freesurfer/dev
source $FREESURFER_HOME/SetUpFreeSurfer.sh

# -------------------------------------------------------
# Paths
# -------------------------------------------------------
TINCEPTION_DIR=/Volumes/Extreme_SSD/payam_data/Tinception
SUBJECTS_DIR=$TINCEPTION_DIR/subjects_fs_dir
FS_SURF_DIR=$SUBJECTS_DIR/fsaverage/surf

MBM_DIR=$TINCEPTION_DIR/MBM
mkdir -p $MBM_DIR

# -------------------------------------------------------
# Convert surfaces
# -------------------------------------------------------
for hemi in lh rh
do

    IN_SURF=$FS_SURF_DIR/${hemi}.white
    OUT_VTK=$MBM_DIR/${hemi}.white.vtk

    if [ ! -f "$IN_SURF" ]; then
        echo "ERROR: missing surface $IN_SURF"
        exit 1
    fi

    echo "Converting $IN_SURF → $OUT_VTK"

    mris_convert $IN_SURF $OUT_VTK

done

echo ""
echo "Done."
echo "Created:"
echo "$MBM_DIR/lh.white.vtk"
echo "$MBM_DIR/rh.white.vtk"