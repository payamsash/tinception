#!/bin/bash

export FREESURFER_HOME=/Applications/freesurfer/dev
export SUBJECTS_DIR=/Users/payamsadeghishabestari/antinomics_clean_codes/dvob_processed/sMRI
source $FREESURFER_HOME/SetUpFreeSurfer.sh

subject_id="0539"
lut_dir=$FREESURFER_HOME/average/AAN/atlas/freeview.lut.txt
structure="AAN"
cd $SUBJECTS_DIR/$subject_id/mri


if [[ "$structure" == "AAN" ]]; then
    freeview -v T1.mgz -v arousalNetworkLabels.v10.mgz:colormap=lut:lut=$lut_dir
fi

if [[ "$structure" == "hippo" ]]; then
    freeview -v nu.mgz -v T1.mgz:sample=cubic \
                -v lh.hippoAmygLabels-T1.v22.mgz:colormap=lut \
                -v rh.hippoAmygLabels-T1.v22.mgz:colormap=lut
fi

if [[ "$structure" == "thalamus" ]]; then
    freeview -v nu.mgz -v ThalamicNuclei.v13.T1.mgz:colormap=lut
fi