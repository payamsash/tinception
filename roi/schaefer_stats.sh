#!/bin/bash

# Written by Payam S. Shabestari, Zurich, 04.2025
# email: payam.sadeghishabestari@uzh.ch
# This script is written mainly for TINCEPTION project. However It could be used for other purposes.

## set Paths
export FREEiSURFER_HOME=/usr/local/freesurfer/8.0.0
export SUBJECTS_DIR=/home/ubuntu/volume/subjects_fs_dir
export LD_LBRARY_PATH=$FREESURFER_HOME/MCRv97/runtime/glnxa64:$FREESURFER_HOME/MCRv97/bin/glnxa64:$FREESURFER_HOME/MCRv97/sys/os/glnxa64:$FREESURFER_HOME/MCRv97/extern/bin/glnxa64
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=60
source $FREESURFER_HOME/SetUpFreeSurfer.sh

export PATH=/usr/lib/mrtrix3/bin:$PATH
export PATH=/home/ubuntu/fsl/bin:$PATH
export PATH=/home/ubuntu/data/src_codes/ants-2.5.4/bin:$PATH
export ANTSPATH=/home/ubuntu/data/src_codes/ants-2.5.4/bin


for subject_path in "$SUBJECTS_DIR"/*/; do
	subject=$(basename "$subject_path")
	if [ -d "$SUBJECTS_DIR/$subject/hist" ]; then
		

		mris_anatomical_stats \
            	-f "$SUBJECTS_DIR/$subject/stats/lh.Schaefer2018_1000Parcels_7Networks.stats" \
            	-b \
            	-a "$SUBJECTS_DIR/$subject/label/lh.Schaefer2018_1000Parcels_7Networks_order.annot" \
            	"$subject" lh

        	mris_anatomical_stats \
            	-f "$SUBJECTS_DIR/$subject/stats/rh.Schaefer2018_1000Parcels_7Networks.stats" \
            	-b \
            	-a "$SUBJECTS_DIR/$subject/label/rh.Schaefer2018_1000Parcels_7Networks_order.annot" \
            	"$subject" rh

		'''
		mris_anatomical_stats "$subject" lh \
  		-f "$SUBJECTS_DIR/$subject/stats/lh.Schaefer2018_1000Parcels_7Networks.stats" \
  		-b -a "$SUBJECTS_DIR/$subject/label/lh.Schaefer2018_1000Parcels_7Networks.annot"

		mris_anatomical_stats "$subject" rh \
  		-f "$SUBJECTS_DIR/$subject/stats/rh.Schaefer2018_1000Parcels_7Networks.stats" \
  		-b -a "$SUBJECTS_DIR/$subject/label/rh.Schaefer2018_1000Parcels_7Networks.annot"
		'''
	fi
done
