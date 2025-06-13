
#!/bin/bash

# Written by Payam S. Shabestari, Zurich, 04.2025
# email: payam.sadeghishabestari@uzh.ch
# This script is written mainly for TINCEPTION project. However It could be used for other purposes.

## set Paths
export FREESURFER_HOME=/usr/local/freesurfer/8.0.0
export SUBJECTS_DIR=/home/ubuntu/volume/subjects_fs_dir
export LD_LIBRARY_PATH=$FREESURFER_HOME/MCRv97/runtime/glnxa64:$FREESURFER_HOME/MCRv97/bin/glnxa64:$FREESURFER_HOME/MCRv97/sys/os/glnxa64:$FREESURFER_HOME/MCRv97/extern/bin/glnxa64
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=42
source $FREESURFER_HOME/SetUpFreeSurfer.sh

export PATH=/usr/lib/mrtrix3/bin:$PATH
export PATH=/home/ubuntu/fsl/bin:$PATH
export PATH=/home/ubuntu/data/src_codes/ants-2.5.4/bin:$PATH
export ANTSPATH=/home/ubuntu/data/src_codes/ants-2.5.4/bin

sch_gcs_dir="$FREESURFER_HOME/gcs"
n_threads=30
hemis=("lh" "rh")

for subject_path in "$SUBJECTS_DIR"/*; do
    subject_id=$(basename "$subject_path")
    if [ -d "$subject_path/hist" ]; then
        for n in 400 1000; do
            for hemi in "${hemis[@]}"; do
                mris_ca_label -l $SUBJECTS_DIR/$subject_id/label/${hemi}.cortex.label \
                                    $subject_id \
                                    ${hemi} \
                                    $SUBJECTS_DIR/$subject_id/surf/${hemi}.sphere.reg \
                                    $sch_gcs_dir/${hemi}.Schaefer2018_${n}Parcels_7Networks.gcs \
                                    $SUBJECTS_DIR/$subject_id/label/${hemi}.Schaefer2018_${n}Parcels_7Networks_order.annot
                    echo "$(date): Schaefer parcellation with $n labels in $hemi completed for subject $subject_id"
            done
        done
    fi
done