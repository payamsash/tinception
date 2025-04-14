#!/bin/bash

# Written by Payam S. Shabestari, Zurich, 04.2025
# email: payam.sadeghishabestari@uzh.ch
# This script is written mainly for TINCEPTION project. However It could be used for other purposes.

## set Paths
export FREESURFER_HOME=/usr/local/freesurfer/8.0.0
export SUBJECTS_DIR=/home/ubuntu/data/subjects_fs_dir
export LD_LIBRARY_PATH=$FREESURFER_HOME/MCRv97/runtime/glnxa64:$FREESURFER_HOME/MCRv97/bin/glnxa64:$FREESURFER_HOME/MCRv97/sys/os/glnxa64:$FREESURFER_HOME/MCRv97/extern/bin/glnxa64
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=25
source $FREESURFER_HOME/SetUpFreeSurfer.sh

export PATH=/usr/lib/mrtrix3/bin:$PATH
export PATH=/home/ubuntu/fsl/bin:$PATH
export PATH=/home/ubuntu/data/src_codes/ants-2.5.4/bin:$PATH
export ANTSPATH=/home/ubuntu/data/src_codes/ants-2.5.4/bin
sch_gcs_dir="$FREESURFER_HOME/gcs"

## the real part
for subject_path in "$SUBJECTS_DIR"/*; do
    subject_id=$(basename "$subject_path")
    echo -e "\e[32msMRI processing of subject $subject_id started at $(date '+%Y-%m-%d %H:%M:%S')"
    start_time=$(date +%s)

    # hippocampal subfields and nuclei of the amygdala
    echo -e "\e[32mSegmentation of hippocampal subfields and nuclei of the amygdala!"
    segmentHA_T1.sh $subject_id 

    # brainstem
    echo -e "\e[32mSegmentation of Brainstem Substructures!"
    segmentBS.sh $subject_id

    # thalamic nuclei
    echo -e "\e[32mSegmentation of thalamic nuclei using only T1 image!"
    segmentThalamicNuclei.sh $subject_id

    # ascending arousal network
    echo -e "\e[32mSegmentations of brainstem nuclei that are part of the Ascending Arousal Network!"
    sudo chmod +x $FREESURFER_HOME/bin/segmentNuclei
    segmentAAN.sh $subject_id

    # hypothalamus
    echo -e "\e[32mSegmentation of the hypothalamus and its associated subunits"
    mri_segment_hypothalamic_subunits --s $subject_id

    # striatum
    echo -e "\e[32mStriatal parcellation!"
    mask_options=("Tight" "Loose")
    for mask_option in "${mask_options[@]}"; do
        mri_vol2vol --mov $SUBJECTS_DIR/$subject_id/mri/norm.mgz \
                    --s $subject_id \
                    --targ $SUBJECTS_DIR/MNI152/choi_atlas/17_network_${mask_option}_mask.nii.gz \
                    --m3z $SUBJECTS_DIR/MNI152/mri/transforms/talairach.m3z \
                    --noDefM3zPath \
                    --o $SUBJECTS_DIR/$subject_id/mri/striatum_17_network_${mask_option}_mask.nii.gz \
                    --inv-morph \
                    --interp nearest
    done
    
    # cerebellum
    echo -e "\e[32mCerebellum parcellation!"
    mask_options=("Tight" "Loose")
    for mask_option in "${mask_options[@]}"; do
        mri_vol2vol --mov $SUBJECTS_DIR/$subject_id/mri/norm.mgz \
                    --s $subject_id \
                    --targ $SUBJECTS_DIR/MNI152/buckner_atlas/17_network_${mask_option}_mask.nii.gz \
                    --m3z $SUBJECTS_DIR/MNI152/mri/transforms/talairach.m3z \
                    --noDefM3zPath \
                    --o $SUBJECTS_DIR/$subject_id/mri/cerebellum_17_network_${mask_option}_mask.nii.gz \
                    --inv-morph \
                    --interp nearest
    done

    # schaefer atlas
    echo -e "\e[32mSchaefer2018 parcellation in individual surface space!"
    mkdir $SUBJECTS_DIR/$subject_id/schaefer
    hemis=("lh" "rh")
    for n in 200 400 600 800 1000; do
        for net_option in 7 17; do
            for hemi in "${hemis[@]}"; do
                mris_ca_label -l $SUBJECTS_DIR/$subject_id/label/${hemi}.cortex.label \
                                $subject_id \
                                ${hemi} \
                                $SUBJECTS_DIR/$subject_id/surf/${hemi}.sphere.reg \
                                $sch_gcs_dir/${hemi}.Schaefer2018_${n}Parcels_${net_option}Networks.gcs \
                                $SUBJECTS_DIR/$subject_id/label/${hemi}.Schaefer2018_${n}Parcels_${net_option}Networks_order.annot
            done
        done
    done

    # histo atlas
    echo -e "\e[32mBayesian Segmentation with Histological Atlas "NextBrain""
    mkdir $SUBJECTS_DIR/$subject_id/hist
    mri_histo_atlas_segment_fireants $SUBJECTS_DIR/$subject_id/mri/T1.mgz \
                                        $SUBJECTS_DIR/$subject_id/hist \
                                        0 25

done