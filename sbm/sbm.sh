#!/bin/bash

# Written by Payam S. Shabestari, Zurich, 04.2025
# email: payam.sadeghishabestari@uzh.ch
# This script is written mainly for TINCEPTION project. However It could be used for other purposes.

## set Paths
export FREESURFER_HOME=/usr/local/freesurfer/8.0.0
export SUBJECTS_DIR=/home/ubuntu/volume/subjects_fs_dir
export LD_LIBRARY_PATH=$FREESURFER_HOME/MCRv97/runtime/glnxa64:$FREESURFER_HOME/MCRv97/bin/glnxa64:$FREESURFER_HOME/MCRv97/sys/os/glnxa64:$FREESURFER_HOME/MCRv97/extern/bin/glnxa64
source $FREESURFER_HOME/SetUpFreeSurfer.sh


## running qcache in parallel ...
for subj in "$SUBJECTS_DIR"/*; do
    if [ -d "$subj" ]; then
        name=$(basename "$subj")
        if [[ "$name" != "fsaverage" && "$name" != "MNI152" ]]; then
            echo "Running qcache on subject: $name"
            recon-all -s "$name" -qcache
        fi
    fi
done


hemis=("lh" "rh")
smoothing_opts=("5" "10" "15")
measures=("volume" "area" "thickness")
sbm_dir="/home/ubuntu/volume/SBM"

## resampling surface or volume source data to a common subject
for hemi in "${hemis[@]}"; do
    for smoothing in "${smoothing_opts[@]}"; do
        for measure in "${measures[@]}"; do
            mris_preproc --fsgd $sbm_dir/FSGD/tinception_fsgd.fsgd \
                        --cache-in {$measure}.fwhm{$smoothing}.fsaverage \
                        --target fsaverage \
                        --hemi $hemi \
                        --out $sbm_dir/results/{$hemi}.{$measure}.{$smoothing}.mgh
        done
    done
done

## GLM analysis in the volume or the surface
for hemi in "${hemis[@]}"; do
    for smoothing in "${smoothing_opts[@]}"; do
        for measure in "${measures[@]}"; do
            mri_glmfit --y $sbm_dir/results/$hemi.$measure.$smoothing.mgh \
                        --fsgd $sbm_dir/FSGD/tinception_fsgd.fsgd \
                        --C $sbm_dir/Contrasts/CO-TI.mtx \
                        --C $sbm_dir/Contrasts/TI-CO.mtx \
                        --surf fsaverage $hemi \
                        --cortex \
                        --glmdir $sbm_dir/results/$hemi.$measure.$smoothing.glmdir
        done
    done
done

## cluster correction (1.3, 2, 2.3, 3)
for hemi in "${hemis[@]}"; do
    for smoothing in "${smoothing_opts[@]}"; do
        for measure in "${measures[@]}"; do
            for dir in "$sbm_dir/results/$hemi.$measure.$smoothing.glmdir"; do
                mri_glmfit-sim --glmdir $dir \
                                --cache 1.3 pos \ 
                                --cwp 0.05 \
                                --2spaces
            done
        done
    done
done


## visualize
surf_dir="/Applications/freesurfer/dev/subjects/fsaverage/surf"
smoothing=10
measure="thickness"
mode="neg"
freeview -f $surf_dir/lh.inflated:overlay=lh.$measure.$smoothing.glmdir/TI-CO-corr/sig.mgh:overlay_threshold=1.3,10  \
            -f $surf_dir/rh.inflated:overlay=rh.$measure.$smoothing.glmdir/TI-CO-corr/sig.mgh:overlay_threshold=1.3,10
