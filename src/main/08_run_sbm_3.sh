#!/bin/bash

## set Paths
export FREESURFER_HOME=/usr/local/freesurfer/8.0.0
export SUBJECTS_DIR=/home/ubuntu/volume/Tinception/subjects_fs_dir
export LD_LIBRARY_PATH=$FREESURFER_HOME/MCRv97/runtime/glnxa64:$FREESURFER_HOME/MCRv97/bin/glnxa64:$FREESURFER_HOME/MCRv97/sys/os/glnxa64:$FREESURFER_HOME/MCRv97/extern/bin/glnxa64
source $FREESURFER_HOME/SetUpFreeSurfer.sh

hemis=("lh" "rh")
smoothing_opts=("10")
measures=("area" "thickness")
sbm_dir="/home/ubuntu/volume/Tinception/SBM"

## Resampling (mris_preproc)
for hemi in "${hemis[@]}"; do
    for smoothing in "${smoothing_opts[@]}"; do
        for measure in "${measures[@]}"; do
            mris_preproc --fsgd $sbm_dir/FSGD/tinception.fsgd \
                        --cache-in ${measure}.fwhm${smoothing}.fsaverage \
                        --target fsaverage \
                        --hemi $hemi \
                        --out $sbm_dir/results/${hemi}.${measure}.${smoothing}.mgh
        done
    done
done

## GLM analysis
for hemi in "${hemis[@]}"; do
    for smoothing in "${smoothing_opts[@]}"; do
        for measure in "${measures[@]}"; do
            mri_glmfit --y $sbm_dir/results/$hemi.$measure.$smoothing.mgh \
                        --fsgd $sbm_dir/FSGD/tinception.fsgd \
                        --C $sbm_dir/FSGD/co_gt_ti.mtx \
                        --C $sbm_dir/FSGD/ti_gt_co.mtx \
                        --surf fsaverage $hemi \
                        --cortex \
                        --glmdir $sbm_dir/results/$hemi.$measure.$smoothing.glmdir
        done
    done
done

## cluster correction
for hemi in "${hemis[@]}"; do
    for smoothing in "${smoothing_opts[@]}"; do
        for measure in "${measures[@]}"; do
            current_glmdir="$sbm_dir/results/$hemi.$measure.$smoothing.glmdir"
            mri_glmfit-sim --glmdir $current_glmdir \
                            --cache 2.0 pos \
                            --cwp 0.05 \
                            --2spaces
        done
    done
done

## RH
mri_segstats --i ../../rh.thickness.10.mgh --seg cache.th20.pos.sig.ocn.mgh --id 1 --avgwf rh.cluster1.txt
mri_segstats --i ../../rh.thickness.10.mgh --seg cache.th20.pos.sig.ocn.mgh --id 2 --avgwf rh.cluster2.txt
mri_segstats --i ../../rh.thickness.10.mgh --seg cache.th20.pos.sig.ocn.mgh --id 3 --avgwf rh.cluster3.txt

## LH
mri_segstats --i ../../lh.thickness.10.mgh --seg cache.th20.pos.sig.ocn.mgh --id 1 --avgwf rh.cluster1.txt
mri_segstats --i ../../lh.thickness.10.mgh --seg cache.th20.pos.sig.ocn.mgh --id 2 --avgwf rh.cluster2.txt


'''
## visualize
surf_dir="/Applications/freesurfer/dev/subjects/fsaverage/surf"
smoothing=10
measure="thickness"
mode="ti_gt_co"
freeview -f $surf_dir/lh.inflated:overlay=lh.$measure.$smoothing.glmdir/$mode/cache.th20.pos.sig.cluster.mgh  \
            -f $surf_dir/rh.inflated:overlay=rh.$measure.$smoothing.glmdir/$mode/cache.th20.pos.sig.cluster.mgh
'''
