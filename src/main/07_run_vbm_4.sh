#!/bin/bash
fslvbm_3_proc
nohup bash -c 'randomise -i GM_mod_merg_s3 -m GM_mask -o fslvbm_s3 -d design.mat -t design.con -T -n 5000; \
                randomise -i GM_mod_merg_s3 -m GM_mask -o fslvbm_s3_80 -m subcortical_mask_thr80 -d design.mat -t design.con -T -n 5000; \
                > out.log 2>&1 &





nohup bash -c '
    # Run 1: Smoothing 3, Standard Mask, Cluster-based thresholding
    randomise -i ../GM_mod_merg_s3 -m ../GM_mask -o fslvbm_s3 -d ../design.mat -t ../design.con -c 2.3 -n 5000; \
    
    # Run 2: Smoothing 3, Subcortical Mask, Cluster-based thresholding
    randomise -i ../GM_mod_merg_s3 -m ../subcortical_mask_thr80 -o fslvbm_s3_80 -d ../design.mat -t ../design.con -c 2.3 -n 5000; \
    
    # Run 3: Smoothing 4, Standard Mask, TFCE (Threshold-Free Cluster Enhancement)
    randomise -i ../GM_mod_merg_s4 -m ../GM_mask -o fslvbm_s4 -d ../design.mat -t ../design.con -T -n 5000; \
    
    # Run 4: Smoothing 4, Subcortical Mask, TFCE
    randomise -i ../GM_mod_merg_s4 -m ../subcortical_mask_thr80 -o fslvbm_s4_80 -d ../design.mat -t ../design.con -T -n 5000' \
    > out.log 2>&1 &