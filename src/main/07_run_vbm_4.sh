#!/bin/bash
fslvbm_3_proc
nohup bash -c 'randomise -i GM_mod_merg_s3 -m GM_mask -o fslvbm_s3 -d design.mat -t design.con -T -n 5000; \
                randomise -i GM_mod_merg_s3 -m GM_mask -o fslvbm_s3_80 -m subcortical_mask_thr80 -d design.mat -t design.con -T -n 5000; \
                > out.log 2>&1 &