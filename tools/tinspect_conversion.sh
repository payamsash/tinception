#!/bin/bash

tinspect_dir="/Volumes/G_USZ_ORL$/Research/ANT/tinception/TINSPECT"
tinspect_subjects_dir="/Users/payamsadeghishabestari/antinomics_clean_codes/subjects_raw"

# for subject_path in "$tinspect_dir"/*; do
#     if [ -d "$subject_path" ] && [[ $(basename "$subject_path") == MR* ]]; then
#         subject_id=$(basename "$subject_path")
#         echo -e "\e[33mConverting subject $subject_id ...\e[0m"
#         t1w_dir=$(find "$subject_path" -maxdepth 1 -type d -name '*_T1W_*' | head -n 1)
#         img_file=$(find "$t1w_dir" -maxdepth 1 -type f -name '*.hdr' | head -n 1)
#         python fs_to_fsl.py --input_fname "$img_file" --output_fname "$tinspect_subjects_dir/$subject_id.nii.gz"
#     fi
# done



####
# for subject_path in "$tinspect_dir"/*; do
#     if [ -d "$subject_path" ] && [[ $(basename "$subject_path") == MR* ]]; then
#         subject_id=$(basename "$subject_path")
#         echo -e "\e[33mConverting subject $subject_id ...\e[0m"
#         t1w_dir=$(find "$subject_path" -maxdepth 1 -type d -name '*_T1W_*' | head -n 1)
#         img_file=$(find "$t1w_dir" -maxdepth 1 -type f -name '*.img' | head -n 1)
#         mri_convert -iid 0 1 0 -ijd 0 0 1 -ikd 1 0 0 $img_file "$tinspect_subjects_dir/$subject_id.nii.gz"
#     fi
# done


for subject_path in "$tinspect_dir"/*; do
    if [ -d "$subject_path" ] && [[ $(basename "$subject_path") == MR* ]]; then
        subject_id=$(basename "$subject_path")
        echo -e "\e[33mConverting subject $subject_id ...\e[0m"
        t1w_dir=$(find "$subject_path" -maxdepth 1 -type d -name '*_T1W_*' | head -n 1)
        img_file=$(find "$t1w_dir" -maxdepth 1 -type f -name '*.hdr' | head -n 1)
        python fs_to_fsl_2.py --input_fname "$img_file" --output_fname "$tinspect_subjects_dir/$subject_id.nii.gz"
    fi
done