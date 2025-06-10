#!/bin/bash

VBM_dir="/home/ubuntu/volume/VBM"
DATA_DIR=$VBM_dir/"struc"

OUTPUT_FILE="tiv_results.txt"
echo -e "SubjectID\tTIV_mm3" > "$OUTPUT_FILE"


for csf_file in "$DATA_DIR"/*_struc_brain_pve_0.nii.gz; do
    fname=$(basename "$csf_file")
    subject_id=${fname%%_struc_brain_pve_0.nii.gz}

    wm_file="${DATA_DIR}/${subject_id}_struc_brain_pve_2.nii.gz"
    gm_file="${DATA_DIR}/${subject_id}_struc_GM.nii.gz"

    ## checks
    if [[ -f "$gm_file" && -f "$wm_file" ]]; then
        csf_vol=$(fslstats "$csf_file" -V | awk '{print $2}')
        gm_vol=$(fslstats "$gm_file" -V | awk '{print $2}')
        wm_vol=$(fslstats "$wm_file" -V | awk '{print $2}')
        tiv=$(echo "$csf_vol + $gm_vol + $wm_vol" | bc)

        # Save to file
        echo -e "${subject_id}\t${tiv}" >> "$OUTPUT_FILE"
        echo "${subject_id} done."
    else
        echo "Missing files for subject ${subject_id}, skipping." >&2
    fi
done

echo "Done. Results saved to $OUTPUT_FILE"