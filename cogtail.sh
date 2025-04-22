src_base="/home/ubuntu/data/tinception/cogtail/CogTAiL/ICCAC"
dest="/home/ubuntu/data/tinception/subjects_raw"

mkdir -p "$dest"

for subj_dir in "$src_base"/*/; do
    subj_id=$(basename "$subj_dir")
    nii_file=$(find "$subj_dir" -maxdepth 1 -name "*.nii" | head -n 1)
    if [ -n "$nii_file" ]; then
        cp "$nii_file" "$dest/${subj_id}.nii"
        echo "Moved $nii_file to $dest/${subj_id}.nii"
    else
        echo "No .nii file found in $subj_dir"
    fi
done
# NT-HA  NT-HL  TI-HA  TI-HL ICCAC
