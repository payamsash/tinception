from pathlib import Path
import numpy as np
import nibabel as nib
from nilearn import datasets

## fetch atlas
atlas_data = datasets.fetch_atlas_harvard_oxford("sub-maxprob-thr0-2mm", symmetric_split=False)
atlas_img = nib.load(atlas_data["filename"])
atlas = atlas_img.get_fdata()
affine = atlas_img.affine

## load atlas
atlas_prob = datasets.fetch_atlas_harvard_oxford('sub-prob-2mm')
mask_data = np.zeros(nib.load(atlas_prob["filename"]).shape[:3])
atlas_img_prob = nib.load(atlas_prob["filename"])
atlas_data_prob = atlas_img_prob.get_fdata()

## thr at 80%
for i in range(atlas_data_prob.shape[3]):
    roi = atlas_data_prob[..., i]
    roi_mask = (roi >= 80).astype(int)
    mask_data += roi_mask
mask_data = (mask_data > 0).astype(np.uint8)

## save it
vbm_dir = "."
mask_img = nib.Nifti1Image(mask_data, affine)
nib.save(mask_img, Path(vbm_dir) / "stats" / "subcortical_mask_thr80.nii.gz")
