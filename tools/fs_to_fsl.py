import os
from pathlib import Path
import numpy as np
import nibabel as nib

def fs_to_fsl(subject):
    """ 
    Reorient Tinspect MRI data in the directory to be used for VBM analysis.
    """
    for fname in os.listdir(directory):
        if fname.endswith(".nii") and fname.startswith("MR"):
            ## load image
            fname = Path(fname)
            img = nib.load(fname)
            data = img.get_fdata()

            ## reorient
            data = np.swapaxes(data, axis1=1, axis2=2)
            data = np.swapaxes(data, axis1=0, axis2=1)

            ## save the image
            new_img = nib.Nifti1Image(data, img.affine)
            nib.save(new_img, fname.parent / "VBM" / fname.name)

