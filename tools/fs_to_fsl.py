import argparse
import numpy as np
import nibabel as nib

def fs_to_fsl(input_fname, output_fname):
    """ 
    Reorient Tinspect MRI data in the directory to be used for VBM analysis.
    """
    if input_fname.endswith(".nii.gz"):
        ## load image
        print("##################################")
        img = nib.load(input_fname)
        data = img.get_fdata()

        ## reorient
        data = np.swapaxes(data, axis1=1, axis2=2)
        data = np.swapaxes(data, axis1=0, axis2=1)

        ## save the image
        new_img = nib.Nifti1Image(data, img.affine)
        nib.save(new_img, output_fname)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reorient MRI data for VBM analysis.")
    parser.add_argument("--input_fname", type=str, required=True, help="Path to the input NIfTI file")
    parser.add_argument("--output_fname", type=str, required=True, help="Path to save the reoriented NIfTI file")

    args = parser.parse_args()
    fs_to_fsl(args.input_fname, args.output_fname)