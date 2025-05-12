import argparse
import numpy as np
import nibabel as nib

def fs_to_fsl(input_fname, output_fname):
    """ 
    Reorient Tinspect MRI data in the directory to be used for VBM analysis.
    """
    if input_fname.endswith(".hdr"):
        ## load image
        print("##################################")
        img = nib.load(input_fname)
        data = img.get_fdata()
        affine = img.affine.copy()

        ## reorient
        data = np.swapaxes(data, axis1=1, axis2=2)
        data = np.swapaxes(data, axis1=0, axis2=1)

        affine[:, [1, 2]] = affine[:, [2, 1]]
        affine[:, [0, 1]] = affine[:, [1, 0]]

        ## update header
        new_header = img.header.copy()
        new_header.set_data_shape(data.shape)

        ## save the image
        new_img = nib.Nifti1Image(data, affine, header=new_header)
        nib.save(new_img, output_fname)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reorient MRI data for VBM analysis.")
    parser.add_argument("--input_fname", type=str, required=True, help="Path to the input NIfTI file")
    parser.add_argument("--output_fname", type=str, required=True, help="Path to save the reoriented NIfTI file")

    args = parser.parse_args()
    fs_to_fsl(args.input_fname, args.output_fname)