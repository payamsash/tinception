import argparse
import numpy as np
import nibabel as nib
from nibabel.orientations import (io_orientation,
                                    axcodes2ornt,
                                    ornt_transform,
                                    apply_orientation,
                                    aff2axcodes)


def fs_to_fsl_2(input_fname, output_fname):
    """ 
    Reorient Tinspect MRI data in the directory to be used for VBM analysis.
    """
    if input_fname.endswith(".hdr"):
        ## load image
        print("##################################")
        img = nib.load(input_fname)
        data = img.get_fdata()
        affine = img.affine
        header = img.header

        ## reorient
        data = np.swapaxes(data, axis1=1, axis2=2)
        data = np.swapaxes(data, axis1=0, axis2=1)
        
        orientation = aff2axcodes(affine)
        print("Orientation (world axes for voxel axes):", orientation)

        current_ornt = io_orientation(affine)
        target_ornt = axcodes2ornt(('R', 'A', 'S'))
        transform = ornt_transform(current_ornt, target_ornt)
        reoriented_data = apply_orientation(data, transform)
        new_affine = img.as_reoriented(transform).affine
        reoriented_img = nib.Nifti1Image(reoriented_data, new_affine, header)

        orientation = aff2axcodes(new_affine)
        print("Orientation (world axes for voxel axes):", orientation)
        nib.save(reoriented_img, output_fname)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reorient MRI data for FS analysis.")
    parser.add_argument("--input_fname", type=str, required=True, help="Path to the input NIfTI file")
    parser.add_argument("--output_fname", type=str, required=True, help="Path to save the reoriented NIfTI file")

    args = parser.parse_args()
    fs_to_fsl_2(args.input_fname, args.output_fname)