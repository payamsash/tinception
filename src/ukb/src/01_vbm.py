from pathlib import Path
import numpy as np
from nilearn.plotting import plot_stat_map

"""
This script visualizes and summarizes VBM (voxel-based morphometry) results.

- Generates statistical maps (TFCE-corrected) for significant group differences.
- Reconstructs deviation maps from normative models for subcortical regions.
- Identifies significant clusters, extracts GM volumes per subject, and computes t-values and Cohen's d.
- Creates publication-ready boxplots of GM volumes across groups and brain regions.
- Saves all figures to a dedicated plots directory.
"""

############ VBM stat ############
tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
fsl_dir = Path("/Users/payamsadeghishabestari/fsl")
vbm_dir = tinception_dir / "vbm_results" / "UKB"
vbm_norm_dir = tinception_dir / "vbm_norm"
plots_dir = tinception_dir / "plots"
plots_dir.mkdir(exist_ok=True)

img_bg = fsl_dir / "data" / "standard" / "MNI152_T1_0.5mm.nii.gz"
img_fname =  vbm_dir / "TI_gt_CO_logp_max_tfce.nii" # TI > CO

cuts = [(-16.07, 8.95, -8.32), (-15.32, -67.59, -44.34)]
pvals = [0.05, 0.005]
for p_thr, cut in zip(pvals, cuts):

    thr = -np.log10(p_thr)

    kwargs = {
                "colorbar": True,
                "cbar_tick_format": "%.2g",
                "annotate": False,
                "draw_cross": False,
                "threshold": thr,
                "radiological": False,
                "cmap": 'autumn',
                "symmetric_cbar": False,
                "vmin": 1.3,
                "vmax": 2.5,
                "dim": -0.3,
                "black_bg": True,
                "cut_coords": cut
            }
    for disp_mode in ["x", "y", "z"]:

        fig = plot_stat_map(
                        stat_map_img=img_fname,
                        bg_img=img_bg,
                        display_mode=disp_mode,
                        **kwargs
                        )
        fig.savefig(plots_dir / "vbm_group" / "ukb" / f"ti_gt_co_{disp_mode}_thr_{p_thr}.pdf", dpi=600, bbox_inches='tight')

