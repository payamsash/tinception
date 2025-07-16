from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import nibabel as nib
from nilearn.plotting import (plot_stat_map,
                                plot_glass_brain,
                                plot_surf_stat_map,
                                plot_roi,
                                show,
                                plot_surf_contours
                                )
import nibabel.freesurfer.io as fsio


############ VBM ############
img_dir = Path.cwd() / "material" / "with_qc" / "VBM" / "stats"
bg_image = "/Users/payamsadeghishabestari/fsl/data/standard/MNI152_T1_2mm.nii.gz"
fname_sub_mask = img_dir / "subcortical_mask_thr80.nii.gz"
fname_1 =  img_dir / "fslvbm_s3_80_tfce_corrp_tstat1.nii.gz"
fname_2 = img_dir / "fslvbm_s3_80_tfce_corrp_tstat2.nii.gz"
bg_image_1 = "/Users/payamsadeghishabestari/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz"
bg_image_2 = "/Users/payamsadeghishabestari/fsl/data/standard/MNI152_T1_0.5mm.nii.gz"
threshold = 0.95
kwargs = {
            "colorbar": True,
            "cbar_tick_format": "%.2g",
            "annotate": False,
            "draw_cross": False,
            "threshold": threshold,
            "radiological": False,
            "cmap": 'autumn',
            "symmetric_cbar": False,
            "vmin": threshold,
            "vmax": 1,
            "dim": 0.1
        }

plot_stat_map(stat_map_img=fname_1, bg_img=bg_image_1, display_mode="xz", cut_coords=None, **kwargs) # CO > TI
plot_stat_map(stat_map_img=fname_2, bg_img=bg_image_1, display_mode="x", cut_coords=(-18, -21, -24, -27, -30), **kwargs) # TI > CO
plot_stat_map(stat_map_img=fname_2, bg_img=bg_image_1, display_mode="z", cut_coords=(-8, -4, 0, 4, 8), **kwargs) # TI > CO


############ SBM ############
fsaverage_dir = Path("/Applications/freesurfer/dev/subjects/fsaverage")
root_dir = Path("./material/with_qc/SBM")

fsaverage_dir = Path("/Applications/freesurfer/dev/subjects/fsaverage")
labels, ctab, names = fsio.read_annot(fsaverage_dir / "label" / f"rh.aparc.annot")

column_names = [
    'ClusterNo', 'Max', 'VtxMax', 'Size_mm2', 
    'MNIX', 'MNIY', 'MNIZ', 
    'CWP', 'CWPLow', 'CWPHi', 
    'NVtxs', 'WghtVtx', 'Annot'
]

for hemi, he in zip(["left", "right"], ["lh", "rh"]):
    sulc = fsaverage_dir / "surf" / f"{he}.sulc"
    inflated = fsaverage_dir / "surf" / f"{he}.inflated"
    
    for mod in ["pos", "neg"]:
        stat_map = root_dir / he / "CO-TI-QC" / f"perm.th23.{mod}.sig.cluster.mgz"
        summery_fname = root_dir / he / "CO-TI-QC" / f"perm.th23.{mod}.sig.cluster.summary"

        if mod == "pos":
            input_fname = stat_map

        if mod == "neg": 
            img = nib.load(stat_map)
            data = img.get_fdata()
            data_neg = data * -1
            new_img = nib.Nifti1Image(data_neg, img.affine, img.header)
            input_fname = stat_map.parent / f"perm.th23.{mod}.sig.cluster_inverted.mgz"
            nib.save(new_img, input_fname)
            
        for view in ["medial", "lateral"]:
            figure = plot_surf_stat_map(
                                        stat_map=input_fname,
                                        surf_mesh=inflated,
                                        hemi=hemi,
                                        view=view,
                                        threshold=1.3,
                                        bg_map=sulc,
                                        darkness=0.7,
                                        bg_on_data=False,
                                        vmin=1.3,
                                        vmax=3.89,
                                        colorbar=False,
                                        cmap="autumn"
                                    )
            
            '''
            df = pd.read_csv(summery_fname, skiprows=40, delim_whitespace=True,
                            comment='#', names=column_names)
            labels = df["Annot"].values.tolist()
            levels = [names.index(f"{item}".encode()) for item in labels]

            plot_surf_contours(
                                surf_mesh=fsaverage_dir / "surf" / f"{he}.inflated",
                                roi_map=fsaverage_dir / "label" / f"{he}.aparc.annot",
                                hemi=hemi,
                                labels=labels,
                                levels=levels,
                                legend=False,
                                figure=figure,
                                colors=["cadetblue"] * len(labels),
                            )
            '''

            figure.figure.set_facecolor('black')
            figure.axes[0].set_facecolor('black')
            figure.frameon = False
            figure.axes[0].axis('off')
            figure.figure.savefig(Path.cwd() / "material" / "with_qc" / "paper_figures" / f"sbm_sig_{mod}_{hemi}_{view}_nolines.pdf",
                                facecolor=figure.figure.get_facecolor(),
                                bbox_inches='tight')
            show()

############ ROI ############
def load_freesurfer_lut_to_df(lut_path):
    columns = ['Index', 'name', 'R', 'G', 'B', 'A']
    data = []
    with open(lut_path, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.split()
            if len(parts) >= 6:
                index = int(parts[0])
                label = parts[1]
                rgb = list(map(int, parts[2:5]))
                alpha = int(parts[5])
                data.append([index, label, *rgb, alpha])
    
    return pd.DataFrame(data, columns=columns)

def rgba_to_hex(r, g, b):
    return '#{0:02x}{1:02x}{2:02x}'.format(r, g, b)

# Load sample subject image 
mri_dir = Path('/Users/payamsadeghishabestari/antinomics_clean_codes/dvob_processed/sMRI/0539/mri')
t1_path = mri_dir / "nu.mgz"
hippo_path = mri_dir / "rh.hippoAmygLabels-T1.v22.mgz"
thalamus_path = mri_dir / "ThalamicNuclei.v13.T1.mgz"
aan_path = mri_dir / "arousalNetworkLabels.v10.mgz"
lut_path = "/Applications/freesurfer/dev/FreeSurferColorLUT.txt"

## make the data label data for nilearn
seg_img = nib.load(aan_path)
seg_data = seg_img.get_fdata().astype(int)
labels = np.unique(seg_data)
label_map = {label: i for i, label in enumerate(labels)}
seg_data_mapped = np.vectorize(label_map.get)(seg_data)
seg_data_4d = seg_data_mapped[..., np.newaxis]
seg_img_4d = nib.Nifti1Image(seg_data_4d, affine=seg_img.affine, header=seg_img.header)
unique_labels = np.unique(seg_data).tolist()

## create lut file for nilearn
lut_df = load_freesurfer_lut_to_df(lut_path)
lut_df['color'] = lut_df.apply(lambda row: rgba_to_hex(row['R'], row['G'], row['B']), axis=1)
lut_df = lut_df.query('Index == @unique_labels')
lut_df.rename(columns={"Index": "index"}, inplace=True)
lut_df = lut_df[["index", "name", "color"]]
lut_df.loc[0, 'name'] = 'Background'
lut_df['index'] = range(len(lut_df))
lut_df.reset_index(drop=True, inplace=True)

figure, ax = plt.subplots(1, 1, figsize=(9, 9), layout="tight")
kwargs = {
        "bg_img": t1_path,                
        "cmap": lut_df,
        "colorbar": False,
        "annotate": False,
        "draw_cross": False,
        "black_bg": True,
        "alpha": 0.7,
        "dim" : -1,
        "figure": figure
        }
plot_roi(
        seg_img_4d,            
        display_mode="y",
        cut_coords=(7,),
        axes=ax,
        **kwargs
        )

legend_handles = [
                Patch(color=row['color'], label=row['name'], alpha=0.7)
                for _, row in lut_df.iterrows()
                ][1:]

fig, ax = plt.subplots()
fig.patch.set_facecolor('black')
ax.set_facecolor('black')

legend = ax.legend(handles=legend_handles, title="Regions", facecolor='black', edgecolor='white', frameon=False)

# Set text colors to white
plt.setp(legend.get_texts(), color='white')
plt.setp(legend.get_title(), color='white')

ax.axis('off')
plt.show()
