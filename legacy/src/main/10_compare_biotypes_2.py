from pathlib import Path
import pandas as pd
import numpy as np
from nilearn import plotting, image, masking
from nilearn.glm.second_level import SecondLevelModel
from nilearn.glm import threshold_stats_img
from nilearn.image import resample_to_img

## compare biotype 1 from Tinnitus with all controls

## lets plot Voxel-wise Contrast
tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
vbm_dir = tinception_dir / "VBM" / "struc"
df_biotype = pd.read_csv(tinception_dir / "biotypes" / "main_TI.csv")
df = pd.read_csv(tinception_dir / "VBM_design" / "covars.csv")
df = df.sort_values("tinception_id").reset_index(drop=True)

## creating covars
df['age'] = df['age'] - df['age'].mean()
site_dummies = pd.get_dummies(df['site'], prefix='site', drop_first=True).astype(int)
df['PTA'] = df['PTA'] - df['PTA'].mean()
df['TIV'] = df['TIV'] - df['TIV'].mean()

df = pd.concat([df, site_dummies], axis=1)
df.drop(columns=["subject_ID", "site", "THI", "distance", "weights", "subclass"], inplace=True)
df.rename(columns={"tinception_id": "subjects"}, inplace=True)

## create design matrix
df = df.merge(df_biotype, on="subjects", how="left")
conditions = [
    (df['group'] == 'CO'),
    (df['group'] == 'TI') & (df['Biotype'] == 0),
    (df['group'] == 'TI') & (df['Biotype'] == 1)
]
choices = ['CO', 'TI_0', 'TI_1']
df['final_label'] = np.select(conditions, choices, default='Unknown')
df.drop(columns=["group", "Biotype"], inplace=True)
df.rename(columns={"final_label": "group"}, inplace=True)
df.sort_values(by="subjects", inplace=True)
df_0 = df.query('group == "CO"')

subjects_to_drop = ["sub-07006", "sub-08016", "sub-09021", "sub-10015", "sub-10017", "sub-10025"] # outliers

for biotype_id in ["TI_0", "TI_1"][:1]:

    df_1 = df.query(f'group == "{biotype_id}"') # its here
    df_design = pd.concat([df_0, df_1], axis=0)
    df_design = df_design.query('subjects != @subjects_to_drop')
    df_design.drop(columns=["subjects", "group"], inplace=True)

    ## loop over files
    group_0 = list(df_0["subjects"].values)
    group_1 = list(df_1["subjects"].values)

    ##############
    drop_set = set(subjects_to_drop)
    group_0 = [s for s in group_0 if s not in drop_set]
    group_1 = [s for s in group_1 if s not in drop_set]
    ##############

    print(f"len of group_0 is {len(group_0)}")
    print(f"len of group_1 is {len(group_1)}")

    files_b0 = []
    files_b1 = []
    for file in sorted(vbm_dir.iterdir()):
        name = file.name
        if name.startswith("."):
            continue
        else:
            if name.endswith("_struc_GM_to_template_GM_mod.nii.gz"):
                subject = name[:9]
                if subject in group_0:
                    files_b0.append(file)
                if subject in group_1:
                    files_b1.append(file)

    ## create a big image
    img_b0 = image.concat_imgs(files_b0)
    img_b1 = image.concat_imgs(files_b1)

    print(f"len of files_b0 is {len(files_b0)}")
    print(f"len of files_b1 is {len(files_b1)}")

    mean_b0 = image.mean_img(img_b0)
    mean_b1 = image.mean_img(img_b1)

    var_b0 = image.new_img_like(mean_b0, np.var(img_b0.get_fdata(), axis=-1))
    var_b1 = image.new_img_like(mean_b1, np.var(img_b1.get_fdata(), axis=-1))
    pooled_std = image.new_img_like(mean_b0, np.sqrt((var_b0.get_fdata() + var_b1.get_fdata()) / 2))

    cohens_d_data = (mean_b0.get_fdata() - mean_b1.get_fdata()) / (pooled_std.get_fdata() + 1e-8)
    cohens_d_img = image.new_img_like(mean_b0, cohens_d_data)

    n_b0 = img_b0.shape[3]
    n_b1 = img_b1.shape[3]

    ## final desgin matrix
    design_matrix = pd.DataFrame(
        [1] * n_b0 + [0] * n_b1,
        columns=['Group_0'] )
    design_matrix['Group_1'] = [0] * n_b0 + [1] * n_b1
    design_matrix = pd.concat([design_matrix, df_design.reset_index(drop=True)], axis=1)

    all_imgs = image.concat_imgs([img_b0, img_b1])
    model = SecondLevelModel(smoothing_fwhm=3.0)
    model.fit(all_imgs, design_matrix=design_matrix)
    z_map = model.compute_contrast('Group_0 - Group_1', output_type='z_score')

    z_map_masked_img = z_map
    d_map_masked_img = cohens_d_img
    d_map_smoothed = image.smooth_img(d_map_masked_img, fwhm=6.0)

    ## FDR
    thresholded_map_1, threshold_1 = threshold_stats_img(
                                                    z_map_masked_img, 
                                                    alpha=0.005, 
                                                    height_control='fdr',
                                                    cluster_threshold=10,
                                                    two_sided=True   
                                                    )
    thresholded_map_2, threshold_2 = threshold_stats_img(
                                                        z_map_masked_img, 
                                                        alpha=0.05, 
                                                        height_control='fdr',
                                                        cluster_threshold=300, 
                                                        two_sided=True,
                                                    )

    data_pos = thresholded_map_2.get_fdata().copy() 
    data_pos[data_pos < 0] = 0


    ## plotting
    fsl_dir = Path("/Users/payamsadeghishabestari/fsl")
    img_bg = fsl_dir / "data" / "standard" / "MNI152_T1_0.5mm.nii.gz"
    brain_mask = masking.compute_brain_mask(img_bg)
    brain_mask_resampled = resample_to_img(brain_mask, thresholded_map_2, interpolation='nearest')
    mask_data = brain_mask_resampled.get_fdata()
    data_pos = data_pos * (mask_data > 0)
    thr_pos = image.new_img_like(thresholded_map_2, data_pos)

    kwargs = {
                "colorbar": False,
                "cbar_tick_format": "%.2g",
                "annotate": False,
                "draw_cross": True,
                "radiological": False,
                "cmap": 'magma',
                "symmetric_cbar": False,
                "vmin": -7,
                "dim": -0.3,
                "black_bg": True,
                "cut_coords": 3
            }


    for disp_mode in ["x", "y", "z"]:
        fig = plotting.plot_stat_map(
                        stat_map_img=thresholded_map_1,
                        threshold=threshold_1,
                        bg_img=img_bg,
                        vmax=-5,
                        display_mode=disp_mode,
                        **kwargs
                        )
        fig.savefig(tinception_dir / "plots" / f"VBM_biotype_{biotype_id}_vs_co_neg_{disp_mode}.pdf", dpi=600, bbox_inches='tight')
        
        fig = plotting.plot_stat_map(
                        stat_map_img=thr_pos,
                        threshold=threshold_2,
                        bg_img=img_bg,
                        vmax=3,
                        display_mode=disp_mode,
                        **kwargs
                        )
        fig.savefig(tinception_dir / "plots" / f"VBM_biotype_{biotype_id}_vs_co_pos_{disp_mode}.pdf", dpi=600, bbox_inches='tight')