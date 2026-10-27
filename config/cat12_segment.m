% CAT12.9 (r2665) volume-only segmentation for Tinception; run by slurm/05_cat12.sbatch through
% cat_standalone.sh, which fills in the input file. The placeholder string must appear ONLY on the
% data line below: cat_standalone.sh takes the first line containing it as the data field.
% Fields not set here keep the cat_defaults.m values: default TPM and the MNI152NLin2009cAsym
% Geodesic Shooting template shipped with CAT12 (standard space, so UKB can be added later).
% Parameters mirror config/analysis.yaml -> cat12.

matlabbatch{1}.spm.tools.cat.estwrite.data = '<UNDEFINED>';

matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.vox = 1.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.lazy = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.ignoreErrors = 1;  % failures leave no mwp1 -> caught by the sbatch
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.verb = 2;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.print = 2;         % volume-only report pdf

% no surfaces (FreeSurfer covers the cortex); outputs in mri/ report/ label/ next to the input
matlabbatch{1}.spm.tools.cat.estwrite.output.surface = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.BIDS.BIDSno = 1;

% ROI volumes (label/catROI_*.xml): whole brain, coarse, and SUIT cerebellum
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.neuromorphometrics = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.cobra = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.suit = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.lpba40 = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.hammers = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.ibsr = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.aal3 = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.mori = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.thalamus = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.anatomy3 = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.julichbrain = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.ownatlas = {''};

% modulated (affine + nonlinear) GM and WM in template space: mwp1*, mwp2*
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.mod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.dartel = 0;

% QC images: native label map (p0*) and the warped bias-corrected T1 (wm*)
matlabbatch{1}.spm.tools.cat.estwrite.output.label.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.label.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.label.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.warped = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.dartel = 0;

matlabbatch{1}.spm.tools.cat.estwrite.output.jacobianwarped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.warps = [0 0];
