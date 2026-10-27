## matlab script to create masks

subjects_dir = '/Volumes/Extreme_SSD/payam_data/Tinception/subjects_fs_dir';
mbm_dir = '/Volumes/Extreme_SSD/payam_data/Tinception/MBM';

hemis = {'lh','rh'};

for i = 1:2

    hemi = hemis{i};

    label_file = fullfile(subjects_dir,'fsaverage','label',[hemi '.cortex.label']);

    label = read_label([], label_file);

    mask = zeros(163842,1); % fsaverage vertex count
    mask(label(:,1)+1) = 1;

    mask_file = fullfile(mbm_dir,[hemi '_mask.txt']);

    dlmwrite(mask_file, mask, 'delimiter',' ');

end