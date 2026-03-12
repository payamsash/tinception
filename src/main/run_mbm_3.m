repoDir = '/Users/payamsadeghishabestari/MBM';
dataDir = '/Volumes/Extreme_SSD/payam_data/Tinception/MBM';

addpath(genpath(repoDir));

designFile = fullfile(dataDir,'thickness_design_ancova.txt');

hemis = {'lh','rh'};

for i = 1:length(hemis)

    hemi = hemis{i};

    MBM = struct();

    % maps
    MBM.maps.anatListFile = fullfile(dataDir,[hemi '_thickness_map_list.txt']);
    MBM.maps.maskFile     = fullfile(dataDir,[hemi '_mask.txt']);

    % statistics
    MBM.stat.test = 'ANCOVA_F';
    MBM.stat.designFile = designFile;
    MBM.stat.nPer = 5000;
    MBM.stat.pThr = 0.1;
    MBM.stat.thres = 0.05;
    MBM.stat.fdr = 'true';

    % eigenmodes
    MBM.eig.nEigenmode = 200;
    MBM.eig.saveResult = 1;
    MBM.eig.resultFile = fullfile(dataDir,[hemi '_mbm_results.mat']);

    % plotting
    MBM.plot.vtkFile = fullfile(dataDir,[hemi '.white.vtk']);
    MBM.plot.visualize = 1;
    MBM.plot.saveFig = 1;
    MBM.plot.figFile = fullfile(dataDir,[hemi '_mbm_result.png']);
    MBM.plot.nInfluentialMode = 10;

    if strcmp(hemi,'lh')
        MBM.plot.hemis = 'left';
    else
        MBM.plot.hemis = 'right';
    end

    fprintf('Running MBM for %s...\n',hemi);

    MBM = mbm_main(MBM);

end