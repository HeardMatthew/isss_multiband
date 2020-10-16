%% coregister_roi
% Coregister ROI masks for all subjects. Doing so ahead of time should make
% the main analysis run faster. 

close all; clear all; clc;
tic 
%% Paths and parameters
dir_batch_roi = pwd;
cd ..
isss_multi_params
cd ..
addpath spm12

spm('defaults', 'FMRI');
spm_jobman('initcfg');

scan(1).name = 'hybrid'; 
scan(1).num_TC_per_run = 180;
scan(1).num_TR = 10;
scan(1).num_runs = 2;

scan(2).name = 'isss'; 
scan(2).num_TC_per_run = 90;
scan(2).num_TR = 5;
scan(2).num_runs = 2;

cd(dir_batch_roi)
rois_nii = dir('*.nii');
if length(rois_nii) ~= 3
    error('Check nii files')
end

numROI = length(rois_nii);

for ROI=1:numROI % For each ROI...
    disp(rois_nii(ROI).name)
    
    for sc = 1:length(scan) % For each scan type... 
        disp(scan(sc).name)
        thisscan = scan(sc); 
        
        for ss = 1:length(subjects) % For each subject...
            thissubj = subjects{ss};        
            disp(thissubj);
            dir_design = fullfile(dir_data, thissubj, 'design', thisscan.name);
            cd(dir_design);

            this_roi = fullfile(rois_nii(ROI).folder, rois_nii(ROI).name);

            cd(dir_batch_roi)
            load SPM_coregister_nn.mat
            matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {fullfile(dir_design, 'mask.nii')}; 
            matlabbatch{1}.spm.spatial.coreg.estwrite.source = {this_roi};
            batchname = fullfile(dir_data, thissubj, 'batch', ['coregister_roi_' rois_nii(ROI).name '.mat']);
            save(batchname, 'matlabbatch')

            spm_jobman('run',matlabbatch);

            mask_old = fullfile(dir_batch_roi, ['r' rois_nii(ROI).name]);
            mask_v = spm_vol(mask_old);
            this_mask = fullfile(dir_batch_roi, 'mask_coreg', [thissubj '_' thisscan.name '_r' rois_nii(ROI).name]);

            movefile(mask_old, this_mask)
            
        end
        
    end
    
end

toc

disp('Done!')

        
        
        
        