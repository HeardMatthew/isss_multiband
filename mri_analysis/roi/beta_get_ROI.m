%% beta_plot_ROI_all
% Extracts LNG betas for congruent/incongruent cluster. Just runs
% hybrid/multiband. 

function beta_get_ROI(varargin)
%% Check input
subj = varargin{1}; 
study = varargin{2}; 
dd = varargin{3}; 
masks = varargin{4}; 
if length(varargin) > 4
    thismask = varargin{5}; 
else
    thismask = '';
end

if ~isstruct(study) || ~isstruct(subj)
    error('subj and study are both struct!')
end

if ~isnumeric(dd)
    error('need to specify which design!')
end

if ~ischar(masks)  
    error('masks need to be strings!')
end

if ~ischar(thismask)  
    error('thismask need to be strings!')
end

%% Pathing
dir_roi = pwd; 
cd ..

comparison = 'hybrid_multiband'; 
scans = {'hybrid', 'multiband'}; 
% dir_masks  = fullfile(dir_roi, masks, comparison); 
dir_masks  = fullfile(dir_roi, masks); 

%% Parameters
numsubjs = length(subj); 

design = study.design(dd); 
cond   = design.cond; 

dir_design_multi = fullfile(study.path, 'data', 'second_level', ['multiband_' design.name], 'LNG_NOI'); 

%% Scan-specific parameters and pathing
numbetas = 8; 
numscans = length(scans); 

%% Load ROI mask
disp('Loading hybrid mask...')

if ~strcmp(thismask, '') 
    fname = fullfile(dir_masks, thismask); 
    fname2 = fullfile(dir_masks, 'hybrid_only_nary.nii'); 
else
    fname = fullfile(dir_masks, 'hybrid_only_nary.nii'); 
    % too lazy to fix this old hard-coding...
end

Vmask = spm_vol(fname); 
map   = spm_read_vols(Vmask);

%% Which ROIs? (10/13/20)
rois = [5, 6, 9, 10];
% anterior thalamus, right stg, right thalamus, brainstem

%% Continue loading
if contains(thismask, '/users/PAS1342/osu9912/fmri/atlases/BNA_PM_3D_246/')
    % if it's a probabilistic map... 
    map(map < 25) = 0; 
    map(map >= 25) = 1; 
    map = logical(map); 
    
    Vmask2 = spm_vol(fname2); 
    map2   = spm_read_vols(Vmask2);
    map2   = logical(map2); 
    
    map = map & map2; 
end

roi_mask = logical(map); 

disp('Cross-referencing with multiband mask...') 
maskname = fullfile(dir_design_multi, 'mask.nii'); 
Vmask = spm_vol(maskname); 
multi_mask  = spm_read_vols(Vmask); multi_mask = logical(multi_mask);

combine_mask = roi_mask & multi_mask; 

map(~combine_mask) = 0; 

temp = unique(map); roi_idx = temp(temp > 0); 
numrois = length(roi_idx); 

for nscan = 1:numscans
    %% Specify scan path
    for nroi = rois %1:numrois % 10/13/20 time-saver
        % For saving data...
        if nroi < 10
            roistr = ['roi0' num2str(nroi)]; 
        else
            roistr = ['roi' num2str(nroi)]; 
        end

        betas_avg_subj_LNG = [];
        betas_std_subj_LNG = [];
        
        betas_avg_subj_NOI = [];
        betas_std_subj_NOI = [];

        %% Get voxels in this ROI
        disp(['ROI #' num2str(nroi)])
        thisroi = roi_idx(nroi); 

        mask_roi = map == thisroi;
        nvox = length(find(mask_roi)); 

        %% Preallocate
        betas_avg_LNG = zeros(numsubjs, numbetas); 
        betas_avg_NOI = zeros(numsubjs, numbetas); 

        for nsubj = 1:numsubjs
            % Subject-specific path and parameters
            thissubj = subj(nsubj); 
            disp(thissubj.name)
            numruns = thissubj.runs(1); % ehh, this works. 
            dir_design = fullfile(study.path, 'data', thissubj.name, ...
                'design', [scans{nscan} '_' design.name]); 

            % Grab some info from SPM.mat
            spmmat = fullfile(dir_design, 'SPM.mat'); 
            load(spmmat)
            beta_names = SPM.xX.name;  

            % Preallocate
            betas_roi_NOI = nan(nvox, numbetas, numruns); 
            betas_roi_OR = nan(nvox, numbetas, numruns); 
            betas_roi_SR = nan(nvox, numbetas, numruns); 
            
            for nrun = 1:numruns
                %% Identify which betas to load for this condition
                thisrun = ['Sn(' num2str(nrun) ') '];
                target_NOI = [thisrun 'NOI'];
                target_OR  = [thisrun 'OR'];
                target_SR  = [thisrun 'SR'];
                these_betas_NOI = contains(beta_names, target_NOI); 
                these_betas_OR  = contains(beta_names, target_OR); 
                these_betas_SR  = contains(beta_names, target_SR); 

                Vbetas_NOI = SPM.Vbeta(these_betas_NOI); 
                Vbetas_OR  = SPM.Vbeta(these_betas_OR); 
                Vbetas_SR  = SPM.Vbeta(these_betas_SR); 
                for nb = 1:numbetas
                    %% Load data for this condition, for each subject
                    Vbetas_NOI(nb).fname = fullfile(dir_design, Vbetas_NOI(nb).fname);
                    Vbetas_OR(nb).fname  = fullfile(dir_design, Vbetas_OR(nb).fname);
                    Vbetas_SR(nb).fname  = fullfile(dir_design, Vbetas_SR(nb).fname);
                    data_NOI = spm_read_vols(Vbetas_NOI(nb));  % loads whole brain
                    data_OR  = spm_read_vols(Vbetas_OR(nb));  % loads whole brain
                    data_SR  = spm_read_vols(Vbetas_SR(nb));  % loads whole brain
                    betas_roi_NOI(:, nb, nrun) = data_NOI(mask_roi); % grabs relevant data
                    betas_roi_OR(:, nb, nrun) = data_OR(mask_roi); % grabs relevant data
                    betas_roi_SR(:, nb, nrun) = data_SR(mask_roi); % grabs relevant data
                end

            end

            %% Get LNG results (average across OR/SR)
            betas_roi_LNG = (betas_roi_OR + betas_roi_SR) / 2; 

            %% Average across runs
            betas_avg_runs_LNG = mean(betas_roi_LNG, 3); 
            betas_avg_runs_roi_LNG = mean(betas_avg_runs_LNG, 1); 
            betas_avg_LNG(nsubj, :) = mean(betas_avg_runs_roi_LNG, 1); 
            
            betas_avg_runs_NOI = mean(betas_roi_NOI, 3); 
            betas_avg_runs_roi_NOI = mean(betas_avg_runs_NOI, 1); 
            betas_avg_NOI(nsubj, :) = mean(betas_avg_runs_roi_NOI, 1); 
        end

        %% Append results onto previous results, ready for table
        betas_avg_subj_LNG = [betas_avg_subj_LNG, mean(betas_avg_LNG, 1)];
        betas_std_subj_LNG = [betas_std_subj_LNG, std(betas_avg_LNG, 1)];
        
        betas_avg_subj_NOI = [betas_avg_subj_NOI, mean(betas_avg_NOI, 1)];
        betas_std_subj_NOI = [betas_std_subj_NOI, std(betas_avg_NOI, 1)];

        %% Save results
        COND = repelem({'LNG', 'NOI'}, numbetas)'; 
        BETA = repmat(1:numbetas, [1 2])'; 

        betas_avg_LNG = betas_avg_subj_LNG'; 
        betas_std_LNG = betas_std_subj_LNG';
        betas_ser_LNG = betas_std_LNG./sqrt(numsubjs); 
        
        betas_avg_NOI = betas_avg_subj_NOI'; 
        betas_std_NOI = betas_std_subj_NOI';
        betas_ser_NOI = betas_std_NOI./sqrt(numsubjs); 
        
        betas_avg = [betas_avg_LNG; betas_avg_NOI]; 
        betas_std = [betas_std_LNG; betas_std_NOI]; 
        betas_ser = [betas_ser_LNG; betas_ser_NOI]; 

        data = table(COND, BETA, betas_avg, betas_std, betas_ser); 

        if contains(thismask, '/users/PAS1342/osu9912/fmri/atlases/BNA_PM_3D_246/')
            % hard-coded
            fname = fullfile(dir_masks, ['betas_LNG_NOI_' scans{nscan} '_CAUDATE.mat']); 
        else
            fname = fullfile(dir_masks, ['betas_LNG_NOI_' scans{nscan} '_' roistr '.mat']); 
        end
        
        save(fname, 'data')
    end

end

end
