%% beta_plot_ROI_all
% Extracts LNG betas for congruent/incongruent cluster. Just runs
% hybrid/multiband. 

function SPA_get_ROI(varargin)
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
dir_masks  = fullfile(dir_roi, masks, comparison); 

%% Parameters
numsubjs = length(subj); 

design = study.design(dd); 
cond   = design.cond; 

dir_design_multi = fullfile(study.path, 'data', 'second_level', ['multiband_' design.name], 'LNG_NOI'); 

%% Scan-specific parameters and pathing
numscans = length(scans); 

%% Load ROI mask
disp('Loading hybrid mask...')

if ~strcmp(thismask, '') 
    fname = thismask; 
    fname2 = fullfile(dir_masks, 'hybrid_only_nary.nii'); 
else
    fname = fullfile(dir_masks, 'hybrid_only_nary.nii'); 
    % too lazy to fix this old hard-coding...
end

Vmask = spm_vol(fname); 
map   = spm_read_vols(Vmask);

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

for nroi = 1:numrois
    % For saving data...
    disp(['ROI #' num2str(nroi)])
    
    if nroi < 10
        roistr = ['roi0' num2str(nroi)]; 
    else
        roistr = ['roi' num2str(nroi)]; 
    end
    
    %% Get voxels in this ROI
    thisroi = roi_idx(nroi); 
    mask_roi = map == thisroi;
    
    % Preallocate
    SPA_avg = []; 
    SPA_std = []; 
        
    for nscan = 1:numscans
        %% Preallocate
        SPA_avg_LNG = zeros(numsubjs, 1); 
        SPA_avg_NOI = zeros(numsubjs, 1); 

        for nsubj = 1:numsubjs
            % Subject-specific path and parameters
            thissubj = subj(nsubj); 
            disp(thissubj.name)
            numruns = thissubj.runs(1); % ehh, this works. 
            dir_design = fullfile(study.path, 'data', thissubj.name, ...
                'design', [scans{nscan} '_' design.name]); 

            % Load SPA for subjects 
            fname = fullfile(dir_design, 'AUE_LNG.nii'); 
            V = spm_vol(fname); 
            data = spm_read_vols(V); 
            SPA_avg_LNG_subj(nsubj) = mean(data(mask_roi)); 
            
            fname = fullfile(dir_design, 'AUE_NOI.nii'); 
            V = spm_vol(fname); 
            data = spm_read_vols(V); 
            SPA_avg_NOI_subj(nsubj) = mean(data(mask_roi)); 
        end

        %% Average across subjects
        SPA_avg_LNG = mean(SPA_avg_LNG_subj);
        SPA_std_LNG = std(SPA_avg_LNG_subj);
        
        SPA_avg_NOI = mean(SPA_avg_NOI_subj); 
        SPA_std_NOI = std(SPA_avg_NOI_subj); 
        
        %% Append results onto previous results, ready for table
        SPA_avg = [SPA_avg; SPA_avg_LNG; SPA_avg_NOI]; 
        SPA_std = [SPA_std; SPA_std_LNG; SPA_std_NOI]; 
    end
    
    %% Save results
    SCAN = repelem({'Hybrid'; 'Multi'}, 2);  
    COND = repmat({'LNG'; 'NOI'}, [2 1]); 

    SPA_ser = SPA_std./sqrt(numsubjs); 

    data = table(SCAN, COND, SPA_avg, SPA_std, SPA_ser); 
    
    if contains(thismask, '/users/PAS1342/osu9912/fmri/atlases/BNA_PM_3D_246/')
        % hard-coded
        fname = fullfile(dir_masks, 'SPA_LNG_NOI_CAUDATE.mat'); 
    else
        fname = fullfile(dir_masks, ['SPA_LNG_NOI_' roistr '.mat']); 
    end
    
    save(fname, 'data')
end

end
