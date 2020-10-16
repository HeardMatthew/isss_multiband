%% get_wb_accuracy
% Calculate mean, SD, and SER across whole brain MVPA results
% 
% MM/DD/YY -- CHANGELOG
% 05/08/20 -- Initialized. 
% 09/14/20 -- Forcing this to call manuscript dataset because I'm at my
%   wit's end. Added "scans" to input
% 09/24/20 -- Forked to check whole brain SPA values across scans. Manually
%   coding the contrast of interest!
% 2. Load each subject's SPA (contrast)
% 3. Compute average SPA 
% 4. Write to file

function compare_wb_spa(subj, study, dd)
%% check inputs
if length(subj) < 2
    error('Submit ALL subjects!')
end

if ~isstruct(subj) || ~isstruct(study)
    error('subj and study are both struct!')
end

if ~isnumeric(dd)
    error('dd specify number of design!')
end

%% Pathing
dir_process = pwd; 
cd ..
dir_batch = pwd; 

%% Parameters
scans = {'hybrid', 'isss', 'multiband'}; 
numsubjs = length(subj); 
design   = study.design(dd); 
numscans = 3; 

%% Load second-level whole brain masks
disp('Loading second level masks...')

% Preallocate
dir_designs = cell(numscans, 1); 

for nscan = 1:numscans
    dir_designs{nscan} = fullfile(study.path, 'data', '050520_manuscript', [scans{nscan} '_' design.name], 'LNG_NOI'); 
    fname = fullfile(dir_designs{nscan}, 'mask.nii'); % WHOLE-BRAIN MASK
    Vmask = spm_vol(fname); [mask, xyz] = spm_read_vols(Vmask); 
    mask = logical(mask);
    if nscan == 1
        brain_mask = mask; 
    else
        brain_mask = brain_mask & mask; 
    end
    
end

voxels_in_mask = xyz(:, brain_mask(:)); % useful for debugging

%% Get accuracy per subject
% Preallocate
SPA_avg_all = []; 
SPA_std_all = []; 
SPA_table = nan(numsubjs, numscans); 

for nscan = 1:numscans
    %% And then just look at significant voxels (this is whole brain...)
%     fname = fullfile(dir_designs{nscan}, 'results.nii');
%     V = spm_vol(fname); 
%     y = spm_read_vols(V); y = logical(y); 
%     sig_mask = brain_mask & y; 
    
    %% Preallocate
    SPA = zeros(numsubjs, 1); 

    for nsubj = 1:numsubjs
        % Subject-specific path and parameters
        thissubj = subj(nsubj); 
        disp(thissubj.name)
        dir_design = fullfile(study.path, 'data', thissubj.name, 'design', [scans{nscan} '_' design.name]); 
        
        % Load SPA for subjects 
        fname = fullfile(dir_design, 'AUE_LNG_NOI.nii'); 
        V = spm_vol(fname); 
        data = spm_read_vols(V); 
        SPA(nsubj) = mean(data(brain_mask)); 
    end
    
    SPA_table(:, nscan) = SPA; 

    %% Average across subjects
    SPA_avg = mean(SPA);
    SPA_std = std(SPA);

    %% Append results onto previous results, ready for table
    SPA_avg_all = [SPA_avg_all; SPA_avg]; 
    SPA_std_all = [SPA_std_all; SPA_std]; 
end

%% Save results
SPA_avg = SPA_avg_all; 
SPA_std = SPA_std_all; 
SPA_ser = SPA_std_all./sqrt(numsubjs); 

data = table(scans', SPA_avg, SPA_std, SPA_ser); 
% scannames = strjoin(scans, '_'); 
fname = fullfile(dir_batch, 'SPA_whole_brain_092420.mat'); 
save(fname, 'data', 'SPA_table')

end