%% get_wb_accuracy
% Calculate mean, SD, and SER across whole brain MVPA results
% 
% MM/DD/YY -- CHANGELOG
% 05/08/20 -- Initialized. 
% 09/14/20 -- Forcing this to call manuscript dataset because I'm at my
%   wit's end. Added "scans" to input

function get_sig_accuracy(subj, study, dd, classifier, scans)
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

if ~ischar(classifier)
    error('please specify which classifier!')
end

%% Pathing
dir_mvpa = pwd; 
cd ..

% dir_masks  = fullfile(dir_roi, masks, comparison); 

%% Parameters
numsubjs = length(subj); 
design   = study.design(dd); 
numscans = length(scans); 
% cond   = design.cond; 

%% Load second-level whole brain masks
disp('Loading second level masks...')

% Preallocate
dir_designs = cell(2, 1); 

for nscan = 1:numscans
    dir_designs{nscan} = fullfile(study.path, 'data', '050520_manuscript', ['MVPA_' scans{nscan} '_' classifier]); 
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
ACC_avg_all = []; 
ACC_std_all = []; 

ACC_table = nan(numsubjs, numscans); 

for nscan = 1:numscans
    %% And then just look at significant voxels (this is whole brain...)
    fname = fullfile(dir_designs{nscan}, 'results.nii');
    V = spm_vol(fname); 
    y = spm_read_vols(V); y = logical(y); 
    sig_mask = brain_mask & y; 
    
    %% Preallocate
    ACC = zeros(numsubjs, 1); 

    for nsubj = 1:numsubjs
        % Subject-specific path and parameters
        thissubj = subj(nsubj); 
        disp(thissubj.name)
        numruns = thissubj.runs(1); % ehh, this works. 
        dir_design = fullfile(study.path, 'data', thissubj.name, 'MVPA'); 
        
        % Load ACC for subjects 
        fname = [thissubj.name '_' scans{nscan} '_' design.name '_beta_' classifier '_rad' num2str(study.mvpa.radius) '.nii']; 
        fname = fullfile(dir_design, fname); 
        V = spm_vol(fname); 
        data = spm_read_vols(V); 
        ACC(nsubj) = mean(data(sig_mask)); 
    end
    
    ACC_table(:, nscan) = ACC; 

    %% Average across subjects
    ACC_avg = mean(ACC);
    ACC_std = std(ACC);

    %% Append results onto previous results, ready for table
    ACC_avg_all = [ACC_avg_all; ACC_avg]; 
    ACC_std_all = [ACC_std_all; ACC_std]; 
end

%% Save results
ACC_avg = ACC_avg_all; 
ACC_std = ACC_std_all; 
ACC_ser = ACC_std_all./sqrt(numsubjs); 

data = table(scans', ACC_avg, ACC_std, ACC_ser); 
scannames = strjoin(scans, '_'); 
fname = fullfile(dir_mvpa, ['ACC_' classifier '_' scannames '_sigclusters_091420.mat']); 

save(fname, 'data', 'ACC_table')

end