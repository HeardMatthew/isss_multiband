%% find_overlapping
% Creates list of overlapping and unique voxels in ROIs across two sets of
% ROIs. 
%
% MM/DD/YY -- CHANGELOG
% 04/21/20 -- Made file, based on Yune's code. Significant overhaul. 
% 04/22/20 -- Added n-ary output. Reworked to allow any number of ROIs

function find_overlapping(study, masks, rois)
%% check inputs
if ~isstruct(study)
    error('subj and study are both struct!')
end

if ~ischar(masks)
    error('masks need to be strings!')
end

if ~iscell(rois)
    error('rois need to be cell!')
end

if length(rois) == 1
    error('needs multiple rois!')
end

%% Pathing
dir_roi = pwd; 
dir_masks = fullfile(dir_roi, masks); 

cd ..

%% Parameters
nrois = length(rois); 

%% Load ROIs
% m{#}: binary mask
% c{#}: n-ary mask
fname = cell(nrois, 1); 
m = cell(nrois, 1);
c = cell(nrois, 1);

for nr = 1:nrois
    fname{nr} = fullfile(dir_masks, [rois{nr} '.nii']); V = spm_vol(fname{nr}); 
    c{nr} = spm_read_vols(V); 
    m{nr} = logical(c{nr});
end

%% Binary comparisons, unique and overlapping
overlapping = cell(nrois, 1); % iterates over comparisons
m_only = cell(nrois, 2);
m_only_nary = cell(nrois, 2);

num_exclusive = nan(nrois, 2); 
num_overlapping_voxels = nan(nrois, 1); 

overlapping_all = true(V.dim); 

pairs = nchoosek([1:nrois], 2); 

for nr = 1:nrois
    p1 = pairs(nr, 1); p2 = pairs(nr, 2); 
    disp(['Comparing ' rois{p1} ' and ' rois{p2}]);
    
    overlapping{nr} = m{p1} & m{p2}; 
    overlapping_all = m{nr} & overlapping_all; % for trinary comparison
    
    m_only{nr, 1} = m{p1} & ~m{p2}; 
    m_only{nr, 2} = m{p2} & ~m{p1}; 
    m_only_nary{nr, 1} = m_only{nr, 1} .* c{p1}; 
    m_only_nary{nr, 2} = m_only{nr, 2} .* c{p2}; 
    
    num_exclusive(nr, 1) = length(find(m_only{nr, 1})); 
    num_exclusive(nr, 2) = length(find(m_only{nr, 2})); 
    num_overlapping_voxels(nr) = length(find(overlapping{nr})); 
end

%% Trinary comparison
exclusive_all = cell(nrois, 1); 
num_exclusive_all = nan(nrois, 1); 
for nr = 1:nrois
    exclusive_all{nr} = m{nr} & ~overlapping_all; 
    num_exclusive_all = length(find(exclusive_all{nr})); 
end

%% Save results
comparison_name = cell(nrois, 1);

for nr = 1:nrois
    p1 = pairs(nr, 1); p2 = pairs(nr, 2); 
    
    name1 = strsplit(rois{p1}, '_'); name1 = name1{1}; 
    name2 = strsplit(rois{p2}, '_'); name2 = name2{1}; 
    dir_new = fullfile(dir_masks, [name1 '_' name2]); 
    if ~exist(dir_new, 'file'); mkdir(dir_new); end
    
    V.fname = fullfile(dir_new, [name1 '_only.nii']);
    spm_write_vol(V, m_only{nr, 1})
    V.fname = fullfile(dir_new, [name1 '_only_nary.nii']);
    spm_write_vol(V, m_only_nary{nr, 1})
    
    V.fname = fullfile(dir_new, [name2 '_only.nii']);
    spm_write_vol(V, m_only{nr, 2})
    V.fname = fullfile(dir_new, [name2 '_only_nary.nii']);
    spm_write_vol(V, m_only_nary{nr, 2})
    
    V.fname = fullfile(dir_new, [name1 '_' name2 '_overlap.nii']);  
    spm_write_vol(V, overlapping{nr});
    
    comparison_name{nr} = [name1 '_' name2];  % for writing table
end

%% Write table with number of overlapping voxels, write exclusive voxels
voxels1 = num_exclusive(:, 1); 
voxels2 = num_exclusive(:, 2); 

T = table(comparison_name, voxels1, voxels2, num_overlapping_voxels); 

fname = fullfile(dir_masks, 'data.mat'); 
save(fname, 'T', 'num_exclusive_all')

for nr = 1:nrois
    name = strsplit(rois{nr}, '_'); name = name{1}; 
    V.fname = fullfile(dir_masks, [name '_exclusive.nii']); 
    spm_write_vol(V, exclusive_all{nr}); 
end

end