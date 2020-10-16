%% roi_snr
% runs ROI analysis on SNR data, prepares for paired t-test across scans

function roi_snr(subj, study, masks, scans)
%% check inputs
if ~isstruct(study) || ~isstruct(subj)
    error('subj and study are both struct!')
end

if ~ischar(masks)  
    error('masks need to be strings!')
end

if ~iscell(scans)
    error('scans need to be cell!')
end

%% Pathing
dir_roi = pwd; 
dir_masks = fullfile(dir_roi, masks); 

cd ..
dir_batch  = pwd; 
comparison = strjoin(scans, '_'); 
dir_masks  = fullfile(dir_roi, masks, comparison); 

%% Parameters
nsubjs = length(subj); 
nscans = length(scans); 

%% Load masks
fnames = cell(nscans, 1); 
Vmask = cell(nscans, 1); 
mask = cell(nscans, 1); 

for ns = 1:nscans
    fnames{ns} = fullfile(dir_masks, [scans{ns} '_only.nii']);
    Vmask{ns} = spm_vol(fnames{ns}); 
    mask{ns} = spm_read_vols(Vmask{ns}); mask{ns} = logical(mask{ns}); 
end

%% Get voxels from each subject
mean_SNR = nan(nscans, nsubjs); 

for ii = 1:nsubjs
    for jj = 1:nscans
        thissubj = subj(ii); 
        dir_snr  = fullfile(study.path, 'data', thissubj.name, 'SNR'); 
        thisfile = fullfile(dir_snr, [study.prefix scans{jj} '_averaged_snr.nii']); 
        
        V = spm_vol(thisfile); data = spm_read_vols(V); 
        data = data(mask{ns}); mean_SNR(jj, ii) = mean(data); 
    end
    
end

%% Save the averages
fname = fullfile(dir_masks, 'mean_SNR_unique.mat'); 
save(fname, 'mean_SNR')

end