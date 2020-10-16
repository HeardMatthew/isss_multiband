%% coregister
% Coregister structural and functional images

% MM/DD/YY: CHANGELOG
% 02/11/20: Forked from isss_multi, turned into universal function. 
% 02/21/20: New dataset with fewer realigned images. Updated!
% 03/30/20: Cloned for isss_multi. New input strategy to identify which run
%   is being analyzed. 

function coregister(varargin)
%% Check input
if length(varargin) < 2
    error('Requires (subj, study) input, (ss) is optional!')
end

subj = varargin{1}; 
study = varargin{2}; 
if length(varargin) > 2
    ss = varargin{3};
else
    ss = 1; 
end

if ~isstruct(subj) || length(subj) ~= 1
    error('Input ("subj") where subj is a SINGLE struct')
end

if ~isstruct(study)
    error('Input ("study") which has experiment info')
end

%% Parameters
scan = study.scan(ss);  
numruns = subj.runs(ss); 

%% Pathing
cd ..
batch = fullfile(pwd, 'matlabbatch', 'SPM_coregister.mat'); 
load(batch)

dir_subj = fullfile(study.path, 'data', subj.name); 
dir_anat = fullfile(dir_subj, 'ANATOMICAL'); 
dir_func = fullfile(dir_subj, 'FUNCTIONAL'); 
dir_thissubj_batch = fullfile(dir_subj, 'batch'); 

%% Poke into matlabbatch
anat = fullfile(dir_anat, study.anat); 
matlabbatch{1}.spm.spatial.coreg.estimate.source = cellstr(anat);

funcfiles = dir(fullfile(dir_func, '*.nii')); 
funcfiles = {funcfiles(:).name}'; 
meanimg = ['^meanu' scan.runname]; 
meanimg = funcfiles(~cellfun(@isempty, regexp(funcfiles, meanimg))); 
if length(meanimg) > 1
    if scan.first == 1
        meanimg = meanimg{1}; 
    elseif length(meanimg) == 2
        meanimg = meanimg{2}; 
    else
        error('Too many meanimg!')
    end
    
end

reference = fullfile(dir_func, meanimg); 
matlabbatch{1}.spm.spatial.coreg.estimate.ref = cellstr(reference);

batchname = fullfile(dir_thissubj_batch, ['coregister_meanu_' scan.runname '.mat']);
save(batchname, 'matlabbatch')

%% Prepare SPM
disp('Loading SPM...')
close all; 
spm('defaults', 'FMRI');
spm_jobman('initcfg');
spm_figure('GetWin','Graphics'); % Thanks Guillaume
disp('Done!')

%% Run job
spm_jobman('run',matlabbatch); 
    
end