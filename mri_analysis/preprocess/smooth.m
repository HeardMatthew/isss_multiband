%% smooth
% Adds Gaussian smoothing to functional data for GLM processing. 
% 
% MM/DD/YY: CHANGELOG
% 02/11/20: Forked from isss_multi, made universal

function smooth(varargin)
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

if length(varargin) > 3
    flag = varargin{4};
    if strcmp(flag, 'MVPA')
        do_mvpa = 1; 
    else
        do_mvpa = 0; 
    end
    
else
    do_mvpa = 0; 
end

if ~isstruct(subj) || length(subj) ~= 1
    error('Input ("subj") where subj is a SINGLE struct')
end

if ~isstruct(study)
    error('Input ("study") which has experiment info')
end

%% Parameters
scan = study.scan(ss);  

%% Pathing
cd ..
dir_batch = pwd; 
batch = fullfile(dir_batch, 'matlabbatch', 'SPM_smooth_GLM.mat'); 
load(batch)

dir_subj = fullfile(study.path, 'data', subj.name); 

dir_func = fullfile(dir_subj, 'FUNCTIONAL'); 
dir_thissubj_batch = fullfile(dir_subj, 'batch'); %% Pathing
dir_glm  = fullfile(dir_subj, 'FUNC_GLM');
dir_mvpa = fullfile(dir_subj, 'FUNC_MVPA'); 

%% Poke values into batch
func_data = {}; 
for rr = 1:subj.runs
    [boldFiles, ~] = spm_select('List', dir_func, ['^wu' scan.runname num2str(rr) '_00*.*\.nii$']);
    boldFiles = [repmat([dir_func filesep], length(boldFiles),1), boldFiles];
    boldFiles = cellstr(boldFiles);
    func_data = vertcat(func_data, boldFiles); 
end

matlabbatch{1}.spm.spatial.smooth.data = func_data;

if do_mvpa
    matlabbatch{1}.spm.spatial.smooth.fwhm = [3 3 3]; 
else
    matlabbatch{1}.spm.spatial.smooth.fwhm = [8 8 8];
end

filename = fullfile(dir_thissubj_batch, ['smooth_GLM_' scan.runname '.mat']); 
save(filename, 'matlabbatch')

%% Prepare SPM
disp('Loading SPM...')
close all; 
spm('defaults', 'FMRI');
spm_jobman('initcfg');
spm_figure('GetWin','Graphics'); % Thanks Guillaume
disp('Done!')

%% Run the job
spm_jobman('run',matlabbatch);
disp('Done with GLM smoothing!')

%% Move files
disp('Moving files to GLM directory')
target = fullfile(dir_func, ['swu' scan.runname '*']);
if do_mvpa
    movefile(target, dir_mvpa); 
else
    movefile(target, dir_glm);
end

end