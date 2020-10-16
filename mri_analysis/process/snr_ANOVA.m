%% second_level_AUE
% Runs 1 sample t-test and/or displays results 

function snr_ANOVA(subj, study)
%% check inputs
if length(subj) < 2
    error('Submit ALL subjects!')
end

if ~isstruct(subj) || ~isstruct(study)
    error('subj and study are both struct!')
end

%% Parameters
scans = study.scan;
numscans = length(scans); 
numsubjs = length(subj); 

%% Prepare SPM
disp('Loading SPM...')
close all; 
spm('defaults', 'FMRI');
spm_jobman('initcfg');
spm_figure('GetWin','Graphics'); % Thanks Guillaume
disp('Done!')

%% Pathing
cd ..
dir_batch  = pwd; 
dir_mlb    = fullfile(dir_batch, 'matlabbatch'); 
dir_second = fullfile(study.path, 'data', 'second_level'); 

batch = fullfile(dir_mlb, 'SPM_anova_snr.mat');
batch2 = fullfile(dir_mlb, 'SPM_estimate_GLM.mat'); 

%% Name the directory 
dir_design = fullfile(dir_second, [study.prefix '_SNR_ANOVA']); 
if ~exist(dir_design, 'file'); mkdir(dir_design); end

spmmat = fullfile(dir_design, 'SPM.mat'); 
if exist(spmmat, 'file'); delete(spmmat); end

load(batch)
matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(dir_design);

for nsubj = 1:numsubjs
    thissubj = subj(nsubj);
    disp(['Loading subj ' thissubj.name ' into batch.'])
    dir_snr = fullfile(study.path, 'data', thissubj.name, 'SNR');
    
    SNRfiles = cell(numscans, 1);
    for nscan = 1:numscans
        %% Load SNR files
        scan = scans(nscan); 
        scanname = scan.runname(1:end-4); 
        
        target = [study.prefix scanname '_averaged_snr.nii']; 
        SNRfiles{nscan} = fullfile(dir_snr, target);
    end
    
    matlabbatch{1}.spm.stats.factorial_design.des.anovaw.fsubject(nsubj).scans = SNRfiles; 
    matlabbatch{1}.spm.stats.factorial_design.des.anovaw.fsubject(nsubj).conds = [1 2 3]; 
end

%% Save the batch
filename = fullfile(dir_design, 'SNR_ANOVA.mat');
save(filename, 'matlabbatch')
disp('SPM file saved')

disp('Starting jobman')
spm_jobman('run', matlabbatch)

load(batch2)
matlabbatch{1}.spm.stats.fmri_est.spmmat = cellstr(spmmat);
spm_jobman('run', matlabbatch)

end