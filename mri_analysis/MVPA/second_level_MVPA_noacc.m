%% second_level_MVPA
% Runs 1 sample t-test and/or displays results 
% 
% MM/DD/YY -- CHANGELOG
% 05/01/20 -- Log started, forked for MVPA

function second_level_MVPA_noacc(subj, study, dd, ss, classifier)
%% check inputs
if length(subj) < 2
    error('Submit ALL subjects!')
end

if ~isstruct(subj) || ~isstruct(study)
    error('subj and study are both struct!')
end

if ~isnumeric(dd) || ~isnumeric(ss)
    error('dd and ss specify number of design and scan!')
end

if ~ischar(classifier)
    error('please specify which classifier!')
end

%% Parameters
scan = study.scan(ss); 
scanname = scan.runname(1:end-4); 

numsubj = length(subj); 
design = study.design(dd); 
numcon = length(design.con.name); 

rad = study.mvpa.radius; 

%% Pathing
cd ..
dir_batch  = pwd; 
dir_mlb    = fullfile(dir_batch, 'matlabbatch'); 
dir_data   = fullfile(study.path, 'data'); 
dir_second = fullfile(dir_data, 'second_level'); 

%% Create model
batch = fullfile(dir_mlb, 'SPM_second_level_1sampleT.mat');
batch2 = fullfile(dir_mlb, 'SPM_estimate_GLM.mat'); 

%% Name the directory 
dir_scanclass = fullfile(dir_second, '091120_noacc', ['MVPA_' scanname '_conn_noacc']); 
if ~exist(dir_scanclass, 'file'); mkdir(dir_scanclass); end

spmmat = fullfile(dir_scanclass, 'SPM.mat'); 
if exist(spmmat, 'file'); delete(spmmat); end

load(batch)
matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(dir_scanclass);

%% Load AUE files
classfiles = cell(numsubj, 1);

for ii = 1:numsubj % each subject
    thissubj = subj(ii);
    disp(['Loading subj ' thissubj.name ' into batch.'])

    dir_subj = fullfile(dir_data, thissubj.name);
    dir_mvpa = fullfile(dir_subj, 'MVPA'); 
    fname = strjoin({thissubj.name, scanname, design.name, 'beta', classifier, ['rad' num2str(rad)]}, '_');
    classfiles{ii} = fullfile(dir_mvpa, [fname '.nii']);
end

matlabbatch{1, 1}.spm.stats.factorial_design.des.t1.scans = classfiles;

%% Load accuracy covariate
% if ~contains(classifier, 'SIL') 
%     if design.groupreg
%         matlabbatch{1}.spm.stats.factorial_design.cov = struct( ...
%             'c', [], 'cname', 'group_accuracy', 'iCFI', 1, 'iCC', 1); 
% 
%         data = fullfile(dir_batch, ['group_level_acc_' scanname '.mat']);
%         load(data)
%         %%% I HAD TO HARD CODE SUBJ 1 OUT OF THE REGRESSOR!
%         if contains(classifier, 'OR_SR')
%             matlabbatch{1}.spm.stats.factorial_design.cov.c = group_level_OR_SR(2:end); 
%         elseif contains(classifier, 'OR_NOI')
%             matlabbatch{1}.spm.stats.factorial_design.cov.c = group_level_OR(2:end); 
%         elseif contains(classifier, 'SR_NOI')
%             matlabbatch{1}.spm.stats.factorial_design.cov.c = group_level_SR(2:end); 
%         else % LNG classifier
%             matlabbatch{1}.spm.stats.factorial_design.cov.c = group_level_acc(2:end); 
%         end
% 
%     end
% 
% end

%% Save the batch
filename = fullfile(dir_scanclass, ['second_level_' classifier '.mat']);
save(filename, 'matlabbatch')
disp('SPM file saved')

%% Prepare SPM
disp('Loading SPM...')
close all; 
spm('defaults', 'FMRI');
spm_jobman('initcfg');
spm_figure('GetWin','Graphics'); % Thanks Guillaume
disp('Done!')

%% Run jobs
disp('Starting jobman')
spm_jobman('run', matlabbatch)

load(batch2)
matlabbatch{1, 1}.spm.stats.fmri_est.spmmat = cellstr(spmmat);
spm_jobman('run', matlabbatch)

end
