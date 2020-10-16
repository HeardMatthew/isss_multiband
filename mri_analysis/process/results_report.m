%% results_report


function results_report(varargin)
%% Check input
subj = varargin{1}; 
study = varargin{2}; 
dd    = varargin{3}; 
ss    = varargin{4}; 

if ~isstruct(subj) || length(subj) ~= 1
    error('Input ("subj") where subj is a SINGLE struct')
end

if ~isstruct(study)
    error('Input ("study") which has experiment info')
end

if ~isnumeric(dd)
    error('Input ("dd") which specifies which designs (#) should be used')
end

if ~isnumeric(ss)
    error('Input ("ss") which specifies scan!')
end

%% Set up SPM
disp('Loading SPM...')
close all; 
spm('defaults', 'FMRI');
spm_jobman('initcfg');
spm_figure('GetWin','Graphics'); % Thanks Guillaume
disp('Done!')

%% Parameters
scan = study.scan(ss); 
scanname = strsplit(scan.runname, '_'); scanname = scanname{1}; 

design = study.design(dd); 
con = design.con.name; 

thresh = study.results.thresh;
extent = study.results.extent;

%% Pathing
cd ..
dir_batch = pwd; 
dir_dest  = fullfile(study.path, 'docs', '041620_single_subj_conditions_p001_k20'); 
dir_mlb   = fullfile(dir_batch, 'matlabbatch'); 

dir_subj  = fullfile(study.path, 'data', subj.name); 
dir_thisdesign  = fullfile(dir_subj, 'design', [scanname '_' design.name]); 
spmmat = fullfile(dir_thisdesign, 'SPM.mat');

filename = regexp(date, '-', 'split');
filename = ['spm_' filename{3} filename{2} filename{1} '.png'];
batch = fullfile(dir_mlb, 'SPM_report_results.mat'); 

for cc = 1:length(con)
    dir_thisdest = fullfile(dir_dest, [scanname '_' design.name], con{cc}); 
    if ~exist(dir_thisdest, 'file')
        mkdir(dir_thisdest)
    end
    
    load(batch)
    matlabbatch{1}.spm.stats.results.export(1:2) = []; % only output PNG
    matlabbatch{1}.spm.stats.results.spmmat = {spmmat};
    matlabbatch{1}.spm.stats.results.conspec.contrasts = cc;
    matlabbatch{1}.spm.stats.results.conspec.thresh = thresh;
    matlabbatch{1}.spm.stats.results.conspec.extent = extent;

    spm_jobman('run', matlabbatch);
    
    

    target = dir([dir_thisdesign filesep '*.png']);
    for ii = 1:length(target)
        thistarget = fullfile(target(ii).folder, target(ii).name); 
        thisdest = fullfile(dir_dest, [scanname '_' design.name], con{cc}, ...
            [subj.name '_' con{cc} '_' num2str(ii) '.png']);
        movefile(thistarget, thisdest);
    end

end
