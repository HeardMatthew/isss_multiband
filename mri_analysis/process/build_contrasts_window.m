%% build_contrasts
% Builds contrasts for 1st level analysis for one subject

% CHANGELOG 
% 11/15/17  Initialized file (hybrid_isss)
% 02/13/20  Forked for YA_FST, made universal
% 03/02/20  Dropping first TR from each event
% 03/24/20  Moving window analysis

function build_contrasts_window(subj, study, design)
%% Check input
if ~isstruct(subj) || length(subj) ~= 1
    error('Input ("subj") where subj is a SINGLE struct')
end

if ~isstruct(study)
    error('Input ("study") which has experiment info')
end

if ~isnumeric(design)
    error('Input ("dd") which specifies which designs (#) should be used')
end

%% Prepare SPM
disp('Loading SPM...')
close all; 
spm('defaults', 'FMRI');
spm_jobman('initcfg');
disp('Done!')

%% Pathing
cd ..
dir_batch = pwd; 
dir_mlb = fullfile(dir_batch, 'matlabbatch'); 
dir_subj = fullfile(study.path, 'data', subj.name); 
dir_design_root = fullfile(dir_subj, 'design');

for dd = design
    thisdesign = study.design(dd); 
    numcond = length(thisdesign.cond); 
    
    %% Load matlabbatch, add SPM.mat
    batch = fullfile(dir_mlb, 'SPM_build_contrasts.mat'); 
    load(batch)
    
    spmmat = fullfile(dir_design_root, thisdesign.name, 'SPM.mat');
    matlabbatch{1}.spm.stats.con.spmmat = cellstr(spmmat);
    
    %% For each contrast...
    idx = 1; 
    for cc = 1:length(thisdesign.con.name)
        for tt = 1:(study.scan.epis-1)
            index = [zeros(1, tt-1), 1, zeros(1, study.scan.epis-1-tt)]; 
            index = repmat(index, [1, numcond]); 
            thisvec = repelem(thisdesign.con.vec(cc, :), study.scan.epis-1); 
            thisvec = thisvec .* index; 
            
            matlabbatch{1}.spm.stats.con.consess{idx} = ... 
                matlabbatch{1}.spm.stats.con.consess{1}; % Preallocate
        
            matlabbatch{1}.spm.stats.con.consess{idx}.tcon.name = ...
                [thisdesign.con.name{cc} '_TR' num2str(tt+1)];
            matlabbatch{1}.spm.stats.con.consess{idx}.tcon.weights = ...
                thisvec;     
            matlabbatch{1}.spm.stats.con.consess{idx}.tcon.sessrep = ...
                'repl'; 
            
            idx = idx + 1; 
        end
        
    end
    
    if size(matlabbatch{1}.spm.stats.con.consess, 2) > idx
        matlabbatch{1}.spm.stats.con.consess(idx+1:end) = [];
    end
    
    matlabbatch{1}.spm.stats.con.delete = 0; 
    
    %% Save files
    filename = fullfile(dir_subj, 'batch', ['build_contrasts_window_' thisdesign.name '.mat']);
    save(filename, 'matlabbatch');
    
    % Run job
    spm_jobman('run', matlabbatch);
end
