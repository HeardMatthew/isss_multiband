%% second_level_AUE
% Runs 1 sample t-test and/or displays results 

function second_level_AUE_noacc(subj, study, dd, ss, create, plot)
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

if ~ismember(create, [0 1]) || ~ismember(plot, [0 1])
    error('flags to specify create and plot should be 0 or 1!')
end

%% Parameters
thresh = study.results.thresh;
extent = study.results.extent;

scan = study.scan(ss); 
scanname = scan.runname(1:end-4); 

numsubj = length(subj); 

design = study.design(dd); 

numcon = length(design.con.name); 

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
dir_data   = fullfile(study.path, 'data'); 
dir_second = fullfile(dir_data, 'second_level'); 
dir_dest   = study.results.dest; 

%% Create model
if create
    batch = fullfile(dir_mlb, 'SPM_second_level_1sampleT.mat');
    batch2 = fullfile(dir_mlb, 'SPM_estimate_GLM.mat'); 
    
    %% Name the directory 
    dir_scandesign = fullfile(dir_second, '091120_noacc', [scanname '_MVPA_conn_noacc']); 
    if ~exist(dir_scandesign, 'file'); mkdir(dir_scandesign); end
    
    for cc = 1:numcon
        thiscon = design.con.name{cc}; 
        dir_con = fullfile(dir_scandesign, thiscon);
        if ~exist(dir_con, 'file'); mkdir(dir_con); end
        
        spmmat = fullfile(dir_con, 'SPM.mat'); 
        if exist(spmmat, 'file'); delete(spmmat); end
        
        load(batch)
        matlabbatch{1, 1}.spm.stats.factorial_design.dir = cellstr(dir_con);
        
        %% Load AUE files
        AUEfiles = cell(numsubj, 1);
        
        if numsubj == 14 % all subjects
            for ii = 1:numsubj % each subject
                thissubj = subj(ii);
                disp(['Loading subj ' thissubj.name ' into batch.'])

                dir_subj = fullfile(dir_data, thissubj.name);
                dir_design = fullfile(dir_subj, 'design', [scanname '_' design.name]);
                AUEfiles{ii} = fullfile(dir_design, ['AUE_' thiscon '.nii']);
            end
            
        else
            for ii = 1:numsubj % each subject
                thissubj = subj(ii);
                disp(['Loading subj ' thissubj.name ' into batch.'])

                dir_subj = fullfile(dir_data, thissubj.name);
                dir_design = fullfile(dir_subj, 'design', [scanname '_groupreg_conn']); % alright, just gonna force it. 
                AUEfiles{ii} = fullfile(dir_design, ['AUE_' thiscon '.nii']);
            end
            
        end
        
%     if cc == 3 % drop DM from OR > SR
%         AUEfiles(contains(AUEfiles, 'DM_05Feb18')) = [];
%     end

        matlabbatch{1, 1}.spm.stats.factorial_design.des.t1.scans = AUEfiles;
        
        %% Save the batch
        filename = fullfile(dir_con, ['second_level_' thiscon '.mat']);
        save(filename, 'matlabbatch')
        disp('SPM file saved')

        disp('Starting jobman')
        spm_jobman('run', matlabbatch)
        
        load(batch2)
        matlabbatch{1, 1}.spm.stats.fmri_est.spmmat = cellstr(spmmat);
        spm_jobman('run', matlabbatch)
    end

end

%% get results
if plot
    batch = fullfile(dir_mlb, 'SPM_build_contrast_second_level.mat'); 
    batch2 = fullfile(dir_mlb, 'SPM_report_results.mat'); 
    th = num2str(thresh);
    
    %% Paths
    dir_scandesign = fullfile(dir_second, [scanname '_' design.name]); 
    dir_destscan = fullfile(dir_dest, [scanname '_' design.name]); 
    
    for cc = 1:numcon
        thiscon = design.con.name{cc}; 
        dir_con = fullfile(dir_scandesign, thiscon);
        spmmat = fullfile(dir_con, 'SPM.mat'); 
        
        dir_destcon = fullfile(dir_destscan, thiscon);
        if ~exist(dir_destcon, 'file'); mkdir(dir_destcon); end
        
        %% Evaluate the contrast
        load(batch)
        matlabbatch{1}.spm.stats.con.spmmat = cellstr(spmmat);
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = [scanname '_' thiscon];
        spm_jobman('run', matlabbatch);
        
        %% Report results
        load(batch2)
        matlabbatch{1}.spm.stats.results.spmmat = {spmmat};
        matlabbatch{1}.spm.stats.results.conspec.contrasts = 1;
        matlabbatch{1}.spm.stats.results.conspec.thresh = thresh;
        matlabbatch{1}.spm.stats.results.conspec.extent = extent;
        matlabbatch{1}.spm.stats.results.export(1:2) = []; % png only
        spm_jobman('run', matlabbatch);
        
        %% Save results
        cd(dir_con)
        files = dir('*.png');
        for ii = 1:length(files)
            target = fullfile(pwd, files(ii).name);
            movefile(target, dir_destcon);
        end
        
        cd(dir_batch)
        
    end
    
end


end
