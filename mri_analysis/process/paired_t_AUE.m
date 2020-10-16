%% second_level_AUE
% Runs paired t-test and/or displays results 

function paired_t_AUE(subj, study, dd, pair, create, plot)
%% check inputs
if length(subj) < 2
    error('Submit ALL subjects!')
end

if ~isstruct(subj) || ~isstruct(study)
    error('subj and study are both struct!')
end

if ~isnumeric(dd)
    error('dd specifies number of design and scan!')
end

if ~ischar(pair)
    error('pair specifies which paired t test to run!')
end

if ~ismember(create, [0 1]) || ~ismember(plot, [0 1])
    error('flags to specify create and plot should be 0 or 1!')
end

%% Parameters
thresh = study.results.thresh;
extent = study.results.extent;

ps = strsplit(pair, '_'); pair1 = ps{1}; pair2 = ps{2}; 
idx1 = contains({study.scan(:).runname}, pair1); 
idx2 = contains({study.scan(:).runname}, pair2); 
scan(1) = study.scan(idx1); scan(2) = study.scan(idx2); 

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
dir_mlb    = fullfile(dir_batch,  'matlabbatch'); 
dir_data   = fullfile(study.path, 'data'); 
dir_second = fullfile(dir_data,   'second_level'); 
dir_dest   = study.results.dest; 

%% Create model
if create
    batch = fullfile(dir_mlb, 'SPM_second_level_pairedT.mat');
    batch2 = fullfile(dir_mlb, 'SPM_estimate_GLM.mat'); 
    
    %% Name the directory 
    dir_design = fullfile(dir_second, [pair '_' design.name]); 
    if ~exist(dir_design, 'file'); mkdir(dir_design); end
    
    for cc = 1:numcon
        thiscon = design.con.name{cc}; 
        dir_con = fullfile(dir_design, thiscon);
        if ~exist(dir_con, 'file'); mkdir(dir_con); end
        
        spmmat = fullfile(dir_con, 'SPM.mat'); 
        if exist(spmmat, 'file'); delete(spmmat); end
        
        load(batch)
        matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(dir_con);
        
        %% Load AUE files
        for ii = 1:numsubj % each subject
            thissubj = subj(ii);
            disp(['Loading subj ' thissubj.name ' into batch.'])    
            dir_thissubj = fullfile(dir_data, thissubj.name);
            AUEfiles = cell(2, 1);
            
            for ss = 1:2 % each pair/scan
                name = strsplit(scan(ss).runname, '_'); name = name{1}; 
                dir_thisdesign = fullfile(dir_thissubj, 'design', [name '_' design.name]); 
                AUEfiles{ss} = fullfile(dir_thisdesign, ['AUE_' thiscon '.nii']);
                
            end
            
            matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(ii).scans = AUEfiles; 
        end
        
%     if cc == 3 % drop DM from OR > SR
%         AUEfiles(contains(AUEfiles, 'DM_05Feb18')) = [];
%     end
        
        %% Remove accuracy covariate from template
        matlabbatch{1}.spm.stats.factorial_design.cov = []; 
%         if design.groupreg
%             matlabbatch{1}.spm.stats.factorial_design.cov = struct( ...
%                 'c', [], 'cname', 'group_accuracy', 'iCFI', 1, 'iCC', 1); 
% 
%             data = fullfile(dir_batch, ['group_level_acc_' scanname '.mat']);
%             load(data)
%             matlabbatch{1}.spm.stats.factorial_design.cov.c = group_level_acc; 
%             matlabbatch{1}.spm.stats.factorial_design.cov.cname = 'group_level_acc'; 
%         end
        
        %% Save the batch
        filename = fullfile(dir_con, ['second_level_ptt_' thiscon '.mat']);
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
    batch  = fullfile(dir_mlb, 'SPM_build_contrast_second_level.mat'); 
    batch2 = fullfile(dir_mlb, 'SPM_report_results.mat'); 
    th = num2str(thresh);
    
    %% Paths
    dir_design = fullfile(dir_second, [pair '_' design.name]); 
    dir_destscan = fullfile(dir_dest, [pair '_' design.name]); 
    
    for cc = 1:numcon
        thiscon = design.con.name{cc}; 
        dir_con = fullfile(dir_design, thiscon);
        spmmat = fullfile(dir_con, 'SPM.mat'); 
        
        dir_destcon = fullfile(dir_destscan, thiscon);
        if ~exist(dir_destcon, 'file'); mkdir(dir_destcon); end
        
        %% Evaluate the contrast
        load(batch)
        name1 = pair; 
        name2 = strsplit(pair, '_'); name2 = [name2{2}, '_', name2{1}]; 
        matlabbatch{1}.spm.stats.con.spmmat = cellstr(spmmat);
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = [name1 '_' thiscon];
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 -1]; 
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'repl'; 
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = [name2 '_' thiscon];
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [-1 1]; 
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'repl'; 
        spm_jobman('run', matlabbatch);
        
        %% Report results (do it manually for now)
%         load(batch2)
%         matlabbatch{1}.spm.stats.results.spmmat = {spmmat};
%         matlabbatch{1}.spm.stats.results.conspec.contrasts = 1;
%         matlabbatch{1}.spm.stats.results.conspec.thresh = thresh;
%         matlabbatch{1}.spm.stats.results.conspec.extent = extent;
%         matlabbatch{1}.spm.stats.results.export(1:2) = []; % png only
%         spm_jobman('run', matlabbatch);
%         
%         %% Save results
%         cd(dir_con)
%         files = dir('*.png');
%         for ii = 1:length(files)
%             target = fullfile(pwd, files(ii).name);
%             movefile(target, dir_destcon);
%         end
%         
%         cd(dir_batch)
        
    end
    
end


end
