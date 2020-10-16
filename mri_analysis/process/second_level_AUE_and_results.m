%% second_level_AUE

function second_level_AUE(subj, 
create = 1;
plot = 0;

thresh = 0.001;
extent = 20;

warning off;
spm('defaults', 'FMRI');
spm_jobman('initcfg');

dir_second = fullfile(dir_data, 'second_level');

%% Create model
if create
    for ss = scantypes % each scan type
        for cc = 9%1:length(cond) % each contrast
            cd(dir_mlb)
            load('SPM_second_level_1sampleT.mat')

%             matlabbatch{1, 1}.spm.stats.factorial_design.dir = cellstr(fullfile(dir_second, [scan{ss} '_dropbetas'], contrasts{cc}));
            matlabbatch{1, 1}.spm.stats.factorial_design.dir = cellstr(fullfile(dir_second, [scan{ss}], contrasts{cc}));

            AUEfiles = cell(numsubj, 1);

            for ii = 1:numsubj % each subject
                thissubj = subjects{ii};
                disp(['Loading subj ' thissubj ' into batch.'])

                dir_subj = fullfile(dir_data, thissubj);
                dir_design = fullfile(dir_subj, 'design');

%                 AUEfiles{ii} = fullfile(dir_design, scan{ss}, ['AUE_' cond{cc} '_dropbetas.nii']);
                AUEfiles{ii} = fullfile(dir_design, scan{ss}, ['AUE_' cond{cc} '.nii']);

            end
            
            if cc == 3 % drop DM from OR > SR
                AUEfiles(contains(AUEfiles, 'DM_05Feb18')) = [];
            end

            matlabbatch{1, 1}.spm.stats.factorial_design.des.t1.scans = AUEfiles;
            disp(['All subjects are in the batch for ' contrasts{cc} '.'])

            mkdir(fullfile(dir_second, scan{ss}, contrasts{cc}))
            filename = fullfile(dir_second, scan{ss}, contrasts{cc}, ['second_level_' contrasts{cc} '_nowrong']);
            save(filename, 'matlabbatch')
            disp('SPM file saved')

            disp('Starting jobman')
            spm_jobman('run', matlabbatch)

            load('GLM_estimate.mat')
%             matlabbatch{1, 1}.spm.stats.fmri_est.spmmat = cellstr(fullfile(dir_second, [scan{ss} '_dropbetas'], contrasts{cc}, 'SPM.mat'));
            matlabbatch{1, 1}.spm.stats.fmri_est.spmmat = cellstr(fullfile(dir_second, [scan{ss}], contrasts{cc}, 'SPM.mat'));
            spm_jobman('run', matlabbatch)

        end

    end
    
end

%% get results
if plot
    th = num2str(thresh);
    % th = th(3:end);
    dir_dest = { ... 
        ['C:\Users\heard.49\Documents\MATLAB\fmri\isss_multiband_preprocessed_04\docs\14subj_07232019\hybrid_' th '_' num2str(extent)] ...
        ['C:\Users\heard.49\Documents\MATLAB\fmri\isss_multiband_preprocessed_04\docs\14subj_07232019\isss_' th '_' num2str(extent)] ...
        ['C:\Users\heard.49\Documents\MATLAB\fmri\isss_multiband_preprocessed_04\docs\14subj_07232019\multi_FIR_' th '_' num2str(extent)] ...
        };
    for ii = scantypes
        mkdir(dir_dest{ii})
    end

    filename = regexp(date, '-', 'split');
    % filename_ps = ['spm_' filename{3} filename{2} filename{1} '*.ps'];
    % filename_pdf = ['spm_' filename{3} filename{2} filename{1} '*.pdf'];
    filename_png = ['spm_' filename{3} filename{2} filename{1} '*.png'];

    for ss = scantypes % each scan type

        for cc = 1:length(cond) % each contrast
            cd(dir_process)
            load('SPM_build_contrasts_secondlevel.mat')
            dir_design = fullfile(dir_second, scan{ss}, contrasts{cc});
            this_spm_mat = fullfile(dir_design, 'SPM.mat'); 
            matlabbatch{1}.spm.stats.con.spmmat = cellstr(this_spm_mat);
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = [scan{ss} '_' contrasts{cc}];
            spm_jobman('run', matlabbatch);

            cd(dir_process)
            load SPM_report_results.mat
            cd(dir_design)
            matlabbatch{1}.spm.stats.results.spmmat = {this_spm_mat};
            matlabbatch{1}.spm.stats.results.conspec.contrasts = 1;
            matlabbatch{1}.spm.stats.results.conspec.thresh = thresh;
            matlabbatch{1}.spm.stats.results.conspec.extent = extent;

            spm_jobman('run', matlabbatch);

    %         target = fullfile(dir_design, filename_ps);
    %         dest = fullfile(dir_dest{ss}, [scan{ss} '_' cond{cc} '.ps']);
    %         movefile(target, dest);

            files = dir('*.png');
            for ii = 1:length(files)
                target = fullfile(dir_design, files(ii).name);
                dest = fullfile(dir_dest{ss}, [scan{ss} '_' cond{cc} '_' num2str(ii) '.png']);
                movefile(target, dest);
            end

    %         target = fullfile(dir_design, filename_pdf);
    %         dest = fullfile(dir_dest{ss}, [scan{ss} '_' cond{cc} '.pdf']);
    %         movefile(target, dest);

        end
        
    end

end
