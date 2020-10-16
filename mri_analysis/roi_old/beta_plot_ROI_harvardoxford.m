%%% beta plot
%%% written by ysl 06/27/2014
%%% updated by mjh 11/22/2017
%%% cloned by mjh 06/14/2018
%%% updated by mjh 05/08/2019
% Change ROIs
% Change how plots are formatted: One plot for each condition, featuring
% each scan type
close all; clear all; clc;
tic 
%% Paths and parameters
dir_batch_roi = pwd;
cd ..
isss_multi_params
cd ..
cd rois
dir_rois = pwd;

% disp('loading SPM')
% spm('defaults', 'FMRI');
% spm_jobman('initcfg');
% disp('SPM loaded')

scan(1).name = 'hybrid'; 
scan(1).num_TC_per_run = 180;
scan(1).num_TR = 10;
scan(1).num_runs = 2;

scan(2).name = 'isss'; 
scan(2).num_TC_per_run = 90;
scan(2).num_TR = 5;
scan(2).num_runs = 2;

scan(3).name = 'multi_FIR'; 
scan(3).num_TC_per_run = 248;
scan(3).num_TR = 10;
scan(3).num_runs = 2;

rois_nii = dir('*.nii');
names = {'LIFG pOp & pTri', 'LSTG post', 'RSTG post'}; % 'LIFG pOp', 'LIFG pTri', 
numROI = length(rois_nii);
conname = {'Noise', 'Silence', 'OR', 'SR'};
numCond = 4; %NOI, SIL, ORA, SRA
betas_scans = cell(1, 3); % hybrid, isss, multi
error_scans = cell(1, 3);

for ROI = 1:numROI % For each ROI...
    disp(rois_nii(ROI).name)
    
    for sc = 1:length(scan) % For each scan type... 
        thisscan = scan(sc);
        disp(thisscan.name)
        betas_total = zeros(thisscan.num_TR, numCond, 2, length(subjects));
        % Dimensions [TR, CON, RUN, SUBJ]

        for ss = 1:length(subjects) % For each subject...
            thissubj = subjects{ss};        
            disp(thissubj);
            dir_design = fullfile(dir_data, thissubj, 'design', thisscan.name);
            cd(dir_design);

            for rr = 1:2 % For each run of the hybrid and ISSS conditions ...
                % This code is a solution I figured out for processing multiple
                % runs of MRI data. The design of each GLM I pulled data from
                % had a variety of regressors and this was put in place to fix
                % problems I had for each run. The variable runIdx specifies
                % the first regressor of a run (THE FIRST RUN IS REPRESENTED BY
                % ZERO).
                if rr == 1
                    runIdx = 0;
                    
                    % Load mask
                    this_mask = fullfile(dir_rois, rois_nii(ROI).name);
                    mask_v = spm_vol(this_mask);
                    [beta_mask, xyz_mask] = spm_read_vols(mask_v);
                    
                elseif rr == 2
                    if any(strcmp(thissubj, {'ZG_03Nov17', 'CC_04Jan18'})) 
                    % ZG and CC have no physio = 4 fewer regressors 
                        runIdx = numCond * thisscan.num_TR + 6;
                    else
                        runIdx = numCond * thisscan.num_TR + 10;
                    end
                    
                end
                
                cd(dir_design)
                       
                for cond = 1:numCond % For each condition... (NOI, SIL, ORA, SRA)
                    imgIdx = (cond-1)*thisscan.num_TR + runIdx;
                    % I use imgIdx to begin the following loop through each TR
                    % of the scan. Its value depends on which condition I am
                    % looking at (NOI, SIL, ORA, SRA) and which run (run1,
                    % run2) is being analyzed. 
                    
                    for TR=1:thisscan.num_TR % For each TR... 
                        thisTR = imgIdx + TR;
                        if thisTR < 10
                            thisTR = ['000' num2str(thisTR)];
                        elseif thisTR < 100
                            thisTR = ['00'  num2str(thisTR)];
                        end
                        
                        beta_V = spm_vol(['beta_' thisTR '.nii']);
                        [beta_all, XYZ] = spm_read_vols(beta_V);  % loads all voxels
                        beta_vec = beta_all(logical(beta_mask)); % gets all voxels within mask
                        beta_vec = beta_vec(~isnan(beta_vec)); % removes voxels from mask outside of betas

                        betas_total(TR, cond, rr, ss)=mean(beta_vec);
                    end
                    
                end
                
            end
            
        end

        %% Calculating mean and std dev
        betas_total_across_run = mean(betas_total, 3); % average across subjects
        betas_total_across_run_sub = mean(betas_total_across_run, 4); % average across runs
        std_dev = zeros(thisscan.num_TR, numCond);
        for cond = 1:numCond
            for TR = 1:thisscan.num_TR
                std_dev(TR,cond) = std(betas_total_across_run(TR, cond, 1, :));
            end
            
        end

        std_error = std_dev/sqrt(length(subjects));
        
        betas_scans{sc} = betas_total_across_run_sub;
        error_scans{sc} = std_error; 

    end % all scans
    
    %% Plotting results
    % I want to make plots for OR, SR, NOI, and SIL across scans
    
    cd(dir_rois)
    for cc = 1:numCond
        figure;
        hold on;
        
        h1 = plot(betas_scans{1}(:, cc), 'r'); % hybrid
        errorbar(betas_scans{1}(:, cc), error_scans{1}(:, cc), 'r');

        h2 = plot([2:2:10], betas_scans{2}(:, cc), 'k'); % isss
        errorbar([2:2:10], betas_scans{2}(:, cc), error_scans{2}(:, cc), 'k');

        h3 = plot(betas_scans{3}(:, cc), 'b'); %ORA
        errorbar(betas_scans{3}(:, cc), error_scans{3}(:, cc), 'b');
        xlim([0 11])
%         ylim([-.05 0.16]);
        
        name_title = [names{ROI} ' ' conname{cc}];
        name_file = regexp(names{ROI}, '\s+', 'split'); name_file = strjoin([name_file, conname{cc}, 'HarvardOxford'], '_');

        title(name_title, 'FontWeight', 'bold', 'Fontsize', 14)

        hleg1=legend([h1, h2, h3],'MISSS','ISSS','CM','Location','NorthEast');

        xlabel('Post Stimulus Window (sec.)', 'Fontsize', 14);
        ylabel('Beta Estimates', 'Fontsize', 14);
        
         %save in png format 
         saveas(gcf, name_file, 'png');

         %save in eps format 
%          saveas(gcf, ['beta_plot_' Region{ROI}], 'eps')

    end

end % all ROIs

toc
