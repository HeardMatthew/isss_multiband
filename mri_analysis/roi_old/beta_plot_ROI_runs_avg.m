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
betas_scans = cell(1, 3); % Each scan type
error_scans = cell(1, 3);

for ROI = 1:numROI % For each ROI...
    disp(rois_nii(ROI).name)
    
    for sc = 1:length(scan) % For each scan type... 
        thisscan = scan(sc);
        disp(thisscan.name)
        betas_total = zeros(thisscan.num_TR, numCond, 2, length(subjects));
        error_total = zeros(thisscan.num_TR, numCond, 2, length(subjects));
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
                        error_total(TR, cond, rr, ss)=std(beta_vec)/sqrt(length(beta_vec));
                    end
                    
                end
                
            end     
            
        end
        
        betas_across_subj = mean(betas_total, 4);
        error_across_subj = std(betas_total, [], 4)/sqrt(size(betas_total, 4));
        
        betas_scans{sc} = betas_across_subj;
        error_scans{sc} = error_across_subj;
    end % all scans
    
    %% Plotting results
    % Each subject gets a plot, with eight lines, four for each run
    
    cd(dir_rois)
    x = {1:10, 2:2:10, 1:10};
    
    for sc = 1:3 % each scan
        figure; hold on;
        h1 = plot(x{sc}, betas_scans{sc}(:, 1, 1), 'b'); % noise run 1
        h2 = plot(x{sc}, betas_scans{sc}(:, 2, 1), 'k'); % silence run 1
        h3 = plot(x{sc}, betas_scans{sc}(:, 3, 1), 'r'); % OR run 1
        h4 = plot(x{sc}, betas_scans{sc}(:, 4, 1), 'm'); % SR run 1
        h5 = plot(x{sc}, betas_scans{sc}(:, 1, 2), 'b--'); % noise run 2
        h6 = plot(x{sc}, betas_scans{sc}(:, 2, 2), 'k--'); % silence run 2
        h7 = plot(x{sc}, betas_scans{sc}(:, 3, 2), 'r--'); % OR run 2
        h8 = plot(x{sc}, betas_scans{sc}(:, 4, 2), 'm--'); % SR run 2
        
        xlim([0 11])
        name_title = [scan(sc).name ' ' names{ROI} ' Average'];
        title(name_title, 'FontWeight', 'bold', 'Fontsize', 14)
        xlabel('Post Stimulus Window (sec.)', 'Fontsize', 14);
        ylabel('Beta Estimates', 'Fontsize', 14);
        hleg1=legend([h1, h2, h7, h8],'Noise1', 'Silent1', 'OR2', 'SR2', 'Location','NorthWest');
        
        name_file = regexp(names{ROI}, '\s+', 'split'); name_file = strjoin([scan(sc).name, name_file, 'avg', 'HarvardOxford', 'nobars'], '_');
        saveas(gcf, name_file, 'png');

        errorbar(x{sc}, betas_scans{sc}(:, 1, 1), error_scans{sc}(:, 1, 1), 'b');
        errorbar(x{sc}, betas_scans{sc}(:, 2, 1), error_scans{sc}(:, 2, 1), 'k');
        errorbar(x{sc}, betas_scans{sc}(:, 3, 1), error_scans{sc}(:, 3, 1), 'r');
        errorbar(x{sc}, betas_scans{sc}(:, 4, 1), error_scans{sc}(:, 4, 1), 'm');
        errorbar(x{sc}, betas_scans{sc}(:, 1, 2), error_scans{sc}(:, 1, 2), 'b--');
        errorbar(x{sc}, betas_scans{sc}(:, 2, 2), error_scans{sc}(:, 2, 2), 'k--');
        errorbar(x{sc}, betas_scans{sc}(:, 3, 2), error_scans{sc}(:, 3, 2), 'r--');
        errorbar(x{sc}, betas_scans{sc}(:, 4, 2), error_scans{sc}(:, 4, 2), 'm--');
        
        hleg1=legend([h1, h2, h7, h8],'Noise1', 'Silent1', 'OR2', 'SR2', 'Location','NorthWest');
        name_file = regexp(names{ROI}, '\s+', 'split'); name_file = strjoin([scan(sc).name, name_file, 'avg', 'HarvardOxford', 'bars'], '_');
        saveas(gcf, name_file, 'png');
         %save in eps format 
%          saveas(gcf, ['beta_plot_' Region{ROI}], 'eps')
        
    end
    
end % all ROIs

toc
