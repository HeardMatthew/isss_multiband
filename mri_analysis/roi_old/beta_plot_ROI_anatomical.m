%%% beta plot
%%% written by ysl 06/27/2014
%%% updated by mjh 11/22/2017
%%% cloned by mjh 06/14/2018
close all; clear all; clc;
tic 
%% Paths and parameters
dir_batch_roi = pwd;
cd ..
isss_multi_params
cd ..
addpath spm12

spm('defaults', 'FMRI');
spm_jobman('initcfg');

scan(1).name = 'hybrid'; 
scan(1).num_TC_per_run = 180;
scan(1).num_TR = 10;
scan(1).num_runs = 2;

scan(2).name = 'isss'; 
scan(2).num_TC_per_run = 90;
scan(2).num_TR = 5;
scan(2).num_runs = 2;

cd(dir_batch_roi)
rois_nii = dir('*.nii');
if length(rois_nii) ~= 3
    error('Check nii files')
end

numROI = length(rois_nii);
numCond = 4; %NOI, SIL, ORA, SRA

for ROI = 1:numROI % For each ROI...
    % As of 06/18/18 I know the first ROI correctly generates. So let's
    % troubleshoot number 2. 
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
                    this_mask = fullfile(dir_batch_roi, 'mask_coreg', [thissubj '_' thisscan.name '_r' rois_nii(ROI).name]);
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

        %% Plotting results
        figure;
        hold on;

        h1=plot(betas_total_across_run_sub(:,1),'b'); %NOI
        errorbar(betas_total_across_run_sub(:,1),std_error(:,1), 'b');

        h2= plot(betas_total_across_run_sub(:,2),'b','LineStyle','--'); %SIL
        errorbar(betas_total_across_run_sub(:,2),std_error(:,2), 'b','LineStyle','--');

        h3=plot(betas_total_across_run_sub(:,3),'r'); %ORA
        errorbar(betas_total_across_run_sub(:,3),std_error(:,3), 'r');

        h4=plot(betas_total_across_run_sub(:,4),'r','LineStyle','--'); %SRA
        errorbar(betas_total_across_run_sub(:,4),std_error(:,4), 'r','LineStyle','--');

    % hold off; 
    % axis off;

       %%for getting EPS without any labeling/ axis comments out below  
        name = regexp(rois_nii(ROI).name, '_', 'split');
        name_title = strjoin(name(2:end-1), '\\_');
        name_file = strjoin(name(2:end-1), '_');
        
        set(gca,'XTick',1:thisscan.num_TR)

        label_x = cell(1, thisscan.num_TR);
        if strcmp(thisscan.name, 'hybrid')
            for ii = 1:thisscan.num_TR
                label_x{ii} = num2str(ii);
            end
        else
            for ii = 1:thisscan.num_TR
                label_x{ii} = num2str(ii*2);
            end
        end
        
        set(gca,'XTickLabel', label_x)

        title(name_title, 'FontWeight', 'bold', 'Fontsize', 14)

        hleg1=legend([h3,h4,h1,h2],'ORA','SRA','NOI','SIL','Location','NorthEast');

        xlabel('Post Stimulus Window (sec.)', 'Fontsize', 14);
        ylabel('Beta Estimates', 'Fontsize', 14);

          cd (dir_batch_roi);

         %save in png format 
         saveas(gcf,['beta_plot_anatomy_' thisscan.name '_' name_file], 'png');

         %save in eps format 
    %      saveas(gcf, ['beta_plot_' Region{ROI}], 'eps')

    %         clf;
    %     

        clear std_dev; 
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

    end % all scans

end % all ROIs

toc
