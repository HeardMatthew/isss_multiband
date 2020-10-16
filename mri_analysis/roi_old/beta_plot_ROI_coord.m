%%% beta plot
%%% written by ysl 06/27/2014
%%% updated by mjh 11/22/2017
close all; clear all; clc;
tic

%% Paths and parameters
[MNI_XYZ,Region] = xlsread([pwd filesep 'multi_roi_list_12subj.xls']);
Region = Region((~cellfun(@isempty, Region)));

dir_batch_roi = pwd;
cd ..
isss_multi_params

scan(1).name = 'hybrid'; 
scan(1).num_TC_per_run = 180;
scan(1).num_TR = 10;
scan(1).num_runs = 2;

scan(2).name = 'isss'; 
scan(2).num_TC_per_run = 90;
scan(2).num_TR = 5;
scan(2).num_runs = 2;

% This is a helper function I run to quickly generate information regarding
% my experiment (e.g. subject list, order of scanning protocols). In
% particular, this script uses data about which subject is being processed.
% Anytime you see code referring to a cell called "subjects", you can
% probably take it out. 

numROI = size(MNI_XYZ,1);
numCond = 4; %NOI, SIL, ORA, SRA

for ROI=1:numROI % For each ROI...
    disp(Region{ROI})
    %% Loading the brain data
    
    this_vox = MNI_XYZ(ROI,:)';
    this_region = Region{ROI};  
    
    for sc = 1:length(scan)
        thisscan = scan(sc);
        disp(thisscan.name)
        numTRs = thisscan.num_TR; 
        betas_total = zeros(numTRs, numCond, 2, length(subjects));
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
                disp(['Run: ', num2str(rr)])
                if rr == 1
                    runIdx = 0;    
                elseif rr == 2
                    if any(strcmp(thissubj, {'ZG_03Nov17', 'CC_04Jan18'})) 
                    % ZG and CC have no physio = 4 fewer regressors 
                        runIdx = numCond * numTRs + 6;
                    else
                        runIdx = numCond * numTRs + 10;
                    end
                    
                end

                for cond = 1:numCond % For each condition... (NOI, SIL, ORA, SRA)
                    imgIdx = (cond-1)*numTRs + runIdx;
                    % I use imgIdx to begin the following loop through each TR
                    % of the scan. Its value depends on which condition I am
                    % looking at (NOI, SIL, ORA, SRA) and which run (run1,
                    % run2) is being analyzed. 
                    
                    for TR=1:numTRs % For each TR... 
                        thisTR = imgIdx + TR;
                        if thisTR < 10
                            thisTR = ['000' num2str(thisTR)];
                        elseif thisTR < 100
                            thisTR = ['00'  num2str(thisTR)];
                        end

                        Vol=spm_vol(['beta_' thisTR '.nii']);
                        [beta, XYZ] = spm_read_vols(Vol);
                        
                        this_vox_idx = find(all(XYZ == this_vox));
                        betas_total(TR, cond, rr, ss) = beta(this_vox_idx);

                    end
                    
                end
                
            end
            
        end

        %% Calculating mean and std dev
        betas_total_across_run = mean(betas_total, 3); % average across subjects
        betas_total_across_run_sub = mean(betas_total_across_run, 4);
        std_dev = zeros(numTRs, numCond);
        for cond = 1:numCond
            for TR = 1:numTRs
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
        name = regexp(Region{ROI}, '_', 'split');
        name_title = strjoin([thisscan.name, name], '\\_');
        name_file = strjoin([thisscan.name, name], '_');
        
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

        xlabel('Post Stimulus (sec.)', 'Fontsize', 14);
        ylabel('Beta Estimates', 'Fontsize', 14);

          cd(dir_batch_roi);

         %save in png format 
         saveas(gcf,['beta_plot_coordinate_' name_file], 'png');

         %save in eps format 
    %      saveas(gcf, ['beta_plot_' Region{ROI}], 'eps')

    %         clf;
    %     

        clear std_dev; 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    end
    
end

toc
