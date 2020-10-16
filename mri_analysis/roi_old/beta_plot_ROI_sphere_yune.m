%%% beta plot
%%% written by ysl 06/27/2014
%%% updated by mjh 11/22/2017
%%% updated to sphere by mjh 06/14/2018
% Uses same sphere radius as MVPA analysis. 

close all; clc;
tic

%% Paths and parameters
dir_roi = 'C:\Users\heard.49\Documents\MATLAB\fmri\isss_multiband_preprocessed_batch_03\rois\yune_coord';
cd(dir_roi);
files_roi = dir('*.nii');
ROI = cell(length(files_roi), 1);
ROI_name = cell(length(files_roi), 1);
for ii = 1:length(files_roi)
    V = spm_vol(files_roi(ii).name);
    ROI{ii} = spm_read_vols(V);    
    ROI_name{ii} = fullfile(files_roi(ii).folder, files_roi(ii).name);
end

dir_docs = 'C:\Users\heard.49\Documents\MATLAB\fmri\isss_multi_docs\11subj_nowrong_beta_plot_sphere_yune';

scan(1).name = 'hybrid_nowrong'; 
scan(1).num_TC_per_run = 180;
scan(1).num_TR = 10;
scan(1).num_runs = 2;

scan(2).name = 'isss_nowrong'; 
scan(2).num_TC_per_run = 90;
scan(2).num_TR = 5;
scan(2).num_runs = 2;

scan(3).name = 'multi_FIR_nowrong';
scan(3).num_TC_per_run = 180;
scan(3).num_TR = 10;
scan(3).num_runs = 2;

numROI = length(ROI);
numCond = 4; %NOI, SIL, ORA, SRA

Region = {'Left IFG', 'Left pSTG'};

for ss = 1:length(subjects)
    for sc = 1:length(scan)
        
        thissubj = subjects{ss};
        dir_subj = fullfile(dir_data, thissubj);
        dir_design = fullfile(dir_subj, 'design');

        cd(fullfile(dir_design, scan(sc).name));
        
        V = spm_vol('mask.nii');
        thismask = spm_read_vols(V);
        if ss == 1 && sc == 1
            mask_master = thismask;
        end
        
        mask_master = mask_master & thismask;
        
    end
    
end

for rr = 1:numROI % For each ROI...
    thisroi = logical(ROI{rr});
    thisroiname = ROI_name{rr};

    for sc = 1:length(scan)
        disp(scan(sc).name)
        numTRs = scan(sc).num_TR; 
        betas_total = zeros(numTRs, numCond, 2, length(subjects));
        % Dimensions [TR, CON, RUN, SUBJ]
        
        for ss = 1:length(subjects) % For each subject...
            thissubj = subjects{ss};        
            disp(thissubj);
%                       
            %% Get data
            
            for rn = 1:2 % For each run of the hybrid and ISSS conditions ...
                % This code is a solution I figured out for processing multiple
                % runs of MRI data. The design of each GLM I pulled data from
                % had a variety of regressors and this was put in place to fix
                % problems I had for each run. The variable runIdx specifies
                % the first regressor of a run (THE FIRST RUN IS REPRESENTED BY
                % ZERO).
                disp(['Run: ', num2str(rn)])
                if rn == 1
                    runIdx = 0;    
                elseif rn == 2
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
                        
                        cd(fullfile(dir_design, scan(sc).name))

                        Vol = spm_vol(['beta_' thisTR '.nii']);
                        beta = spm_read_vols(Vol);
                        
                        betas_this_sphere = beta(thisroi);
                        
                        betas_total(TR, cond, rn, ss) = mean(betas_this_sphere);

                    end
                    
                end
                
            end
            
        end

        %% Calculating mean and std error
        betas_total_across_run = mean(betas_total, 3); % average across runs
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
        name = strsplit(Region{rr}, ' ');
        name_scan = regexp(scan(sc).name, '_', 'split'); name_scan = strjoin(name_scan, '\\_');
        name_title = strjoin([name_scan, name], '\\_');
        name_file = strjoin([scan(sc).name, name], '_');
        
        xlim([0 scan(sc).num_TR+1]) 
        set(gca,'XTick',1:scan(sc).num_TR)

        label_x = cell(1, scan(sc).num_TR);
%         label_x{1} = 0;
        if strcmp(scan(sc).name, 'isss')
            for ii = 1:scan(sc).num_TR
                label_x{ii} = num2str(ii*2);
            end
        else
            for ii = 1:scan(sc).num_TR
                label_x{ii} = num2str(ii);
            end
        end
        
        set(gca,'XTickLabel', label_x)

        title(name_title, 'FontWeight', 'bold', 'Fontsize', 14)

        hleg1=legend([h3,h4,h1,h2],'ORA','SRA','NOI','SIL','Location','NorthEast');

        xlabel('Post Stimulus (TR)', 'Fontsize', 14);
        ylabel('Beta Estimates', 'Fontsize', 14);

        cd(dir_docs);

         %save in svg format 
         saveas(gcf,['beta_plot_sphere_6mm_' name_file], 'png');

         %save in eps format 
    %      saveas(gcf, ['beta_plot_' Region{ROI}], 'eps')

    %         clf;
    %     

        clear std_dev; 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    end
    
end

toc
