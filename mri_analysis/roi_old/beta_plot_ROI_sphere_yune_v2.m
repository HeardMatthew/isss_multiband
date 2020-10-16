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
condition = {'NOI', 'SIL', 'ORA', 'SRA'};

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

results_betas = cell(2, 3); % ROI by SCAN
results_error = cell(2, 3); % ROI by SCAN
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
            
            dir_subj = fullfile(dir_data, thissubj);
            dir_design = fullfile(dir_subj, 'design');
            cd(fullfile(dir_design, scan(sc).name))
            load SPM.mat
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
                
                for cond = 1:numCond % For each condition... (NOI, SIL, ORA, SRA)
                    
                    thistag = ['Sn(' num2str(rn) ') ' condition{cond}];
                    thesebetas = contains(SPM.xX.name, thistag);
                    vbetas = SPM.Vbeta(thesebetas);
                    
                    for TR = 1:length(vbetas)
                        beta = spm_read_vols(vbetas(TR));
                        betas_this_sphere = beta(thisroi & mask_master);
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
                foo = reshape(betas_total_across_run(TR, cond, 1, :), [1 11]);
                std_dev(TR,cond) = std(foo);
            end
            
        end

        std_error = std_dev/sqrt(length(subjects));
        
        results_betas{rr, sc} = betas_total_across_run_sub; 
        results_error{rr, sc} = std_error; 

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
         saveas(gcf,['beta_plot_yune_6mmsphere_' name_file], 'png');

         %save in eps format 
    %      saveas(gcf, ['beta_plot_' Region{ROI}], 'eps')

    %         clf;
    %     

        clear std_dev; 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    end
    
end

%% make things pretty
contrast = repelem({'noise'; 'silence'; 'or'; 'sr'}, 10);
betas = [results_betas{1, 1}(:, 1); results_betas{1, 1}(:, 2); results_betas{1, 1}(:, 3); results_betas{1, 1}(:, 4)];
T_lifg_hybrid = table(contrast, betas); 

contrast = repelem({'noise'; 'silence'; 'or'; 'sr'}, 5);
betas = [results_betas{1, 2}(:, 1); results_betas{1, 2}(:, 2); results_betas{1, 2}(:, 3); results_betas{1, 2}(:, 4)];
T_lifg_isss = table(contrast, betas); 

contrast = repelem({'noise'; 'silence'; 'or'; 'sr'}, 10);
betas = [results_betas{1, 3}(:, 1); results_betas{1, 3}(:, 2); results_betas{1, 3}(:, 3); results_betas{1, 3}(:, 4)];
T_lifg_multi = table(contrast, betas); 

contrast = repelem({'noise'; 'silence'; 'or'; 'sr'}, 10);
betas = [results_betas{2, 1}(:, 1); results_betas{2, 1}(:, 2); results_betas{2, 1}(:, 3); results_betas{2, 1}(:, 4)];
T_lpstg_hybrid = table(contrast, betas); 

contrast = repelem({'noise'; 'silence'; 'or'; 'sr'}, 5);
betas = [results_betas{2, 2}(:, 1); results_betas{2, 2}(:, 2); results_betas{2, 2}(:, 3); results_betas{2, 2}(:, 4)];
T_lpstg_isss = table(contrast, betas); 

contrast = repelem({'noise'; 'silence'; 'or'; 'sr'}, 10);
betas = [results_betas{2, 3}(:, 1); results_betas{2, 3}(:, 2); results_betas{2, 3}(:, 3); results_betas{2, 3}(:, 4)];
T_lpstg_multi = table(contrast, betas); 

cd 'C:\Users\heard.49\Documents\MATLAB\fmri'

writetable(T_lifg_hybrid, 'LIFG_hybrid.csv')
writetable(T_lifg_isss, 'LIFG_isss.csv')
writetable(T_lifg_multi, 'LIFG_multi.csv')
writetable(T_lpstg_hybrid, 'LPSTG_hybrid.csv')
writetable(T_lpstg_isss, 'LPSTG_isss.csv')
writetable(T_lpstg_multi, 'LPSTG_multi.csv')

%% 
contrast = repelem({'noise'; 'silence'; 'or'; 'sr'}, 10);
betas = [results_error{1, 1}(:, 1); results_error{1, 1}(:, 2); results_error{1, 1}(:, 3); results_error{1, 1}(:, 4)];
T_lifg_hybrid = table(contrast, betas); 

contrast = repelem({'noise'; 'silence'; 'or'; 'sr'}, 5);
betas = [results_error{1, 2}(:, 1); results_error{1, 2}(:, 2); results_error{1, 2}(:, 3); results_error{1, 2}(:, 4)];
T_lifg_isss = table(contrast, betas); 

contrast = repelem({'noise'; 'silence'; 'or'; 'sr'}, 10);
betas = [results_error{1, 3}(:, 1); results_error{1, 3}(:, 2); results_error{1, 3}(:, 3); results_error{1, 3}(:, 4)];
T_lifg_multi = table(contrast, betas); 

contrast = repelem({'noise'; 'silence'; 'or'; 'sr'}, 10);
betas = [results_error{2, 1}(:, 1); results_error{2, 1}(:, 2); results_error{2, 1}(:, 3); results_error{2, 1}(:, 4)];
T_lpstg_hybrid = table(contrast, betas); 

contrast = repelem({'noise'; 'silence'; 'or'; 'sr'}, 5);
betas = [results_error{2, 2}(:, 1); results_error{2, 2}(:, 2); results_error{2, 2}(:, 3); results_error{2, 2}(:, 4)];
T_lpstg_isss = table(contrast, betas); 

contrast = repelem({'noise'; 'silence'; 'or'; 'sr'}, 10);
betas = [results_error{2, 3}(:, 1); results_error{2, 3}(:, 2); results_error{2, 3}(:, 3); results_error{2, 3}(:, 4)];
T_lpstg_multi = table(contrast, betas); 

cd 'C:\Users\heard.49\Documents\MATLAB\fmri'

writetable(T_lifg_hybrid, 'error_LIFG_hybrid.csv')
writetable(T_lifg_isss, 'error_LIFG_isss.csv')
writetable(T_lifg_multi, 'error_LIFG_multi.csv')
writetable(T_lpstg_hybrid, 'error_LPSTG_hybrid.csv')
writetable(T_lpstg_isss, 'error_LPSTG_isss.csv')
writetable(T_lpstg_multi, 'error_LPSTG_multi.csv')

toc
