%%% beta plot
%%% written by ysl 06/27/2014
%%% updated by mjh 11/22/2017
%%% updated to sphere by mjh 06/14/2018
% V4 -- Added TPM gray matter mask. MJH
% V5 -- Removed TPM, added LNG, added (mostly fixed) SPA. Combined LIFG. MJH
% SPAfigure -- generating SPA from the betas plotted on the figures. 
% V6 -- Down to 9 ROIs, only comparing hybrid and multiband. No
%   normalization. SPA metric generated from the mean betas calculated for
%   the ROI. 
% V7 -- Added isss back in
close all; clc;
tic

%% Flags
%%%%%%%%%%%%%%%%
do_extract = 0;%
             %%%
do_mask  = 0;%
do_betas = 1;% 
do_spa   = 1;%
%%%%%%%%%%%%%
do_plot = 1;%
            %%%%%
plot_by_roi = 1;%
plot_by_con = 0;%
plot_extra  = 0;%
        %%%%%%%%%
png = 0;%
svg = 1;%
eps = 0;% Kills color for some reason?
%%%%%%%%%%%%
do_ptt = 0;%
%%%%%%%%%%%%

these_cond = [1 2 3 4 5]; 
these_rois = [3 5]; 
% these_rois = 1:9;  
% was 14 total, now 13 with combined LIFG 
% Down to 9 for the manuscript -- 08/09/19 MJH

%% Paths and parameters
home = pwd; 
dir_rois = 'C:\Users\heard.49\Documents\MATLAB\fmri\isss_multiband_preprocessed_04\rois\v4_080119\mask_final'; 
dir_docs = 'C:\Users\heard.49\OneDrive\Documents\isss_multiband\manuscript\figures\workbench\plots';
% called in master (paired_T_norm_ROI.m) 
dir_docs_betas = fullfile(dir_docs, 'beta_plots_v6'); 

scan(1).name = 'hybrid_nowrong'; 
scan(1).num_TC_per_run = 180;
scan(1).num_TR = 10;
scan(1).num_runs = 2;
scan(1).full = 'Fast Silent:'; 

scan(2).name = 'isss_nowrong'; 
scan(2).num_TC_per_run = 90;
scan(2).num_TR = 5;
scan(2).num_runs = 2;
scan(2).full = 'Slow Silent:'; 

scan(3).name = 'multi_FIR_nowrong';
scan(3).num_TC_per_run = 180;
scan(3).num_TR = 10;
scan(3).num_runs = 2;
scan(3).full = 'Fast Loud:'; 

condition = { ...
    'NOI', ... 
    'SIL', ... 
    'LNG', ... 
    'ORA', ...
    'SRA', ...
    };


cond_full = { ...
    '1ch Vocoded Noise:', ...
    'Silence:', ...
    'Correct sentences:', ...
    'Correct OR sentences:', ...
    'Correct SR sentences:', ...
    }; 

numCond = length(condition); 

files_roi = dir(fullfile(dir_rois, '*.nii')); % Is now nine. 

%% Extract beta data
if do_extract    
    %% Prepare ROIs
    folders_roi = repelem({dir_rois}, length(files_roi))';
    rois = fullfile(folders_roi, {files_roi.name}'); 
    
    if do_mask
        mask_master = cell(1, length(rois)); 
        numvox      = nan(1, length(rois)); 
        
        for rr = these_rois
            V = spm_vol(rois{rr});
            [mask_roi, xyz] = spm_read_vols(V); mask_roi = logical(mask_roi); 

            for ss = 1:length(subjects)
                for sc = 1:length(scan)

                    thissubj = subjects{ss};
                    dir_subj = fullfile(dir_data, thissubj);
                    dir_thisdesign = fullfile(dir_subj, 'design');

                    cd(fullfile(dir_thisdesign, scan(sc).name));

                    V = spm_vol('mask.nii');
                    thismask = spm_read_vols(V);
                    if ss == 1 && sc == 1
                        mask_brain = thismask;
                    else
                        mask_brain = mask_brain & thismask;
                    end
                    
                end

            end

            mask_master{rr} = mask_roi & mask_brain; 
            numvox(rr) = length(find(mask_master{rr})); 
        end
        
        filename = fullfile(dir_rois, 'roi_masks_v7.mat'); 
        save(filename, 'mask_master', 'numvox');      
    else
        load(fullfile(dir_rois, 'roi_masks_v7.mat'))
    end
    
    %% Extract ROI betas
    if do_betas
        betas_mean_all = cell(length(scan), length(rois));
        betas_std_all  = cell(length(scan), length(rois));
        betas_se_all   = cell(length(scan), length(rois));
        % Dimensions [SCAN, ROI]

        betas_each_all = cell(length(rois), length(scan), 10, numCond, 2, length(subjects));
        % Dimensions [ROI, SCAN, TR, CON, RUN, SUBJ]
                
        for rr = these_rois 
            for sc = 1:length(scan)
                disp(scan(sc).name)
                numTRs = scan(sc).num_TR; 

                betas_total = nan(numTRs, numCond, 2, length(subjects));
                % Dimensions [TR, CON, RUN, SUBJ]

                for ss = 1:length(subjects) % For each subject...
                    disp(subjects{ss});
                    thisdir = fullfile(dir_data, subjects{ss}, 'design', scan(sc).name);
                    thisspm = fullfile(thisdir, 'SPM.mat');
                    load(thisspm)

                    %% Get data
                    if strcmp(subjects{ss}, 'JV_30Aug17')
                        maxrun = 1;
                    else
                        maxrun = 2;
                    end

                    for rn = 1:maxrun % For each run 
                        disp(['Run: ', num2str(rn)])

                        for cc = these_cond % For each condition... (NOI, SIL, ORA, SRA)
                            betas_this_roi = nan(numvox(rr), scan(sc).num_TR); 
                            % Dimensions [voxel, TR]
                            
                            if cc ~= 3 % NOI SIL ORA SRA
                                thistag = ['Sn(' num2str(rn) ') ' condition{cc}];
                                thesebetas = contains(SPM.xX.name, thistag);
                                vbetas = SPM.Vbeta(thesebetas);
                                
                                for TR = 1:length(vbetas)
                                    thisv = vbetas(TR); thisv.fname = fullfile(thisdir, thisv.fname);
                                    beta = spm_read_vols(thisv);
                                    betas_this_roi(:, TR) = beta(mask_master{rr});
                                    betas_each_all{rr, sc, TR, cc, rn, ss} = betas_this_roi(:, TR);
                                    betas_total(TR, cc, rn, ss) = mean(betas_this_roi(:, TR));
                                end
                            
                            else % LNG
                                tag1 = ['Sn(' num2str(rn) ') ORA'];
                                tag2 = ['Sn(' num2str(rn) ') SRA'];
                                betas1 = contains(SPM.xX.name, tag1);
                                betas2 = contains(SPM.xX.name, tag2);
                                vbetas1 = SPM.Vbeta(betas1);
                                vbetas2 = SPM.Vbeta(betas2);
                                
                                for TR = 1:length(vbetas1)
                                    thisv = vbetas1(TR); thisv.fname = fullfile(thisdir, thisv.fname);
                                    beta1 = spm_read_vols(thisv);
                                    
                                    thisv = vbetas2(TR); thisv.fname = fullfile(thisdir, thisv.fname);
                                    beta2 = spm_read_vols(thisv);
                                    
                                    beta = (beta1 + beta2) ./ 2; % because mean() doesn't work?

                                    betas_this_roi(:, TR) = beta(mask_master{rr});
                                    betas_each_all{rr, sc, TR, cc, rn, ss} = betas_this_roi(:, TR);
                                    betas_total(TR, cc, rn, ss) = mean(betas_this_roi(:, TR));
                                end
                                
                            end

                        end

                    end

                    if strcmp(subjects{ss}, 'JV_30Aug17') % quick and dirty fix as subj has 1 run
                        betas_total(:, :, 2, ss) = betas_total(:, :, 1, ss); 
                        betas_each_all(rr, sc, :, :, 2, ss) = betas_each_all(rr, sc, :, :, 1, ss); 
                    end

                end

                %% Get mean and std error (betas)
                betas_total_across_run = mean(betas_total, 3); % average across runs
                betas_total_across_run_sub = mean(betas_total_across_run, 4);
                betas_std_dev = nan(numTRs, numCond);
                for cc = these_cond
                    for TR = 1:numTRs
                        foo = reshape(betas_total_across_run(TR, cc, 1, :), [1 length(subjects)]);
                        betas_std_dev(TR,cc) = std(foo);
                    end

                end

                betas_std_error = betas_std_dev/sqrt(length(subjects));

                betas_mean_all{sc, rr} = betas_total_across_run_sub; 
                betas_std_all{sc, rr}  = betas_std_dev; 
                betas_se_all{sc, rr}   = betas_std_error; 

            end

        end
    
        cd(dir_rois)
        save('roi_betas_v7.mat', 'betas_mean_all', 'betas_std_all', 'betas_se_all', 'betas_each_all')
    else
        cd(dir_rois)
        load('roi_betas_v7.mat')
    end
    
    %% Calculate (and export) SPA
    % 1. get betas for each voxel in an roi
    % 2. calculate mean betas in the ROI
    % 3. calculate SPA of mean betas in the roi per run
    % 4. average SPA across runs
    % 5. average SPA across subjects
    if do_spa
        SPA_mean_all = cell(length(scan), 1); 
        SPA_std_all = cell(length(scan), 1); 
        SPA_se_all = cell(length(scan), 1); 
        SPA_each_all = cell(length(scan), 1); 
        
        for sc = 1:length(scan)
            numTRs = scan(sc).num_TR; 
            % Dimensions [ROI, SUB, CON]
            SPA_all = nan(numCond, length(rois), length(subjects)); 
            for ss = 1:length(subjects)
                if strcmp(subjects{ss}, 'JV_30Aug17')
                    maxrun = 1;
                else
                    maxrun = 2;
                end      
                
                for rr = these_rois
                    % Verify size of ROI is consistent!
                    numvox = betas_each_all(rr, :, :, :, :, :); 
                    numvox = cellfun(@length, numvox); numvox = reshape(numvox, [], 1); 
                    numvox = numvox(numvox ~= 0); 
                    if ~all(numvox(1) == numvox)
                        error('Something is wrong with size!!!')
                    else
                        numvox = numvox(1); 
                    end
                    
                    %% New SPA calculation
                    for cc = these_cond
                        thesebetas = betas_each_all(rr, sc, :, cc, :, ss); 
                        % Recall, dimensions [ROI, SCAN, TR, CON, RUN, SUBJ]
                        thesebetas = reshape(thesebetas, [10 2]); 
                        % Now dimensions [TR, RUN]
                        
                        % Dimensions [VOX, RUN]
                        SPA = nan(2, 1); 
                        % Dimensions [RUN]
                        
                        for rn = 1:maxrun
                            betas = thesebetas(:, rn); % still a cell
                            betas = [betas{:}]; % now a matrix, yay!
                            % Dimensions [VOX, TR]
                            
                            betas_roi = mean(betas, 1); 
                            % Dimensions [TR]
                            
                            area_all = 0; 
                            for TR = 1:scan(sc).num_TR-1 % minus one because last sample calls on TR end-1 and end
                                y1=betas_roi(TR); y2=betas_roi(TR+1);

                                % If both are above zero...
                                if y1 > 0 && y2 > 0 % trapezoid shape
                                    area = (y1+ y2)/2;

                                % If both are below zero...
                                elseif y1 < 0 && y2 < 0
                                    % area= -1*(y1+ y2)/2;
                                    area = 0;

                                % If 1 is above and 2 is below
                                elseif y1 > 0 && y2 < 0
                                    b = y1;
                                    a = y2 - y1;
                                    x = -b/a;
                                    area = (x*y1/2);%+ ((1-x)*y2/2);

                                elseif y1 < 0 && y2 > 0
                                    b=y1;
                                    a=y2-y1;
                                    x=-b/a;
                                    %  area=(x*y1/2)+ ((1-x)*y2/2);
                                    area=((1-x)*y2/2);

                                end

                                area_all = area_all+area;
                            end
                            
                            SPA(rn) = area_all; 
                            if maxrun == 1
                                SPA(2) = SPA(1); 
                            end
                            
                        end 
                        
                        SPA_all(cc, rr, ss) = mean(SPA); 
                    end % each condition
                    
                end  % each roi
                
            end % each subject
            
            % Average across subjects
            SPA_across_subjects = mean(SPA_all, 3); 
            
            SPA_std_dev = nan(numCond,length(rois)); 
            for cc = these_cond
                for rr = 1:length(rois)
                    SPA_std_dev(cc, rr) = std(SPA_all(cc, rr, :)); 
                end
                
            end
            
            SPA_std_error = SPA_std_dev/sqrt(length(subjects)); 
            
            SPA_each_all{sc} = SPA_all; 
            SPA_mean_all{sc} = SPA_across_subjects; 
            SPA_std_all{sc} = SPA_std_dev;
            SPA_se_all{sc} = SPA_std_error;
        end % each scan
        
        cd(dir_rois)
        save('roi_SPA_v7.mat', 'SPA_mean_all', 'SPA_std_all', 'SPA_se_all', 'SPA_each_all')
    end
    
else
    cd(dir_rois)
    load('roi_betas_v7.mat')
    load('roi_SPA_v7.mat')
end

%% Plot results
if do_plot    
    cd(dir_docs); %cd beta_plots_v7
    
    roi_title = cell(length(files_roi), 1);
    roi_file = cell(length(files_roi), 1);
    
    for rr = 1:length(files_roi)
        name = strsplit(files_roi(rr).name, '_'); name = [name(1:end-1) strsplit(name{end}, '.')]; 
        roi_name = name(4:end-1);
        if length(roi_name) == 3
            roi_name = roi_name(1:2); 
        end
        
        if contains(roi_name{2}, 'preSMA')
            roi_name{2} = [roi_name{2}(1:3) '-' roi_name{2}(4:end)]; 
        end
        
        roi_name{1} = [upper(roi_name{1}(1)) lower(roi_name{1}(2:end))]; 
        roi_title{rr} = strjoin(roi_name, ' ');
        roi_file{rr} = strjoin(roi_name, '_');
    end
        
    if plot_by_roi
        for rr = these_rois
            for sc = 1:length(scan)
                numTR = scan(sc).num_TR; 
                %% Load results                
                betas_mean_noise = betas_mean_all{sc, rr}(:, 1); 
                betas_se_noise   = betas_se_all{sc, rr}(:, 1);
                
                betas_mean_silence = betas_mean_all{sc, rr}(:, 2); 
                betas_se_silence   = betas_se_all{sc, rr}(:, 2);
                
                betas_mean_ora = betas_mean_all{sc, rr}(:, 4); 
                betas_se_ora   = betas_se_all{sc, rr}(:, 4);
                
                betas_mean_sra = betas_mean_all{sc, rr}(:, 5); 
                betas_se_sra   = betas_se_all{sc, rr}(:, 5);
                
                SPA_mean_noise = SPA_mean_all{sc}(1, rr);
                SPA_se_noise   = SPA_se_all{sc}(1, rr);
                
                SPA_mean_silence = SPA_mean_all{sc}(2, rr);
                SPA_se_silence   = SPA_se_all{sc}(2, rr);
                
                SPA_mean_ora = SPA_mean_all{sc}(4, rr);
                SPA_se_ora   = SPA_se_all{sc}(4, rr);
                
                SPA_mean_sra = SPA_mean_all{sc}(5, rr);
                SPA_se_sra   = SPA_se_all{sc}(5, rr);
                
                %% Plotting results
                figure;
                hold on;

                h1 = plot(betas_mean_noise,'b'); %NOI
                errorbar(betas_mean_noise, betas_se_noise, 'b');

                h2 = plot(betas_mean_silence,'b','LineStyle','--'); %SIL
                errorbar(betas_mean_silence, betas_se_silence, 'b','LineStyle','--');

                h3 = plot(betas_mean_ora,'r'); %ORA
                errorbar(betas_mean_ora, betas_se_ora, 'r');

                h4 = plot(betas_mean_sra,'r','LineStyle','--'); %SRA
                errorbar(betas_mean_sra, betas_se_sra, 'r','LineStyle','--');
                
                %% Formatting axes
                xlim([0 numTR+1]) 
                if sc == 2
                    ylim([-0.08 0.20])
                else
                    ylim([-0.04 0.10])
                end
                
                set(gca,'XTick', 1:numTR)
                set(gca,'XTickLabel', 1:numTR)
                xlabel('TR', 'Fontsize', 14);
                ylabel('Beta Estimates', 'Fontsize', 14);
                
                %% Add SPA                
                spatxt = {'SPA values:'; ...
                    ['1ch Noise = ' num2str(SPA_mean_noise, 2) ' \pm ' num2str(SPA_se_noise, 2)]; ...
                    ['Silence = ' num2str(SPA_mean_silence, 2) ' \pm ' num2str(SPA_se_silence, 2)]; ...
                    ['OR = ' num2str(SPA_mean_ora, 2) ' \pm ' num2str(SPA_se_ora, 2)]; ...
                    ['SR = ' num2str(SPA_mean_sra, 2) ' \pm ' num2str(SPA_se_sra, 2)]; ...
                }; 
                if sc == 2
                    txt = text(0.2, 0.17, spatxt, 'HorizontalAlignment', 'left'); 
                else
                    txt = text(0.2, 0.085, spatxt, 'HorizontalAlignment', 'left'); 
                end
                
                %% Legend
                hleg = legend([h1, h2, h3, h4], ...
                    '1ch Vocoded Noise','Silence','Correct OR','Correct SR','Location','NorthEast');                
                
                %% Titles
                name_title = strjoin([scan(sc).full, roi_title(rr)], ' ');
                name_file = strjoin([scan(sc).name, roi_file(rr)], '_');
                title(name_title, 'FontWeight', 'bold', 'Fontsize', 14)

                %% Save as PNG or SVG or EPS
                if png
                    cd png; saveas(gcf,['beta_plot_' name_file], 'png'); cd ..
                end
                
                if svg
%                     cd svg; 
                    cd(dir_docs)
                    saveas(gcf,['beta_plot_' name_file], 'svg'); cd ..
                end
                
                if eps
                    cd eps; saveas(gcf,['beta_plot_' name_file], 'eps'); cd ..
                end

            end
        
        end
        
    end
    
    if plot_by_con
        for rr = these_rois
            for cc = these_cond
                %% Load results                
                betas_mean_hybrid = betas_mean_all{1, rr}(:, cc); 
                betas_se_hybrid   = betas_se_all{1, rr}(:, cc);
                betas_std_hybrid  = betas_std_all{1, rr}(:, cc);
                
                betas_mean_isss = betas_mean_all{1, rr}(:, cc); 
                betas_se_isss   = betas_se_all{1, rr}(:, cc);
                betas_std_isss  = betas_std_all{1, rr}(:, cc);
                
                betas_mean_multi = betas_mean_all{2, rr}(:, cc); 
                betas_se_multi   = betas_se_all{2, rr}(:, cc);
                betas_std_multi  = betas_std_all{2, rr}(:, cc);
                
                SPA_mean_hybrid = SPA_mean_all{1}(cc, rr);
                SPA_se_hybrid   = SPA_se_all{1}(cc, rr);
                
                SPA_mean_multi = SPA_mean_all{2}(cc, rr);
                SPA_se_multi   = SPA_se_all{2}(cc, rr);
                
                %% SPA percent change
                SPA = [SPA_mean_hybrid, SPA_mean_multi]; 
                dSPA = -1*diff(SPA); 
                pSPA = -100*diff(SPA)/SPA(2); % calculate percent change relative to multiband (the default)
                
                %% Plotting results
                figure;
                hold on;
                
                h1 = plot(betas_mean_hybrid, 'b');
                errorbar(betas_mean_hybrid, betas_se_hybrid, 'b');

                h3 = plot(betas_mean_multi,'r'); 
                errorbar(betas_mean_multi, betas_se_multi, 'r');
                
                %% Formatting axes
                xlim([0 11]) 
                ylim([-0.04 0.1])
                set(gca,'XTick', 1:10)
                set(gca,'XTickLabel', 1:10)
                xlabel('TR', 'Fontsize', 14);
                ylabel('Beta Estimates', 'Fontsize', 14);
                
                %% Add SPA
                if dSPA < 0 % multi greater
                    color = '\color{red}'; 
                elseif dSPA > 0 % hybrid greater
                    dSPA = abs(dSPA); 
                elseif dSPA == 0 % no difference
                    pSPA = 0; 
                end
                
                spatxt = {'SPA values:'; ...
                    ['Interleaved = ' num2str(SPA(1), 3) '\pm ' num2str(SPA_se_hybrid, 3)]; ...
                    ['Continuous = ' num2str(SPA(2), 3) '\pm ' num2str(SPA_se_multi, 3)]; ...
                    }; 
                
                if dSPA < 0 % multi greater
                    spatxt = [spatxt(1:2); {['\color{red}' spatxt{3}]}]; 
                elseif dSPA > 0 % hybrid greater
                    spatxt = [{[spatxt{1} '\color{blue}']}; {[spatxt{2} '\color{black}']}; spatxt(3)]; 
                end
                
                txt = text(0.5, 0.09, spatxt, 'HorizontalAlignment', 'left'); 
                
                %% Legend
                hleg = legend([h1, h3],'Fast Interleaved','Fast Continuous','Location','NorthEast');                
                
                %% Titles
                name_title = strjoin([cond_full{cc}, roi_title(rr)], ' ');
                name_file = strjoin([condition{cc}, roi_file(rr)], '_');
                title(name_title, 'FontWeight', 'bold', 'Fontsize', 14)

                %% Save as PNG or SVG or EPS
                if png
                    cd png; saveas(gcf,['beta_plot_' name_file], 'png'); cd ..
                end
                
                if svg
                    cd svg; saveas(gcf,['beta_plot_' name_file], 'svg'); cd ..
                end
                
                if eps
                    cd eps; saveas(gcf,['beta_plot_' name_file], 'eps'); cd ..
                end

            end
            
        end
        
    end
    
    if plot_extra
        % add plot showing right pSTG, noise and silence, hybrid and multi
        % on the same figure
        
        thisroi = 8; 
        %% Load results                
        betas_mean_hybrid_noise = betas_mean_all{1, thisroi}(:, 1); 
        betas_se_hybrid_noise   = betas_se_all{1, thisroi}(:, 1);

        betas_mean_hybrid_silence = betas_mean_all{1, thisroi}(:, 2); 
        betas_se_hybrid_silence   = betas_se_all{1, thisroi}(:, 2);

        betas_mean_multi_noise = betas_mean_all{2, thisroi}(:, 1); 
        betas_se_multi_noise   = betas_se_all{2, thisroi}(:, 1);

        betas_mean_multi_silence = betas_mean_all{2, thisroi}(:, 2); 
        betas_se_multi_silence  = betas_se_all{2, thisroi}(:, 2);

        SPA_mean_hybrid_noise = SPA_mean_all{1}(1, thisroi);
        SPA_se_hybrid_noise   = SPA_se_all{1}(1, thisroi);

        SPA_mean_hybrid_silence = SPA_mean_all{1}(2, thisroi);
        SPA_se_hybrid_silence   = SPA_se_all{1}(2, thisroi);

        SPA_mean_multi_noise = SPA_mean_all{2}(1, thisroi);
        SPA_se_multi_noise   = SPA_se_all{2}(1, thisroi);

        SPA_mean_multi_silence = SPA_mean_all{2}(2, thisroi);
        SPA_se_multi_silence   = SPA_se_all{2}(2, thisroi);

        %% Plotting results
        figure;
        hold on;

        h1 = plot(betas_mean_hybrid_noise,'b'); %hybrid NOI
        errorbar(betas_mean_hybrid_noise, betas_se_hybrid_noise, 'b');

        h2 = plot(betas_mean_hybrid_silence,'b','LineStyle','--'); %hybrid SIL
        errorbar(betas_mean_hybrid_silence, betas_se_hybrid_silence, 'b','LineStyle','--');

        h3 = plot(betas_mean_multi_noise,'r'); %multi noi
        errorbar(betas_mean_multi_noise, betas_se_multi_noise, 'r');

        h4 = plot(betas_mean_multi_silence,'r','LineStyle','--'); %multi sil
        errorbar(betas_mean_multi_silence, betas_se_multi_silence, 'r','LineStyle','--');

        %% Formatting axes
        xlim([0 11]) 
        ylim([-0.03 0.06])
        set(gca,'XTick', 1:10)
        set(gca,'XTickLabel', 1:10)
        xlabel('TR', 'Fontsize', 14);
        ylabel('Beta Estimates', 'Fontsize', 14);

        %% Add SPA                
%         spatxt = {'SPA values:'; ...
%             ['Interleaved 1ch Noise = ' num2str(SPA_mean_hybrid_noise, 2) ' \pm ' num2str(SPA_se_hybrid_noise, 2)]; ...
%             ['Interleaved Silence = ' num2str(SPA_mean_hybrid_silence, 2) ' \pm ' num2str(SPA_se_hybrid_silence, 2)]; ...
%             ['Continuous 1ch Noise = ' num2str(SPA_mean_multi_noise, 2) ' \pm ' num2str(SPA_se_multi_noise, 2)]; ...
%             ['Continuous Silence = ' num2str(SPA_mean_multi_silence, 2) ' \pm ' num2str(SPA_se_multi_silence, 2)]; ...
%         }; 
%         txt = text(0.2, 0.085, spatxt, 'HorizontalAlignment', 'left'); 

        %% Legend
        hleg = legend([h1, h2, h3, h4], ...
            'Interleaved 1ch Noise','Interleaved Silence','Continuous 1ch Noise','Continuous Silence','Location','NorthEast');                

        %% Titles
        name_file = strjoin(['hybrid_multi', roi_file(thisroi)], '_');
        title('Fast Interleaved and Continuous: 1ch Noise and Silence', 'FontWeight', 'bold', 'Fontsize', 14)

        %% Save as PNG or SVG or EPS
        if png
            cd png; saveas(gcf,['beta_plot_' name_file], 'png'); cd ..
        end

        if svg
            cd svg; saveas(gcf,['beta_plot_' name_file], 'svg'); cd ..
        end

        if eps
            cd eps; saveas(gcf,['beta_plot_' name_file], 'eps'); cd ..
        end
        
        
    end
 
end

%% Paired t-tests
if do_ptt
    % Compare SPA across scans in each ROI
    cd(dir_roi)
    load roi_SPA_v6.mat
    data_hybrid = SPA_each_all{1};
    data_multi  = SPA_each_all{2};
    % Dimensions [CON, ROI, SUB]
    
    hypothesis = nan(length(roi), numCond); 
    pvalue = nan(length(roi), numCond); 
    
    stats = struct('tstat', [], 'df', [], 'sd', []); 
    num_com = length(roi) * 5; 
    stats(num_com).tstat = []; 
    
    idx = 1; 
    for rr = 1:length(roi)
        thishybrid = reshape(data_hybrid(:, rr, :),[5 14]); 
        thismulti  = reshape(data_multi(:, rr, :), [5 14]); 
        % Dimensions [CON, SUB]
        
        for cc = 1:5
            hybrid = thishybrid(cc, :); 
            multi  = thismulti(cc, :); 
            
            [hypothesis(rr, cc), pvalue(rr, cc), ~, stats(idx)] = ttest(hybrid, multi); 
            idx = idx + 1; 
        end
        
    end
    
    significant = pvalue < 0.05; 
    trending    = pvalue < 0.1; 
    
    [x, y] = find(significant); 
    for ii = 1:length(find(significant))
        disp(['Region ' roi{x(ii)} ' condition ' condition{y(ii)} ' is significant'])
        stats(sub2ind(size(significant), x(ii), y(ii)))
    end
    
    [x, y] = find(trending); 
    for ii = 1:length(find(trending))
        disp(['Region ' roi{x(ii)} ' condition ' condition{y(ii)} ' is trending'])
        stats(sub2ind(size(significant), x(ii), y(ii)))
    end
    
    noi = pvalue(:, 1);
    sil = pvalue(:, 2);
    lng = pvalue(:, 3);
    or = pvalue(:, 4);
    sr = pvalue(:, 5);
    
    t = table(roi, noi, sil, lng, or, sr)
    
end

cd(home)
toc
