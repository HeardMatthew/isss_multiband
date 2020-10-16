%%% beta plot
%%% written by ysl 06/27/2014
%%% updated by mjh 11/22/2017
%%% updated to sphere by mjh 06/14/2018
% V4 -- Added TPM gray matter mask. MJH

close all; clc;
tic

%% Flags
do_extract  = 1;
do_mask     = 0; 
do_betas    = 1; 
do_spa      = 0;

do_plot     = 0;
plot_by_roi = 0;
plot_by_con = 1;
png = 1;
svg = 1;

these_rois = 1:14; % 14 total. 
% these_rois = [10 11 14]; % These are the significant ones

%% Paths and parameters
dir_rois = 'C:\Users\heard.49\Documents\MATLAB\fmri\isss_multiband_preprocessed_04\rois\v4_080119\mask_final'; 
% dir_docs = 'C:\Users\heard.49\Documents\MATLAB\fmri\isss_multiband_preprocessed_04\docs\14subj_07232019\paired_T_clusters_tpm';
% called in master (paired_T_norm_ROI.m) 
dir_docs_betas = fullfile(dir_docs, 'beta_plots'); 

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

numCond = 5; %NOI, SIL, ORA, SRA
condition = {'NOI', 'SIL', 'LNG', 'ORA', 'SRA'};

files_roi = dir(fullfile(dir_rois, '*.nii')); 

%% Extract beta data
if do_extract    
    %% Prepare ROIs
    folders_roi = repelem({dir_rois}, length(files_roi))';
    rois = fullfile(folders_roi, {files_roi.name}'); 
    
    betas_all     = cell(3, length(rois));
    betas_std_all = cell(3, length(rois));
    betas_se_all  = cell(3, length(rois));
    SPA_all     = cell(3, length(rois));
    SPA_std_all = cell(3, length(rois));
    SPA_se_all  = cell(3, length(rois));
    
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
        
        filename = fullfile(dir_rois, 'roi_masks_small.mat'); 
        save(filename, 'mask_master', 'numvox');      
    else
        load(fullfile(dir_rois, 'roi_masks_small.mat'))
    end
    
    %% Extract ROI betas
    if do_betas
        betas_total_all = cell(length(rois), length(scan), 10, numCond, 2, length(subjects));
        % Dimensions [ROI, SCAN, TR, CON, RUN, SUBJ]
        
        beta_global_max = zeros(length(rois), length(scan));
        % Dimensions [ROI, SCAN]
                
        for rr = these_rois %these_rois
            for sc = 1:length(scan)
                disp(scan(sc).name)
                numTRs = scan(sc).num_TR; 

                betas_total = nan(numTRs, numCond, 2, length(subjects));
                % Dimensions [TR, CON, RUN, SUBJ]

                SPA_total = nan(numCond, 2, length(subjects));
                % Dimensions [CON, RUN, SUBJ]

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

                        for cc = 1:numCond % For each condition... (NOI, SIL, ORA, SRA)
                            thistag = ['Sn(' num2str(rn) ') ' condition{cc}];
                            thesebetas = contains(SPM.xX.name, thistag);
                            vbetas = SPM.Vbeta(thesebetas);

                            betas_this_roi = nan(numvox(rr), scan(sc).num_TR); 
                            % dimensions (voxel, TR)

                            for TR = 1:length(vbetas)
                                thisv = vbetas(TR); thisv.fname = fullfile(thisdir, thisv.fname);
                                beta = spm_read_vols(thisv);
                                if max(max(max(beta))) > beta_global_max(rr, sc)
                                    beta_global_max(rr, sc) = max(max(max(beta))); 
                                    disp(['Current max beta is ' num2str(beta_global_max(rr, sc))])
                                end

                                betas_this_roi(:, TR) = beta(mask_master{rr});
                                betas_total_all{rr, sc, TR, cc, rn, ss} = betas_this_roi(:, TR);
                                betas_total(TR, cc, rn, ss) = mean(betas_this_roi(:, TR));
                            end

                        end

                    end

                    if strcmp(subjects{ss}, 'JV_30Aug17') % quick and dirty fix as subj has 1 run
                        betas_total(:, :, 2, ss) = betas_total(:, :, 1, ss); 
                        betas_total_all(rr, sc, :, :, 2, ss) = betas_total_all(rr, sc, :, :, 1, ss); 
                    end

                end

                %% Normalize, get mean and std error (betas)
                betas_norm = betas_total ./ beta_global_max(rr, sc); 
                betas_total_across_run = mean(betas_norm, 3); % average across runs
                betas_total_across_run_sub = mean(betas_total_across_run, 4);
                betas_std_dev = nan(numTRs, numCond);
                for cc = 1:numCond
                    for TR = 1:numTRs
                        foo = reshape(betas_total_across_run(TR, cc, 1, :), [1 length(subjects)]);
                        betas_std_dev(TR,cc) = std(foo);
                    end

                end

                betas_std_error = betas_std_dev/sqrt(length(subjects));

                betas_all{sc, rr}     = betas_total_across_run_sub; 
                betas_std_all{sc, rr} = betas_std_dev; 
                betas_se_all{sc, rr}  = betas_std_error; 

            end

        end
    
        cd(dir_rois)
        save('roi_betas_v2_small.mat', 'betas_total_all', 'betas_all', 'betas_std_all', 'betas_se_all', 'beta_global_max')
    end
    
    %% Calculate SPA
    if do_spa
        SPA_total_all = cell(length(rois), length(scan), numCond, 2, length(subjects));
        % Dimensions [ROI, SCAN, CON, RUN, SUB]
        
        SPA_global_max = zeros(length(rois), length(scan)); 
        % Dimensions [ROI, SCAN]
        
        if ~do_betas
            cd(dir_rois)
            load('roi_betas_v2_small.mat')
        end
        
        for rr = these_rois %these_rois
            for sc = 1:length(scan)
                
                disp(scan(sc).name)
                numTRs = scan(sc).num_TR; 

                SPA_total = nan(numCond, 2, length(subjects));
                % Dimensions [CON, RUN, SUBJ]
                
                thesebetas = betas_total_all(rr, sc, :, :, :, :); 
                newsize = size(thesebetas); newsize = newsize(3:end); 
                thesebetas = reshape(thesebetas, newsize); 
                % Dimensions [TR, CON, RUN, SUB]
                
                thisnorm = beta_global_max(rr, sc); 
                thesebetas_norm = cellfun((@(x) x ./ thisnorm), thesebetas, 'UniformOutput', false); 
                
                for ss = 1:length(subjects)
                    for cc = 1:numCond
                        for rn = 1:2
                            SPA = nan(numvox(rr), 2); 
                            betas = thesebetas_norm(:, cc, rn, ss); betas = [betas{:}]; 

                            for vv = 1:numvox(rr)
                                area_all = 0;

                                for TR = 1:scan(sc).num_TR-1
                                    y1 = betas(vv, TR); y2 = betas(vv,TR+1);

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
                                        a = y2-y1;
                                        x = -b/a;
                                        area = (x*y1/2);
                                    elseif y1 < 0 && y2 > 0
                                        b = y1;
                                        a = y2-y1;
                                        x = -b/a;
                                        area=((1-x)*y2/2);
                                    end

                                    area_all = area_all+area;
                                end

                                SPA(vv, rn) = area_all;
                                

                                if max(max(max(SPA))) > SPA_global_max(rr, sc)
                                    SPA_global_max(rr, sc) = max(max(max(SPA))); 
                                    disp(['Current max SPA is ' num2str(SPA_global_max(rr, sc))])
                                end

                            end
                        
                            SPA_total_all{rr, sc, cc, rn, ss} = SPA;
                            SPA_total(cc, rn, ss) = mean(SPA(:, rn)); 
                            
                        end
                        
                    end

                end
                
                %% Calculating mean and std error (SPA)
%                 SPA_norm = SPA_total ./ SPA_global_max(rr, sc); 

                % Apply same scaling factor as we calculated earlier as we
                % are NOT working with scaled betas!
                
                
                SPA_total_across_run = mean(SPA_total, 2); % average across runs
                SPA_total_across_run_sub = mean(SPA_total_across_run, 3);
                SPA_std_dev = nan(1, numCond);
                for cc = 1:numCond
                    foo = reshape(SPA_total_across_run(cc, 1, :), [1 length(subjects)]);
                    SPA_std_dev(cc) = std(foo);
                end

                SPA_std_error = SPA_std_dev/sqrt(length(subjects));

                SPA_all{sc, rr}     = SPA_total_across_run_sub; 
                SPA_std_all{sc, rr} = SPA_std_dev; 
                SPA_se_all{sc, rr}  = SPA_std_error; 
            end
            
        end
        
        save('roi_SPA_small.mat', 'SPA_total_all', 'SPA_all', 'SPA_std_all', 'SPA_se_all', 'SPA_global_max')
    end
        
else
    cd(dir_rois)
    load('roi_betas_v2_small.mat')
    load('roi_SPA_small.mat')
end

%% Plot results
if do_plot    
    cd(dir_docs);
    cd beta_plots_formatted
    
    roi_title = cell(length(files_roi), 1);
    roi_file = cell(length(files_roi), 1);
    
    for rr = 1:length(files_roi)
        name = strsplit(files_roi(rr).name, '_'); name = [name(1:end-1) strsplit(name{end}, '.')]; 
        roi_name = name(2:end-1);
        roi_title{rr} = strjoin(roi_name, '\\_');
        roi_file{rr} = strjoin(roi_name, '_');
    end
        
    if plot_by_roi
        for rr = 1:length(files_roi)
            for sc = 1:length(scan)
                %% Plotting results
                figure;
                hold on;

                h1=plot(betas_all{sc, rr}(:,1),'b'); %NOI
                errorbar(betas_all{sc, rr}(:,1),se_all{sc, rr}(:,1), 'b');

                h2= plot(betas_all{sc, rr}(:,2),'b','LineStyle','--'); %SIL
                errorbar(betas_all{sc, rr}(:,2),se_all{sc, rr}(:,2), 'b','LineStyle','--');

                h3=plot(betas_all{sc, rr}(:,3),'r'); %ORA
                errorbar(betas_all{sc, rr}(:,3),se_all{sc, rr}(:,3), 'r');

                h4=plot(betas_all{sc, rr}(:,4),'r','LineStyle','--'); %SRA
                errorbar(betas_all{sc, rr}(:,4),se_all{sc, rr}(:,4), 'r','LineStyle','--');

                %%for getting EPS without any labeling/ axis comments out below  
                name_scan = regexp(scan(sc).name, '_', 'split'); name_scan = strjoin(name_scan, '\\_');
                name_title = strjoin([name_scan, roi_title(rr)], '\\_');
                name_file = strjoin([scan(sc).name, roi_file(rr)], '_');

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
                ylabel('Normalized Beta Estimates', 'Fontsize', 14);

                 %save in svg format 
                if png
                    saveas(gcf,['beta_plot_' name_file], 'png');
                elseif svg
                    saveas(gcf,['beta_plot_' name_file], 'svg');
                end

            end
        
        end
        
    end
    
    if plot_by_con
        for rr = these_rois% 1:length(roi_title)
            for cc = 3%1:length(condition)
                %% Load results
                beta_mean_hybrid = betas_all{1, rr}(:, cc); 
                beta_se_hybrid   = betas_se_all{1, rr}(:, cc);
                beta_std_hybrid  = betas_std_all{1, rr}(:, cc);
                
                SPA_mean_hybrid = SPA_all{1, rr}(cc); 
                SPA_se_hybrid   = SPA_se_all{1, rr}(cc); 
                SPA_std_hybrid  = SPA_std_all{1, rr}(cc); 
                
                mean_multi = betas_all{3, rr}(:, cc); 
                se_multi   = betas_se_all{3, rr}(:, cc);
                std_multi  = betas_std_all{3, rr}(:, cc);
                
                SPA_mean_multi = SPA_all{3, rr}(cc); 
                SPA_se_multi   = SPA_se_all{3, rr}(cc); 
                SPA_std_multi  = SPA_std_all{3, rr}(cc); 
                
                %% Plotting results
                figure;
                hold on;
                
                h1 = plot(beta_mean_hybrid, 'b');
                errorbar(beta_mean_hybrid, beta_se_hybrid, 'b');

                h3 = plot(mean_multi,'r'); 
                errorbar(mean_multi, se_multi, 'r');
                
                %% Add SPA
                yval = ylim; 
                yloc = min(yval) + diff(yval)/16; % 1/10th from bottom
                xloc = 10; 
                
                SPAtxt = { ...
                    ['Fast Interleaved SPA \pm SE = ' num2str(SPA_mean_hybrid, 3) ' \pm ' num2str(SPA_se_hybrid, 3)],  ...
                    ['Fast Continuous SPA \pm SE = ' num2str(SPA_mean_multi, 3) ' \pm ' num2str(SPA_se_multi, 3)] ...
                    }; 
                text(xloc, yloc, SPAtxt, 'HorizontalAlignment', 'right'); 
                
                %% Formatting axes
                xlim([0 11]) 
                set(gca,'XTick', 1:10)
                set(gca,'XTickLabel', 1:10)
                xlabel('Post Stimulus (TR)', 'Fontsize', 14);
                ylabel('Normalized Beta Estimates', 'Fontsize', 14);
                
                %% Legend
                hleg1 = legend([h1, h3],'Fast Interleaved','Fast Continuous','Location','NorthEast');                
                
                %% Titles
                name_title = strjoin([condition{cc}, roi_title(rr)], '\\_');
                name_file = strjoin([condition{cc}, roi_file(rr)], '_');
                title(name_title, 'FontWeight', 'bold', 'Fontsize', 14)

                %% Save as PNG or SVG
                if png
                    saveas(gcf,['beta_plot_' name_file], 'png');
                elseif svg
                    saveas(gcf,['beta_plot_' name_file], 'svg');
                end

            end
            
        end
        
    end
 
end

toc
