%%% beta plot
%%% written by ysl 06/27/2014
%%% updated by mjh 11/22/2017
%%% updated to sphere by mjh 06/14/2018

close all; clc;
tic

%% Flags
do_extract  = 0;
do_plot     = 1;
plot_by_roi = 0;
plot_by_con = 1;
png = 1;
svg = 0;

%% Paths and parameters
dir_rois = 'C:\Users\heard.49\Documents\MATLAB\fmri\isss_multiband_preprocessed_04\rois\v3_072319\mask_cluster';
dir_docs = 'C:\Users\heard.49\Documents\MATLAB\fmri\isss_multiband_preprocessed_04\docs\14subj_07232019\paired_T_clusters';

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

numCond = 4; %NOI, SIL, ORA, SRA
condition = {'NOI', 'SIL', 'ORA', 'SRA'};

files_roi = dir(fullfile(dir_rois, '*.nii')); 

%% Load ROIs
if do_extract
    folders_roi = repelem({dir_rois}, length(files_roi))';
    rois = fullfile(folders_roi, {files_roi.name}'); 
    betas_all = cell(3, length(rois));
    se_all = cell(3, length(rois));
    
    for rr = 1:length(rois)
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
                end

                mask_brain = mask_brain & thismask;
            end

        end

        mask_master = mask_roi & mask_brain; 

        for sc = 1:length(scan)
            disp(scan(sc).name)
            numTRs = scan(sc).num_TR; 
            betas_total = nan(numTRs, numCond, 2, length(subjects));
            beta_global_max = 0;
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

                    for cond = 1:numCond % For each condition... (NOI, SIL, ORA, SRA)
                        thistag = ['Sn(' num2str(rn) ') ' condition{cond}];
                        thesebetas = contains(SPM.xX.name, thistag);
                        vbetas = SPM.Vbeta(thesebetas);

                        for TR = 1:length(vbetas)
                            thisv = vbetas(TR); thisv.fname = fullfile(thisdir, thisv.fname);
                            beta = spm_read_vols(thisv);
                            if max(max(max(beta))) > beta_global_max
                                beta_global_max = max(max(max(beta))); 
                                disp(['Current max beta is ' num2str(beta_global_max)])
                            end

                            betas_this_roi = beta(mask_master);
                            betas_total(TR, cond, rn, ss) = mean(betas_this_roi);
                        end

                    end

                end

                if strcmp(subjects{ss}, 'JV_30Aug17') % quick and dirty fix as subj has 1 run
                    betas_total(:, :, 2, ss) = betas_total(:, :, 1, ss); 
                end

            end

            %% Normalize SPA
            betas_norm = betas_total ./ beta_global_max; 

            %% Calculating mean and std error
            betas_total_across_run = mean(betas_norm, 3); % average across runs
            betas_total_across_run_sub = mean(betas_total_across_run, 4);
            std_dev = nan(numTRs, numCond);
            for cond = 1:numCond
                for TR = 1:numTRs
                    foo = reshape(betas_total_across_run(TR, cond, 1, :), [1 length(subjects)]);
                    std_dev(TR,cond) = std(foo);
                end

            end

            std_error = std_dev/sqrt(length(subjects));

            betas_all{sc, rr} = betas_total_across_run_sub; 
            se_all{sc, rr} = std_error; 
        end
        
    end

    cd(dir_rois)
    save('roi_betas.mat', 'betas_all', 'se_all')
else
    cd(dir_rois)
    load('roi_betas.mat')
end

if do_plot    
    cd(dir_docs);
    if plot_by_roi
        roi_title = cell(length(files_roi), 1);
        roi_file = cell(length(files_roi), 1);
        for rr = 1:length(files_roi)
            name = strsplit(files_roi(rr).name, '_'); name = [name(1:end-1) strsplit(name{end}, '.')]; 
            roi_name = name(2:end-1);
            roi_title{rr} = strjoin(roi_name, '\\_');
            roi_file{rr} = strjoin(roi_name, '_');
        end
    
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
        for rr = 1:length(roi_title)
            for cc = 1:length(condition)
                %% Plotting results
                figure;
                hold on;

                h1 = plot(betas_all{1, rr}(:, cc), 'b'); %hybrid
                errorbar(betas_all{1, rr}(:, cc), se_all{1, rr}(:, cc), 'b');

    %             h2= plot(betas_all{2}(:,2),'b','LineStyle','--'); %isss
    %             errorbar(betas_all{sc}(:,2),se_all{sc}(:,2), 'b','LineStyle','--');

                h3 = plot(betas_all{3, rr}(:, cc),'r'); %multi
                errorbar(betas_all{3, rr}(:, cc), se_all{3, rr}(:, cc), 'r');

                %%for getting EPS without any labeling/ axis comments out below  
                name_title = strjoin([condition{cc}, roi_title(rr)], '\\_');
                name_file = strjoin([condition{cc}, roi_file(rr)], '_');
                
                xlim([0 scan(sc).num_TR+1]) 
                set(gca,'XTick',1:scan(sc).num_TR)

                label_x = cell(1, scan(sc).num_TR);
                for ii = 1:scan(sc).num_TR
                    label_x{ii} = num2str(ii);
                end

                set(gca,'XTickLabel', label_x)

                title(name_title, 'FontWeight', 'bold', 'Fontsize', 14)

                hleg1=legend([h1, h3],'Fast Interleaved','Continuous','Location','NorthEast');

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
 
end

toc
