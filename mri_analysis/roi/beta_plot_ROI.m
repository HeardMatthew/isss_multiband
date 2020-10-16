%% beta_plot_ROI
% Plots betas extracted from ROIs
% 
% MM/DD/YY -- CHANGELOG
% 06/27/14 -- Written by Yune
% 11/22/17 -- Updated by myself
% 04/23/20 -- Cloned for new isss_multiband, made universal

function beta_plot_ROI(subj, study, dd, masks, thesescans, domake, doplot)
close all;

%% check inputs
if ~isstruct(study) || ~isstruct(subj)
    error('subj and study are both struct!')
end

if ~isnumeric(dd)
    error('need to specify which design!')
end

if ~ischar(masks)  
    error('masks need to be strings!')
end

if ~iscell(thesescans)
    error('scans need to be cell!')
end

%% Pathing
dir_roi = pwd; 

cd ..
dir_batch  = pwd; 
comparison = strjoin(thesescans, '_'); 
dir_masks  = fullfile(dir_roi, masks, comparison); 

%% Parameters
nsubjs = length(subj); 
nscans = length(thesescans); 

design = study.design(dd); 
cond   = design.cond; 
nconds = length(cond); 

scannames = {study.scan(:).runname}; 

if domake
    %% Preallocate
    fnames = cell(nscans, 1); 
    Vmask  = cell(nscans, 1); 
    map    = cell(nscans, 1); 
    mask   = cell(nscans, 1); 
    roi_num = nan(nscans, 1); 
    roi_idx = cell(nscans, 1); 

    for ns = 1:nscans
        %% Scan-specific parameters and pathing
        thisscan = contains(scannames, thesescans{ns}); 
        scanname = thesescans{ns}; 
        nbetas = study.scan(thisscan).epis - (2/study.scan(thisscan).TR); 

        %% Load ROI mask
        disp(['Loading ' thesescans{ns} '...'])

        fnames{ns} = fullfile(dir_masks, [thesescans{ns} '_only_nary.nii']);
        Vmask{ns} = spm_vol(fnames{ns}); 
        map{ns} = spm_read_vols(Vmask{ns}); mask{ns} = logical(map{ns}); 

        temp = unique(map{ns}); roi_idx{ns} = temp(temp > 0); 
        roi_num(ns) = length(roi_idx{ns}); 

        % Make code more readable
        nrois = roi_num(ns); 
        thismap = map{ns}; 
        thisidx = roi_idx{ns}; 

        for nr = 1:nrois
            % For saving data...
            if nr < 10
                roistr = ['roi0' num2str(nr)]; 
            else
                roistr = ['roi' num2str(nr)]; 
            end
            
            betas_avg_subj = [];
            betas_std_subj = [];
            betas_ser_subj = []; 

            %% Get voxels in this ROI
            disp(['ROI #' num2str(nr)])
            thisroi = thisidx(nr); 

            mask_roi = thismap == thisroi;
            nvox = length(find(mask_roi)); 

            for nc = 1:nconds
                % Preallocate
                betas_avg = zeros(nsubjs, nbetas); 
                betas_std = zeros(nsubjs, nbetas); 

                % Again, legibility
                thiscond = cond{nc}; 
                disp(thiscond)

                for nu = 1:nsubjs
                    % Subject-specific path and parameters
                    thissubj = subj(nu); 
                    disp(thissubj.name)
                    nruns = thissubj.runs(thisscan); 
                    dir_design = fullfile(study.path, 'data', thissubj.name, ...
                        'design', [thesescans{ns} '_' design.name]); 

                    % Grab some info from SPM.mat
                    spmmat = fullfile(dir_design, 'SPM.mat'); 
                    load(spmmat)
                    beta_names = SPM.xX.name; 

                    % Preallocate
                    betas_roi = nan(nvox, nbetas, nruns); 

                    for nn = 1:nruns
                        %% Identify which betas to load for this condition
                        thisrun = ['Sn(' num2str(nn) ') '];
                        target  = [thisrun thiscond];
                        these_betas = contains(beta_names, target); 

                        Vbetas = SPM.Vbeta(these_betas); 
                        for nb = 1:nbetas
                            %% Load data for this condition, for each subject
                            Vbetas(nb).fname = fullfile(dir_design, Vbetas(nb).fname);
                            data = spm_read_vols(Vbetas(nb));  % loads whole brain
                            betas_roi(:, nb, nn) = data(mask_roi); % grabs relevant data
                        end

                    end

                    %% Average across runs
                    betas_avg_runs = mean(betas_roi, 3); 
                    betas_avg_runs_roi = mean(betas_avg_runs, 1); 

                    betas_avg(nu, :) = mean(betas_avg_runs_roi, 1); 
                    
                end
                
                %% Append results onto previous results, ready for table
                betas_avg_subj = [betas_avg_subj, mean(betas_avg, 1)];
                betas_std_subj = [betas_std_subj, std(betas_avg, 1)];

            end
            
            %% Save results
            COND = repelem(cond, nbetas)'; 
            BETA = repmat(1:nbetas, [1 3])'; 
            
            betas_avg = betas_avg_subj'; 
            betas_std = betas_std_subj';
            betas_ser = betas_std./sqrt(nsubjs); 
            
            data = table(COND, BETA, betas_avg, betas_std, betas_ser);  %#ok<NASGU>

            fname = fullfile(dir_masks, ['betas_' scanname '_' roistr '.mat']); 
            save(fname, 'data')

        end

    end

end

%% Plot the data
if doplot
    dir_docs = fullfile(study.path, 'docs', '042420_roi_plots'); 
    if ~exist(dir_docs, 'file'); mkdir(dir_docs); end
    for ns = 1:nscans
        %% Scan-specific parameters and pathing
        thisscan = contains(scannames, thesescans{ns}); 
        scanname = thesescans{ns}; 
        nbetas = study.scan(thisscan).epis - (2/study.scan(thisscan).TR); 
        otherscan = thesescans(~strcmp(thesescans, scanname)); otherscan = otherscan{1}; 

        %% Get number of ROIs
        target = fullfile(dir_masks, '*.mat'); 
        files = dir(target); 
        filenames = {files(:).name}';
        nrois = length(find(contains(filenames, ['betas_' scanname])));

        disp(['Plotting ' scanname '...'])
        for nr = 1:nrois
            %% Load data
            disp(['ROI #' num2str(nr)])
            if nr < 10
                roistr = ['roi0' num2str(nr)]; 
            else
                roistr = ['roi' num2str(nr)]; 
            end
            
            fname = fullfile(dir_masks, ['betas_' scanname '_' roistr '.mat']); 
            load(fname)
            
            fname = fullfile(dir_masks, ['unique_' scanname '_against_' otherscan '.txt']); 
            fid = fopen(fname); 
            idx = 1; 
            while 1
                temp = fgetl(fid); 
                if temp == -1; break; end
                txt{idx} = temp; idx = idx + 1; 
            end
            fclose(fid); 
            
            roi_names = txt';
            
            %% Create figure
            data_NOI = data(strcmp(data.COND, 'NOI'), :);
            data_OR  = data(strcmp(data.COND, 'OR'), :);  
            data_SR  = data(strcmp(data.COND, 'SR'), :);
            
            BETAS = data_NOI.BETA; 
            
            figure;
            hold on;

            h1 = plot(BETAS, data_NOI.betas_avg, 'b'); %NOI
            errorbar(data_NOI.betas_avg, data_NOI.betas_ser, 'b');

            h2 = plot(BETAS, data_OR.betas_avg, 'r'); %ORA
            errorbar(data_OR.betas_avg, data_OR.betas_ser, 'r');

            h3 = plot(BETAS, data_SR.betas_avg, 'r', 'LineStyle', '--'); %SRA
            errorbar(data_SR.betas_avg, data_SR.betas_ser, 'r', 'LineStyle','--');

        % hold off; 
        % axis off;

            %%for getting EPS without any labeling/ axis comments out below  
            name_title = [scanname ' vs. ' otherscan ': ' roi_names{nr}];
            title(name_title, 'FontWeight', 'bold', 'Fontsize', 14)
            temp = strsplit(roi_names{nr}, ' '); temp = strjoin(temp, '_'); 
            temp = strsplit(temp, '/'); temp = strjoin(temp, '_'); 
            name_file = strjoin({scanname, 'against', otherscan, temp}, '_'); 
            name_file = fullfile(dir_docs, name_file); 

            xlim([0 nbetas+1]) 
            set(gca,'XTick',1:nbetas)
            
%             label_x = cell(1, nbetas);
%             if strcmp(thisscan, 'isss')
%                 for ii = 1:scan(sc).num_TR
%                     label_x{ii} = num2str(ii*2);
%                 end
%             else
%                 for ii = 1:scan(sc).num_TR
%                     label_x{ii} = num2str(ii);
%                 end
%             end
% 
%             set(gca,'XTickLabel', label_x)

            % Legend
            hleg1 = legend([h2, h3 ,h1],'ORA','SRA','NOI','Location','NorthEast');

            xlabel('Beta #', 'Fontsize', 14);
            ylabel('Beta Estimates', 'Fontsize', 14); 
            
            %save as png for now
            saveas(gcf, name_file, 'png');
%             saveas(gcf,['beta_swau_plot_sphere_' num2str(sphere_radius) '_' name_file], 'svg');
        end
        
    end
    
end

close all
end
