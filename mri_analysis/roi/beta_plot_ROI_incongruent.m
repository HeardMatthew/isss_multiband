%% beta_plot_ROI
% Plots betas extracted from ROIs
% 
% MM/DD/YY -- CHANGELOG
% 06/27/14 -- Written by Yune
% 11/22/17 -- Updated by myself
% 04/23/20 -- Cloned for new isss_multiband, made universal
% 04/24/20 -- Incongruent version

function beta_plot_ROI_incongruent(subj, study, dd, masks, thesescans, domake, doplot)
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
        thatscan = thesescans(~strcmp(thesescans, scanname)); thatscan = thatscan{1}; 
        temp = contains(scannames, thatscan); 
        nbetas = study.scan(temp).epis - (2/study.scan(temp).TR); 

        %% Load ROI mask
        disp(['Loading ' thesescans{ns} '...'])

        fnames{ns} = fullfile(dir_masks, [scanname '_only_nary.nii']);
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

            %% Get voxels in this ROI
            disp(['ROI #' num2str(nr)])
            thisroi = thisidx(nr); 

            mask_roi = thismap == thisroi;

            for nc = 1:nconds
                % Preallocate
                betas_avg_in_roi = zeros(nsubjs, nbetas); 
                betas_std_in_roi = zeros(nsubjs, nbetas); 
                
                % Again, legibility
                thiscond = cond{nc}; 
                disp(thiscond)

                for nu = 1:nsubjs
                    % Subject-specific path and parameters (INCONGRUENT)
                    thissubj = subj(nu); 
                    disp(thissubj.name)
                    nruns = thissubj.runs(thisscan); 
                    dir_design = fullfile(study.path, 'data', thissubj.name, ...
                        'design', [thatscan '_' design.name]); 

                    % Grab some info from SPM.mat
                    spmmat = fullfile(dir_design, 'SPM.mat'); 
                    load(spmmat)
                    beta_names = SPM.xX.name; 
                    
                    % Load mask for this design
                    temp = fullfile(dir_design, 'mask.nii'); temp = spm_vol(temp); 
                    mask_scan  = spm_read_vols(temp); mask_scan = logical(mask_scan); 

                    mask_roi_scan = mask_roi & mask_scan; 
                    nvox = length(find(mask_roi_scan)); 
                    
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
                            betas_roi(:, nb, nn) = data(mask_roi_scan); % grabs relevant data
                        end

                    end

                    %% Average across runs
                    betas_avg_runs = mean(betas_roi, 3); 
                    betas_avg_runs_roi = mean(betas_avg_runs, 1); 
                    betas_avg_in_roi(nu, :) = mean(betas_avg_runs_roi, 1); 
                    
                end
                
                %% Append results onto previous results, ready for table
                betas_avg_subj = [betas_avg_subj, mean(betas_avg_in_roi, 1)];
                betas_std_subj = [betas_std_subj, std(betas_avg_in_roi, 1)];

            end
            
            %% Save results
            COND = repelem(cond, nbetas)'; 
            BETA = repmat(1:nbetas, [1 3])'; 
            
            betas_avg = betas_avg_subj'; 
            betas_std = betas_std_subj';
            betas_ser = betas_std./sqrt(nsubjs); 
            
            data = table(COND, BETA, betas_avg, betas_std, betas_ser);  %#ok<NASGU>

            fname = fullfile(dir_masks, ['betas_incongruent_' thatscan '_' roistr '.mat']); 
            save(fname, 'data')

        end

    end

end

%% Plot the data
if doplot
    dir_docs = fullfile(study.path, 'docs', '042420_roi_plots_incongruent'); 
    if ~exist(dir_docs, 'file'); mkdir(dir_docs); end
    for ns = 1:nscans
        %% Scan-specific parameters and pathing
        thisscan = contains(scannames, thesescans{ns}); 
        scanname = thesescans{ns}; 
        nbetas = study.scan(thisscan).epis - (2/study.scan(thisscan).TR); 
        thatscan = thesescans(~strcmp(thesescans, scanname)); thatscan = thatscan{1}; 
        
        temp = contains(scannames, thatscan); 
        nbetas = study.scan(temp).epis - (2/study.scan(temp).TR); 

        %% Get number of ROIs
        target = fullfile(dir_masks, '*.mat'); 
        files = dir(target); 
        filenames = {files(:).name}';
        nrois = length(find(contains(filenames, ['incongruent_' thatscan])));

        disp(['Plotting ' scanname '...'])
        for nr = 1:nrois
            %% Load data
            disp(['ROI #' num2str(nr)])
            if nr < 10
                roistr = ['roi0' num2str(nr)]; 
            else
                roistr = ['roi' num2str(nr)]; 
            end
            
            fname = fullfile(dir_masks, ['betas_incongruent_' thatscan '_' roistr '.mat']); 
            load(fname)
            
            fname = fullfile(dir_masks, ['unique_' scanname '_against_' thatscan '.txt']); 
            fid = fopen(fname); 
            idx = 1; 
            while 1
                temp = fgetl(fid); 
                if temp == -1; break; end
                txt{idx} = temp; idx = idx + 1; 
            end
            fclose(fid); 
            
            roi_names = txt';
            
            %% Create figure (remember, it's "that" data in "this" ROI)
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
            try
                name_title = ['INCON. ' thatscan ' vs. ' scanname ': ' roi_names{nr}];
            catch
                err
            end
            title(name_title, 'FontWeight', 'bold', 'Fontsize', 14)
            temp = strsplit(roi_names{nr}, ' '); temp = strjoin(temp, '_'); 
            temp = strsplit(temp, '/'); temp = strjoin(temp, '_'); 
            name_file = strjoin({'incongruent', scanname, 'against', thatscan, temp}, '_'); 
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
