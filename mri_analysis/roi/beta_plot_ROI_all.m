%% beta_plot_ROI_all
% Plot betas for clusters, on the same axis. Just does hybrid_multiband

function beta_plot_ROI_all(subj, study, dd, masks, domake, doplot)
close all

%% Check input
if ~isstruct(study) || ~isstruct(subj)
    error('subj and study are both struct!')
end

if ~isnumeric(dd)
    error('need to specify which design!')
end

if ~ischar(masks)  
    error('masks need to be strings!')
end

%% Pathing
dir_roi = pwd; 
cd ..

comparison = 'hybrid_multiband'; 
dir_masks  = fullfile(dir_roi, masks, comparison); 

%% Parameters
nsubjs = length(subj); 

design = study.design(dd); 
cond   = design.cond; 

scannames = {study.scan(:).runname}; 

if domake
    %% Scan-specific parameters and pathing
    nbetas = 8; 

    %% Load ROI mask
    disp('Loading hybrid mask...')

    fname = fullfile(dir_masks, 'hybrid_only_nary.nii');
    Vmask = spm_vol(fname); 
    map   = spm_read_vols(Vmask); % mask = logical(map); 

    temp = unique(map); roi_idx = temp(temp > 0); 
    nrois = length(roi_idx); 

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
        thisroi = roi_idx(nr); 

        mask_roi = map == thisroi;
        nvox = length(find(mask_roi)); 

        %% Preallocate
        betas_avg = zeros(nsubjs, nbetas); 

        for ns = 1:nsubjs
            % Subject-specific path and parameters
            thissubj = subj(ns); 
            disp(thissubj.name)
            nruns = thissubj.runs(1); % thankfully everyone has same number of runs across scans!
            dir_design = fullfile(study.path, 'data', thissubj.name, ...
                'design', [thesescans{ns} '_' design.name]); 

            % Grab some info from SPM.mat
            spmmat = fullfile(dir_design, 'SPM.mat'); 
            load(spmmat)
            beta_names = SPM.xX.name; 

            % Preallocate
            betas_roi_OR = nan(nvox, nbetas, nruns); 
            betas_roi_SR = nan(nvox, nbetas, nruns); 
            for nu = 1:nruns
                %% Identify which betas to load for this condition
                thisrun = ['Sn(' num2str(nu) ') '];
                target_OR  = [thisrun 'OR'];
                target_SR  = [thisrun 'SR'];
                these_betas_OR = contains(beta_names, target_OR); 
                these_betas_SR = contains(beta_names, target_SR); 

                Vbetas_OR = SPM.Vbeta(these_betas_OR); 
                Vbetas_SR = SPM.Vbeta(these_betas_SR); 
                for nb = 1:nbetas
                    %% Load data for this condition, for each subject
                    Vbetas_OR(nb).fname = fullfile(dir_design, Vbetas_OR(nb).fname);
                    Vbetas_SR(nb).fname = fullfile(dir_design, Vbetas_SR(nb).fname);
                    data_OR = spm_read_vols(Vbetas_OR(nb));  % loads whole brain
                    data_SR = spm_read_vols(Vbetas_SR(nb));  % loads whole brain
                    betas_roi_OR(:, nb, nu) = data_OR(mask_roi); % grabs relevant data
                    betas_roi_SR(:, nb, nu) = data_SR(mask_roi); % grabs relevant data
                end

            end

            %% Average across conditions
            betas_roi = (betas_roi_OR + betas_roi_SR) / 2; 
            
            %% Average across runs
            betas_avg_runs = mean(betas_roi, 3); 
            betas_avg_runs_roi = mean(betas_avg_runs, 1); 
            betas_avg(ns, :) = mean(betas_avg_runs_roi, 1); 

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

    data = table(COND, BETA, betas_avg, betas_std, betas_ser); 

    fname = fullfile(dir_masks, ['betas_' scanname '_' roistr '.mat']); 
    save(fname, 'data')

end


%% Plot the data
dir_docs = fullfile(study.path, 'docs', '042420_roi_plots'); 
if ~exist(dir_docs, 'file'); mkdir(dir_docs); end

for ns = 1 % quick and dirty solution
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
        data_hybrid = data; 

        fname = fullfile(dir_masks, ['betas_incongruent_' otherscan '_' roistr '.mat']); 
        load(fname)
        data_multi = data; 
        
        %% Get ROI name
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
        data_hybrid_NOI = data_hybrid(strcmp(data_hybrid.COND, 'NOI'), :);
        data_hybrid_OR  = data_hybrid(strcmp(data_hybrid.COND, 'OR'), :);  
        data_hybrid_SR  = data_hybrid(strcmp(data_hybrid.COND, 'SR'), :);
        
        data_multi_NOI = data_hybrid(strcmp(data_multi.COND, 'NOI'), :);
        data_multi_OR  = data_hybrid(strcmp(data_multi.COND, 'OR'), :);  
        data_multi_SR  = data_hybrid(strcmp(data_multi.COND, 'SR'), :);

        BETAS = data_hybrid_NOI.BETA; 

        figure;
        hold on;

        h1 = plot(BETAS, data_hybrid_NOI.betas_avg, 'r', 'Marker', 'o', 'LineStyle', '--'); %NOI
        errorbar(data_hybrid_NOI.betas_avg, data_hybrid_NOI.betas_ser, 'r', 'LineStyle', '--');

        h2 = plot(BETAS, data_hybrid_OR.betas_avg, 'r', 'Marker', '+', 'LineStyle', '--'); %ORA
        errorbar(data_hybrid_OR.betas_avg, data_hybrid_OR.betas_ser, 'r', 'LineStyle', '--');

        h3 = plot(BETAS, data_hybrid_SR.betas_avg, 'r', 'Marker', '+', 'LineStyle', '--'); %SRA
        errorbar(data_hybrid_SR.betas_avg, data_hybrid_SR.betas_ser, 'r', 'LineStyle','--');

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

close all
end



