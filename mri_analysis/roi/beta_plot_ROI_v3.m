%% beta_plot_ROI
% Plots betas extracted from ROIs
% 
% MM/DD/YY -- CHANGELOG
% 06/27/14 -- Written by Yune
% 11/22/17 -- Updated by myself
% 04/23/20 -- Cloned for new isss_multiband, made universal
% 05/01/20 -- Simplified for hybrid > multiband comparison
% 10/13/20 -- Here we go again... Version 3 made, runs 4 ROIs in 4 colors
% on 1 figure. 

function beta_plot_ROI_v3(varargin)
close all;

%% Check input
subj = varargin{1}; 
study = varargin{2}; 
dd = varargin{3}; 
masks = varargin{4}; 
if length(varargin) > 4
    thismask = varargin{5}; 
else
    thismask = '';
end

if ~isstruct(study) || ~isstruct(subj)
    error('subj and study are both struct!')
end

if ~isnumeric(dd)
    error('need to specify which design!')
end

if ~ischar(masks)  
    error('masks need to be strings!')
end

if ~ischar(thismask)  
    error('thismask need to be strings!')
end

%% Pathing
dir_roi = pwd; 

cd ..
comparison = 'hybrid_multiband'; 
% dir_masks  = fullfile(dir_roi, masks, comparison); 
dir_masks  = fullfile(dir_roi, masks); 

%% Path to results
dir_docs = fullfile(study.path, 'docs', '101320_hybrid_multi_plots'); 
if ~exist(dir_docs, 'file'); mkdir(dir_docs); end

%% Get number of ROIs
% More hard-coding
if contains(thismask, '/users/PAS1342/osu9912/fmri/atlases/BNA_PM_3D_246/')
    % If it's probabilistic...
    target = fullfile(dir_masks, 'betas_LNG_NOI*CAUDATE.mat');
    files  = dir(target); 
    filenames = {files(:).name}';
    numscans = length(filenames); 
    
    % Label of ROI
    roi_names = {'Left Caudate'}; 
    numrois   = 1; 
else
%     %% Get labels
%     fname = fullfile(dir_masks, 'unique_hybrid_against_multiband.txt'); 
%     fid = fopen(fname); 
%     idx = 1; 
%     while 1
%         temp = fgetl(fid); 
%         if temp == -1; break; end
%         txt{idx} = temp; idx = idx + 1; 
%     end
%     
%     fclose(fid); 
%     roi_names = txt';
%     numrois = length(txt); 
    %%% 101320 edits
    roi_names = cell(1, 10);
    roi_names{5} = 'anterior thalamus'; 
    roi_names{6} = 'right stg'; 
    roi_names{9} = 'right thalamus'; 
    roi_names{10} = 'brainstem'; 
    rois = [5, 6, 9, 10];
    
    rois_cort = [6, 9];
    rois_subc = [5, 10]; 
    % Hard-coded
    %%% 
end

disp('Plotting data...')

% just copy-paste hybrid and multiband separately, just get it done
colors = {'r', 'b'}; 

%% Hybrid cortical ROIs
figure;
hold on;

roi_idx = 1; 
h_NOI = cell(1, 2);
h_LNG = cell(1, 2);
for nroi = rois_cort %1:numrois
    %% Load data
    disp(['ROI #' num2str(nroi)])
    if nroi < 10
        roistr = ['roi0' num2str(nroi)]; 
    else
        roistr = ['roi' num2str(nroi)]; 
    end
    
    %% Hybrid data
    if contains(thismask, '/users/PAS1342/osu9912/fmri/atlases/BNA_PM_3D_246/')
        fname = fullfile(dir_masks, filenames{1}); 
    else
        fname = fullfile(dir_masks, ['betas_LNG_NOI_hybrid_' roistr '.mat']); 
    end
    
    load(fname)
    
    data_NOI = data(strcmp(data.COND, 'NOI'), :);
    data_LNG = data(strcmp(data.COND, 'LNG'), :);  

    BETAS = data_NOI.BETA; 

    h_NOI{roi_idx} = plot(BETAS, data_NOI.betas_avg, colors{roi_idx}, 'LineStyle', '--'); % NOI
    errorbar(data_NOI.betas_avg, data_NOI.betas_ser, colors{roi_idx}, 'LineStyle', '--');

    h_LNG{roi_idx} = plot(BETAS, data_LNG.betas_avg, colors{roi_idx}); % LNG
    errorbar(data_LNG.betas_avg, data_LNG.betas_ser, colors{roi_idx});
    
    roi_idx = roi_idx + 1; 
end

%%for getting EPS without any labeling/ axis comments out below  
name_title = 'hybrid cortical';
title(name_title, 'FontWeight', 'bold', 'Fontsize', 14)
% temp = strsplit(roi_names{nroi}, ' '); temp = strjoin(temp, '_'); 
% temp = strsplit(temp, '/'); temp = strjoin(temp, '_'); 
name_file = fullfile(dir_docs, 'hybrid_cortical'); 

xlim([0 9]) 
set(gca,'XTick',1:8)
ylim([-.03 .081])

% Legend
hleg1 = legend([h_LNG{1}, h_NOI{1}, h_LNG{2}, h_NOI{2}], 'right stg LNG', 'right stg NOI', 'right insula LNG', 'right insula NOI', 'Location', 'NorthEast');

xlabel('Beta #', 'Fontsize', 14);
ylabel('Beta Estimates', 'Fontsize', 14); 

%save as png for now
saveas(gcf, name_file, 'pdf');

%% Hybrid SUBcortical ROIs
figure;
hold on;

roi_idx = 1; 
h_NOI = cell(1, 2);
h_LNG = cell(1, 2);
for nroi = rois_subc %1:numrois
    %% Load data
    disp(['ROI #' num2str(nroi)])
    if nroi < 10
        roistr = ['roi0' num2str(nroi)]; 
    else
        roistr = ['roi' num2str(nroi)]; 
    end
    
    %% Hybrid data
    if contains(thismask, '/users/PAS1342/osu9912/fmri/atlases/BNA_PM_3D_246/')
        fname = fullfile(dir_masks, filenames{1}); 
    else
        fname = fullfile(dir_masks, ['betas_LNG_NOI_hybrid_' roistr '.mat']); 
    end
    
    load(fname)
    
    data_NOI = data(strcmp(data.COND, 'NOI'), :);
    data_LNG = data(strcmp(data.COND, 'LNG'), :);  

    BETAS = data_NOI.BETA; 

    h_NOI{roi_idx} = plot(BETAS, data_NOI.betas_avg, colors{roi_idx}, 'LineStyle', '--'); % NOI
    errorbar(data_NOI.betas_avg, data_NOI.betas_ser, colors{roi_idx}, 'LineStyle', '--');

    h_LNG{roi_idx} = plot(BETAS, data_LNG.betas_avg, colors{roi_idx}); % LNG
    errorbar(data_LNG.betas_avg, data_LNG.betas_ser, colors{roi_idx});
    
    roi_idx = roi_idx + 1; 
end

%%for getting EPS without any labeling/ axis comments out below  
name_title = 'hybrid subcortical';
title(name_title, 'FontWeight', 'bold', 'Fontsize', 14)
% temp = strsplit(roi_names{nroi}, ' '); temp = strjoin(temp, '_'); 
% temp = strsplit(temp, '/'); temp = strjoin(temp, '_'); 
name_file = fullfile(dir_docs, 'hybrid_subcortical'); 

xlim([0 9]) 
set(gca,'XTick',1:8)
ylim([-.03 .081])

% Legend
hleg1 = legend([h_LNG{1}, h_NOI{1}, h_LNG{2}, h_NOI{2}], 'thalamus LNG', 'thalamus NOI', 'brainstem LNG', 'brainstem NOI', 'Location', 'NorthEast');

xlabel('Beta #', 'Fontsize', 14);
ylabel('Beta Estimates', 'Fontsize', 14); 

%save as png for now
saveas(gcf, name_file, 'pdf');

%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%

%% Multiband cortical ROIs
figure;
hold on;

roi_idx = 1; 
h_NOI = cell(1, 2);
h_LNG = cell(1, 2);
for nroi = rois_cort %1:numrois
    %% Load data
    disp(['ROI #' num2str(nroi)])
    if nroi < 10
        roistr = ['roi0' num2str(nroi)]; 
    else
        roistr = ['roi' num2str(nroi)]; 
    end
    
    %% Hybrid data
    if contains(thismask, '/users/PAS1342/osu9912/fmri/atlases/BNA_PM_3D_246/')
        fname = fullfile(dir_masks, filenames{1}); 
    else
        fname = fullfile(dir_masks, ['betas_LNG_NOI_multiband_' roistr '.mat']); 
    end
    
    load(fname)
    
    data_NOI = data(strcmp(data.COND, 'NOI'), :);
    data_LNG = data(strcmp(data.COND, 'LNG'), :);  

    BETAS = data_NOI.BETA; 

    h_NOI{roi_idx} = plot(BETAS, data_NOI.betas_avg, colors{roi_idx}, 'LineStyle', '--'); % NOI
    errorbar(data_NOI.betas_avg, data_NOI.betas_ser, colors{roi_idx}, 'LineStyle', '--');

    h_LNG{roi_idx} = plot(BETAS, data_LNG.betas_avg, colors{roi_idx}); % LNG
    errorbar(data_LNG.betas_avg, data_LNG.betas_ser, colors{roi_idx});
    
    roi_idx = roi_idx + 1; 
end

%%for getting EPS without any labeling/ axis comments out below  
name_title = 'multiband cortical';
title(name_title, 'FontWeight', 'bold', 'Fontsize', 14)
% temp = strsplit(roi_names{nroi}, ' '); temp = strjoin(temp, '_'); 
% temp = strsplit(temp, '/'); temp = strjoin(temp, '_'); 
name_file = fullfile(dir_docs, 'multiband_cortical'); 

xlim([0 9]) 
set(gca,'XTick',1:8)
ylim([-.03 .081])

% Legend
hleg1 = legend([h_LNG{1}, h_NOI{1}, h_LNG{2}, h_NOI{2}], 'right stg LNG', 'right stg NOI', 'right insula LNG', 'right insula NOI', 'Location', 'NorthEast');

xlabel('Beta #', 'Fontsize', 14);
ylabel('Beta Estimates', 'Fontsize', 14); 

%save as png for now
saveas(gcf, name_file, 'pdf');

%% Multiband SUBcortical ROIs
figure;
hold on;

roi_idx = 1; 
h_NOI = cell(1, 2);
h_LNG = cell(1, 2);
for nroi = rois_subc %1:numrois
    %% Load data
    disp(['ROI #' num2str(nroi)])
    if nroi < 10
        roistr = ['roi0' num2str(nroi)]; 
    else
        roistr = ['roi' num2str(nroi)]; 
    end
    
    %% Hybrid data
    if contains(thismask, '/users/PAS1342/osu9912/fmri/atlases/BNA_PM_3D_246/')
        fname = fullfile(dir_masks, filenames{1}); 
    else
        fname = fullfile(dir_masks, ['betas_LNG_NOI_multiband_' roistr '.mat']); 
    end
    
    load(fname)
    
    data_NOI = data(strcmp(data.COND, 'NOI'), :);
    data_LNG = data(strcmp(data.COND, 'LNG'), :);  

    BETAS = data_NOI.BETA; 

    h_NOI{roi_idx} = plot(BETAS, data_NOI.betas_avg, colors{roi_idx}, 'LineStyle', '--'); % NOI
    errorbar(data_NOI.betas_avg, data_NOI.betas_ser, colors{roi_idx}, 'LineStyle', '--');

    h_LNG{roi_idx} = plot(BETAS, data_LNG.betas_avg, colors{roi_idx}); % LNG
    errorbar(data_LNG.betas_avg, data_LNG.betas_ser, colors{roi_idx});
    
    roi_idx = roi_idx + 1; 
end

%%for getting EPS without any labeling/ axis comments out below  
name_title = 'multiband subcortical';
title(name_title, 'FontWeight', 'bold', 'Fontsize', 14)
% temp = strsplit(roi_names{nroi}, ' '); temp = strjoin(temp, '_'); 
% temp = strsplit(temp, '/'); temp = strjoin(temp, '_'); 
name_file = fullfile(dir_docs, 'multiband_subcortical'); 

xlim([0 9]) 
set(gca,'XTick',1:8)
ylim([-.03 .081])

% Legend
hleg1 = legend([h_LNG{1}, h_NOI{1}, h_LNG{2}, h_NOI{2}], 'thalamus LNG', 'thalamus NOI', 'brainstem LNG', 'brainstem NOI', 'Location', 'NorthEast');

xlabel('Beta #', 'Fontsize', 14);
ylabel('Beta Estimates', 'Fontsize', 14); 

%save as png for now
saveas(gcf, name_file, 'pdf');

end
