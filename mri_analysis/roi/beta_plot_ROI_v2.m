%% beta_plot_ROI
% Plots betas extracted from ROIs
% 
% MM/DD/YY -- CHANGELOG
% 06/27/14 -- Written by Yune
% 11/22/17 -- Updated by myself
% 04/23/20 -- Cloned for new isss_multiband, made universal
% 05/01/20 -- Simplified for hybrid > multiband comparison
% 10/13/20 -- Here we go again...

function beta_plot_ROI_v2(varargin)
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
    % Hard-coded
    %%% 
end

disp('Plotting data...')
for nroi = rois %1:numrois
    %% Load data
    disp(['ROI #' num2str(nroi)])
    if nroi < 10
        roistr = ['roi0' num2str(nroi)]; 
    else
        roistr = ['roi' num2str(nroi)]; 
    end

    %% Open figure
    figure;
    hold on;
    
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

    h1 = plot(BETAS, data_NOI.betas_avg, 'r', 'LineStyle', '--'); % NOI
    errorbar(data_NOI.betas_avg, data_NOI.betas_ser, 'r', 'LineStyle', '--');

    h2 = plot(BETAS, data_LNG.betas_avg, 'r'); % LNG
    errorbar(data_LNG.betas_avg, data_LNG.betas_ser, 'r');

    %% Multiband data
    if contains(thismask, '/users/PAS1342/osu9912/fmri/atlases/BNA_PM_3D_246/')
        fname = fullfile(dir_masks, filenames{2}); % wait, won't this work below??
    else
        fname = fullfile(dir_masks, ['betas_LNG_NOI_multiband_' roistr '.mat']); 
    end
    
    load(fname)
    
    data_NOI = data(strcmp(data.COND, 'NOI'), :);
    data_LNG = data(strcmp(data.COND, 'LNG'), :);  

    h3 = plot(BETAS, data_NOI.betas_avg, 'k', 'LineStyle', '--'); % NOI
    errorbar(data_NOI.betas_avg, data_NOI.betas_ser, 'k', 'LineStyle', '--');

    h4 = plot(BETAS, data_LNG.betas_avg, 'k'); % LNG
    errorbar(data_LNG.betas_avg, data_LNG.betas_ser, 'k');

    %%for getting EPS without any labeling/ axis comments out below  
    name_title = roi_names{nroi};
    title(name_title, 'FontWeight', 'bold', 'Fontsize', 14)
    temp = strsplit(roi_names{nroi}, ' '); temp = strjoin(temp, '_'); 
    temp = strsplit(temp, '/'); temp = strjoin(temp, '_'); 
    name_file = fullfile(dir_docs, temp); 

    xlim([0 9]) 
    set(gca,'XTick',1:8)
    
    % Legend
    hleg1 = legend([h2, h1, h4, h3], 'Hybrid LNG', 'Hybrid NOI', 'Multi LNG', 'Multi NOI', 'Location', 'NorthEast');

    xlabel('Beta #', 'Fontsize', 14);
    ylabel('Beta Estimates', 'Fontsize', 14); 

    %save as png for now
    saveas(gcf, name_file, 'png');
%             saveas(gcf,['beta_swau_plot_sphere_' num2str(sphere_radius) '_' name_file], 'svg');
end

close all
end
