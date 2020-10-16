%% beta_plot_ROI
% Plots betas extracted from ROIs
% 
% MM/DD/YY -- CHANGELOG
% 06/27/14 -- Written by Yune
% 11/22/17 -- Updated by myself
% 04/23/20 -- Cloned for new isss_multiband, made universal
% 05/01/20 -- Simplified for hybrid > multiband comparison

function SPA_plot_ROI(varargin)
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
% This is hard-coded
comparison = 'hybrid_multiband'; 
dir_masks  = fullfile(dir_roi, masks, comparison); 

%% Plot the data
% This is hard-coded
dir_docs = fullfile(study.path, 'docs', '050420_hybrid_multi_plots'); 
if ~exist(dir_docs, 'file'); mkdir(dir_docs); end

%% Get number of ROIs
% More hard-coding
if contains(thismask, '/users/PAS1342/osu9912/fmri/atlases/BNA_PM_3D_246/')
    % If it's probabilistic...
    target = fullfile(dir_masks, 'SPA_LNG_NOI_CAUDATE.mat');
    files  = dir(target); 
    filenames = {files(:).name}';
    numrois = length(filenames); 
    
    % Label of ROI
    roi_names = {'Left Caudate'}; 
else
    target = fullfile(dir_masks, '*.mat'); 
    files  = dir(target); 
    filenames = {files(:).name}';
    numrois = length(find(contains(filenames, ['betas_LNG_NOI_hybrid'])));
    
    %% Get labels
    fname = fullfile(dir_masks, 'unique_hybrid_against_multiband.txt'); 
    fid = fopen(fname); 
    idx = 1; 
    while 1
        temp = fgetl(fid); 
        if temp == -1; break; end
        txt{idx} = temp; idx = idx + 1; 
    end
    
    fclose(fid); 
    roi_names = txt';
end

disp('Plotting data...')
for nroi = 1:numrois
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
    
    %% Load data
    if contains(thismask, '/users/PAS1342/osu9912/fmri/atlases/BNA_PM_3D_246/')
        fname = fullfile(dir_masks, filenames{nroi}); % wait, won't this work below??
    else
        fname = fullfile(dir_masks, ['SPA_LNG_NOI_' roistr '.mat']); 
    end
    
    load(fname)
    
    data_LNG = data(strcmp(data.COND, 'LNG'), :);  
    data_NOI = data(strcmp(data.COND, 'NOI'), :);
    
    avgdata = [data_LNG.SPA_avg'; data_NOI.SPA_avg']; 
    serdata = [data_LNG.SPA_ser'; data_NOI.SPA_ser']; 
    
%     X = categorical({'LNG','NOI'});
%     X = reordercats(X,{'LNG','NOI'});

    b = bar(avgdata, 'FaceColor', 'flat'); 
    % spits out structure with 2 values. 
    % Value 1 is hybrid, value 2 is multiband
    b(1).FaceColor = 'r';
    b(2).FaceColor = 'k';
    
    % Make error bars
    ngroups = size(serdata, 1);
    nbars = size(serdata, 2);
    % Calculating the width for each bar group
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    for i = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x, avgdata(:,i), serdata(:,i), '.', 'Color', [0.5 0.5 0.5]);
            
    end
    
    % Set X labels
    set(gca, 'XTick', [1 2])
    set(gca, 'XTickLabel', {'LNG', 'NOI'})
    
    % Legend
    l = legend(b, 'Hybrid', 'Multiband'); 

    % Title and fname
    name_title = roi_names{nroi};
    title(name_title, 'FontWeight', 'bold', 'Fontsize', 14)
    temp = strsplit(roi_names{nroi}, ' '); temp = strjoin(temp, '_'); 
    temp = strsplit(temp, '/'); temp = strjoin(['BAR', temp], '_'); 
    name_file = fullfile(dir_docs, temp); 
    
    xlabel('Condition', 'Fontsize', 14);
    ylabel('SPA value', 'Fontsize', 14); 
    
    hold off
    
    %save as png for now
    saveas(gcf, name_file, 'png');
%             saveas(gcf,['beta_swau_plot_sphere_' num2str(sphere_radius) '_' name_file], 'svg');
end

close all
end
