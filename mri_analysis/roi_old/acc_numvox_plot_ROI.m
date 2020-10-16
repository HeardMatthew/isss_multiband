%% acc_plot_ROI
% Make bar plots of MVPA accuracy

% MM/DD/YY  --  CHANGELOG
% 10/01/19  --  File inception. MJH
clearvars; clc; close all; 
%% FLAGS
do_extract = 1; 
do_calc = 1; 
do_plot = 1; 

%% paths and parameters
cd ..
isss_multi_params_v2_MVPA

scan(1).name = 'hybrid'; 
scan(1).full = 'hybrid_nowrong_unsmoothed_BETA_1_10'; 
scan(1).short = 'hybrid_nowrong_unsmoothed_BETA'; 
scan(1).second = 'GNB_hybrid_nowrong_unsmoothed_lng_noi_BETAS_1_10_2vox'; 

scan(2).name = 'isss'; 
scan(2).full = 'isss_nowrong_unsmoothed_BETA_1_5'; 
scan(2).short = 'isss_nowrong_unsmoothed_BETA'; 
scan(2).second = 'GNB_isss_nowrong_unsmoothed_lng_noi_BETAS_1_5_2vox'; 

scan(3).name = 'multi'; 
scan(3).full = 'multi_FIR_nowrong_unsmoothed_BETA_1_10'; 
scan(3).short = 'multi_FIR_nowrong_unsmoothed_BETA'; 
scan(3).second = 'GNB_multi_FIR_nowrong_unsmoothed_lng_noi_BETAS_1_10_2vox'; 

dir_thisdoc = 'C:\Users\heard.49\Documents\MATLAB\fmri\isss_multiband_preprocessed_04\docs\14subj_07232019\MVPA_numvox_plot'; 

filename = fullfile(dir_thisdoc, 'mvpa_numvox.mat'); 

%% 1. get number of voxels
if do_extract
    disp('Making mask...')
    dir_second = fullfile(dir_data, 'second_level', {scan.second}); 
    
    for sc = 1:length(scan)
        % load scan mask
        thismask = fullfile(dir_second{sc}, 'mask.nii'); 
        Vmask = spm_vol(thismask); 
        ymask = spm_read_vols(Vmask); 
        
        if sc == 1
            mask_all = logical(ymask); 
        else
            mask_all = mask_all & logical(ymask); 
        end
        
    end
    
    % combine masks into master mask
    disp('done! now getting voxels for subjects')
    
    % preallocate data storage
    numscan = length(scan);
    numsubj = length(subjects);
    
    data = nan(numscan, numsubj); % subject, scan, voxel
    
    for sc = 1:length(scan)
        for ss = 1:length(subjects) % 1-13, one subject only has one run!
            disp(subjects{ss})
            
            % load mvpa results for each subject
            thisfile = fullfile(dir_data, subjects{ss}, scan(sc).full, ...
                [subjects{ss}(1:2) '_' scan(sc).short '_search_GNB_lng_noi_rad2.nii']);
            Vacc = spm_vol(thisfile); 
            yacc = spm_read_vols(Vacc); 
            
            % get common voxels
            voxels = yacc(mask_all); 
            % count number above chance
            voxels_above_chance = length(find(voxels > 0.3)); 
            % poke it into the results matrix
            data(sc, ss) = voxels_above_chance; 
        end
        
    end
    
    % save as .mat to save time
    save(filename, 'data')
    
else
    load(filename)
end

%% 2. calculate average and within-subject error
if do_calc
    data_subjavg = mean(data, 2); 
    data_grandavg = mean(data_subjavg, 1); 
    data_new = data - data_subjavg + data_grandavg; 
    data_new_mean = mean(data_new, 2); 
    data_new_sdev = std(data_new, 0, 2); 
    data_new_serr = data_new_sdev / sqrt(length(subjects)); 
    
    save(filename, 'data', 'data_new_mean', 'data_new_serr')
else
    load(filename)
    if ~exist('data_new_mean', 'var')
        error('have not saved mean/standard error accuracy!')
    end
    
end

%% 3. make a pretty picture
if do_plot
    bar(data_new_mean)
    hold on
    er = errorbar(1:3, data_new_mean, data_new_serr); 
    er.Color = [0 0 0];
    er.LineStyle = 'none'; 
    line([0, 4], [0.5, 0.5], 'Color', 'red', 'LineStyle', '--')
%     set(gca, 'YLim', [22000 29000])
    xticklabels({'Fast Silent', 'Slow Silent', 'Fast Loud'})
    set(gca, 'XLim', [0.5 3.5])
    saveas(gca, 'acc_barplot', 'svg')
end

%% 4. new figure: # significant voxels





