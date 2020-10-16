%% acc_plot_ROI
% Make bar plots of MVPA accuracy

% MM/DD/YY  --  CHANGELOG
% 10/01/19  --  File inception. MJH
clearvars; clc; close all; 
%% FLAGS
do_extract = 0; 
do_calc = 1; 
do_plot = 0; 

%% paths and parameters
cd ..
isss_multi_params_v2_MVPA

scan(1).name = 'hybrid'; 
scan(1).full = 'hybrid_nowrong_unsmoothed_BETA_1_10'; 
scan(1).short = 'hybrid_nowrong_unsmoothed_BETA'; 

scan(2).name = 'isss'; 
scan(2).full = 'isss_nowrong_unsmoothed_BETA_1_5'; 
scan(2).short = 'isss_nowrong_unsmoothed_BETA'; 

scan(3).name = 'multi'; 
scan(3).full = 'multi_FIR_nowrong_unsmoothed_BETA_1_10'; 
scan(3).short = 'multi_FIR_nowrong_unsmoothed_BETA'; 

dir_roi = 'C:\Users\heard.49\Documents\MATLAB\fmri\isss_multiband_preprocessed_04\rois\lng_noi_overlap_0930219\masks_all'; 

%% 1. get accuracy within overlapping voxels
filename = fullfile(dir_roi, 'mvpa_acc_overlap.mat'); 

if do_extract
    % Load mask
    file_mask = fullfile(dir_roi, 'mask_mvpa.nii'); 
    Vmask = spm_vol(file_mask); 
    [ymask, xyz] = spm_read_vols(Vmask); 
    mask = logical(ymask); 
    
    numvox = length(find(ymask ~= 0)); 
    numscan = length(scan);
    numsubj = length(subjects);
    
    data = nan(numvox, numscan, numsubj); % subject, scan, voxel
    
    for ss = 1:length(subjects) % 1-13, one subject only has one run!
        disp(subjects{ss})
        
        for sc = 1:length(scan)
            disp(scan(sc).name)
            
            % load mvpa results for each subject
            thisfile = fullfile(dir_data, subjects{ss}, scan(sc).full, ...
                [subjects{ss}(1:2) '_' scan(sc).short '_search_GNB_lng_noi_rad2.nii']);
            Vacc = spm_vol(thisfile); 
            yacc = spm_read_vols(Vacc); 
            
            % poke it into the results matrix
            data(:, sc, ss) = yacc(mask); 
        end
        
    end
    
    % save as .mat to save time
    save(filename, 'data')
else
    load(filename)
end

%% 2. calculate average and within-subject error
if do_calc
    acc_acrossvox = mean(data, 1); 
    acc_acrossvox_subjavg = mean(acc_acrossvox, 2); 
    acc_acrossvox_grandavg = mean(acc_acrossvox_subjavg, 3); 
    acc_acrossvox_new = acc_acrossvox - acc_acrossvox_subjavg + acc_acrossvox_grandavg; 
    acc_acrossvox_new_mean = mean(acc_acrossvox, 3); 
    acc_acrossvox_new_sdev = std(acc_acrossvox, 0, 3); 
    acc_acrossvox_new_serr = acc_acrossvox_new_sdev / sqrt(length(subjects)); 
    
    save(filename, 'data', 'acc_acrossvox_new_mean', 'acc_acrossvox_new_serr')
else
    load(filename)
    if ~exist('acc_acrossvox_new_mean', 'var')
        error('have not saved mean/standard error accuracy!')
    end
    
end

%% 3. make a pretty picture
if do_plot
    acc = acc_acrossvox_new_mean + 0.5; 
    bar(acc)
    hold on
    er = errorbar(1:3, acc, acc_acrossvox_new_serr); 
    er.Color = [0 0 0];
    er.LineStyle = 'none'; 
    line([0, 4], [0.5, 0.5], 'Color', 'red', 'LineStyle', '--')
    set(gca, 'YLim', [0.45 0.65])
    xticklabels({'Fast Silent', 'Slow Silent', 'Fast Loud'})
    set(gca, 'XLim', [0.5 3.5])
    saveas(gca, 'acc_barplot', 'svg')
end

%% 4. new figure: # significant voxels
data = [964, 187, 966]; 




