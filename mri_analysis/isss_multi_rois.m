%% isss_multi_rois
% Runs all ROI analyses
%
% MM/DD/YY -- CHANGELOG
% 04/23/20 -- Started the file. Moved ROI analysis from _process to here. 
% 10/13/20 -- Sisyphus was a graduate student. Three years in and I'm back
% here again: 
% 1. Make a list of hybrid ROIs (n-ary)
% 2. Identify which ones are the right STG, right anterior insula, anterior 
%    thalamus, and inferior colliculus
% 3. Grab betas from all subjects (language and 1ch noise)
% 4. Average betas across subjects, get the standard error
% 5. Make it look pretty. 
clc; clearvars 

%% Pathing
dir_batch = pwd;
dir_roi = fullfile(dir_batch, 'roi');
isss_multi_params

%% Make overlapping and unique, binary and n-ary masks
% disp('Making unique/overlapping ROIs')
% cd(dir_roi); find_overlapping(study, 'mask_lng_noi_042020', {'hybrid_lng_noi_rois', 'isss_lng_noi_rois', 'multiband_lng_noi_rois'})
% disp('Done!')

%% SNR ROI analysis
%     disp('SNR ROI analysis')
%     cd(dir_roi); roi_snr(subj, study, 'mask_lng_noi_042020', {'hybrid', 'isss'})
%     cd(dir_roi); roi_snr(subj, study, 'mask_lng_noi_042020', {'hybrid', 'multiband'})
%     disp('Done!')

%% Plot ROIs
% disp('Plotting ROIs...')
% cd(dir_roi); beta_plot_ROI(subj, study, 4, 'mask_lng_noi_042020', {'hybrid', 'isss'}, 0, 1)
% cd(dir_roi); beta_plot_ROI(subj, study, 4, 'mask_lng_noi_042020', {'hybrid', 'multiband'}, 0, 1)
% 
% cd(dir_roi); beta_plot_ROI_incongruent(subj, study, 4, 'mask_lng_noi_042020', {'hybrid', 'isss'}, 0, 1)
% cd(dir_roi); beta_plot_ROI_incongruent(subj, study, 4, 'mask_lng_noi_042020', {'hybrid', 'multiband'}, 0, 1)

% cd(dir_roi); beta_get_ROI(subj, study, 4, 'mask_lng_noi_042020')
% cd(dir_roi); beta_plot_ROI_v2(subj, study, 4, 'mask_lng_noi_042020')

% maskname = '/users/PAS1342/osu9912/fmri/atlases/BNA_PM_3D_246/r227.nii'; 
% cd(dir_roi); beta_get_ROI(subj, study, 4, 'mask_lng_noi_042020', maskname)
% cd(dir_roi); beta_plot_ROI_v2(subj, study, 4, 'mask_lng_noi_042020', maskname)

% cd(dir_roi); SPA_get_ROI(subj, study, 4, 'mask_lng_noi_042020', maskname)
% cd(dir_roi); SPA_plot_ROI(subj, study, 4, 'mask_lng_noi_042020', maskname)

maskname = 'hybrid_lng_noi_rois.nii'; 
% cd(dir_roi); beta_get_ROI(subj, study, 4, 'lng_noi_again_101320', maskname)
cd(dir_roi); beta_plot_ROI_v3(subj, study, 4, 'lng_noi_again_101320', maskname)

disp('Done!')
