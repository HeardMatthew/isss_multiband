%% isss_multi_process
% Builds models for each subject
% CHANGELOG (MM/DD/YY)
% 11/13/17  Initialized file -- MH
% 02/28/17  Running an analysis to see which acquisition window is best
% 05/30/18  Supercomputer analysis. Forking from original. 
% 02/12/20  Forked for YA_FST. Another timing error found, should be
%   resolved with how I calculated onsets. Adding universal functions.
% 02/17/20  Modifying 1st level design so it specifies design with 1 run
%   for visual inspection. Re-preprocessed with 3 x 3 x 3.5mm voxels. 
% 02/18/20  Dropping physio.
% 02/21/20  New timing scheme, dropped first 5 images
% 03/10/20  New design with CONN regressors
% 04/01/20  Cloned for isss_multiband
% 04/27/20  New design, noise, silent, and language events
% 09/11/20  New analysis, second level without accuracy regressor. Easiest
%   to just throw something together. 
% 09/24/20  Comparing whole brain SPA

clc; clearvars 
do_first  = 0; 
do_second = 1; 

%% Pathing
dir_batch = pwd;
dir_process = fullfile(dir_batch, 'process');
isss_multi_params

if do_first
    for ii = 1:length(subj)
        thissubj = subj(ii);
        disp(['Starting to process ' thissubj.name '.'])

        %% Extract timing 
%         disp('Extracting timing from behav for lang...')
%         cd(dir_process); extract_timing_all(thissubj, study, 1, 6); % HYBRID
%         cd(dir_process); extract_timing_all(thissubj, study, 2, 6); % ISSS
%         cd(dir_process); extract_timing_all(thissubj, study, 3, 6); % MULTIBAND
%         disp('Done!')

        %% Specify and estimate GLM using FIR 
%         disp(['Specifying 1st level GLM for subject ' thissubj.name '.'])
%         % subj, study, design, scan, first run, interactive
%         
%         cd(dir_process); clear_spm_mat(thissubj, study, 6, 1) % HYBRID
%         cd(dir_process); spec_est_GLM(thissubj,  study, 6, 1, 0, 0)
% 
%         cd(dir_process); clear_spm_mat(thissubj, study, 6, 2) % ISSS
%         cd(dir_process); spec_est_GLM(thissubj,  study, 6, 2, 0, 0)
%         
%         cd(dir_process); clear_spm_mat(thissubj, study, 6, 3) % MULTIBAND
%         cd(dir_process); spec_est_GLM(thissubj,  study, 6, 3, 0, 0)
% 
%         disp('Done!')

        %% Build contrasts
        disp(['Building contrasts for subject ' thissubj.name '.'])
        cd(dir_process); build_contrasts(thissubj, study, 6, 1) % HYBRID
        cd(dir_process); build_contrasts(thissubj, study, 6, 2) % ISSS
        cd(dir_process); build_contrasts(thissubj, study, 6, 3) % MULTIBAND
        % subj, study, design, scan
    %     cd(dir_process); build_contrasts_window(thissubj, study, 6)
        disp('Done!')

        %% Calculate AUE 
        disp('Calculating AUE...')
        cd(dir_process); SPA_calculate(thissubj,  study, 6, 1) % HYBRID
        cd(dir_process); SPA_manipulate(thissubj, study, 6, 1)

        cd(dir_process); SPA_calculate(thissubj,  study, 6, 2) % ISSS
        cd(dir_process); SPA_manipulate(thissubj, study, 6, 2)

        cd(dir_process); SPA_calculate(thissubj,  study, 6, 3) % MULTIBAND
        cd(dir_process); SPA_manipulate(thissubj, study, 6, 3)
        disp('Done!')  

        %% Check 1st level results
%         disp('Checking results...')
%         cd(dir_process); results_report(thissubj, study, 4, 1) % HYBRID
%         cd(dir_process); results_report(thissubj, study, 4, 2) % ISSS
%         cd(dir_process); results_report(thissubj, study, 4, 3) % MULTIBAND
%         disp('Done!')    

    end
    
end

if do_second
    %% Create second-level accuracy regressor
%     disp('Creating 2nd level accuracy regressor')
%     cd(dir_process); group_acc_reg(subj, study, 1) % HYBRID
%     cd(dir_process); group_acc_reg(subj, study, 2) % ISSS
%     cd(dir_process); group_acc_reg(subj, study, 3) % MULTIBAND
%     % subjects, study, scan
%     disp('Done!')

    %% Build second-level GLMs
%     disp('Building, evaluating, and generating results for second-level GLMs...')
% %     subj = subj(2:end); % manually drop first one!
%     cd(dir_process); second_level_AUE(subj, study, 4, 1, 1, 0) % HYBRID
%     cd(dir_process); second_level_AUE(subj, study, 4, 2, 1, 0) % ISSS
%     cd(dir_process); second_level_AUE(subj, study, 4, 3, 1, 0) % MULTIBAND
%     % subjects, study, design, scan, create, plot
%     disp('Done!')
    
%     subj = subj(2:end); % manually drop first one!
%     cd(dir_process); second_level_AUE_noacc(subj, study, 5, 1, 1, 0) % HYBRID
%     cd(dir_process); second_level_AUE_noacc(subj, study, 5, 2, 1, 0) % ISSS
%     cd(dir_process); second_level_AUE_noacc(subj, study, 5, 3, 1, 0) % MULTIBAND
%     % subjects, study, design, scan, create, plot
%     disp('Done!')

%     %% Build paired t-test
%     disp('Running paired t-tests...')
%     cd(dir_process); paired_t_AUE(subj, study, 6, 'hybrid_isss',      1, 0)
%     cd(dir_process); paired_t_AUE(subj, study, 6, 'hybrid_multiband', 1, 0)
%     disp('Done!')
% 
%     %% Check SNR
% %     disp('Running SNR ANOVA...')
%     cd(dir_process); combine_snr(subj, study, 6, 1) % HYBRID
%     cd(dir_process); combine_snr(subj, study, 6, 2) % ISSS
%     cd(dir_process); combine_snr(subj, study, 6, 3) % MULTIBAND
%     cd(dir_process); snr_ANOVA(subj, study)
%     disp('Done!')

    %% Compare whole brain SPA
%     disp('Comparing whole brain SPA...')
%     cd(dir_process); compare_wb_spa(subj, study, 4) % all three scans
%     disp('Done!')
    
    %% Anatomical compartments
    disp('Anatomical SPA comparison...')
    cd(dir_process); get_anatomy_spa(subj, study, 4) % all three scans
    disp('Done!')
end
