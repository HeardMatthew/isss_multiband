%% isss_multi_MVPA.m
% Batch script to run all MVPA analysis on supercomputer. 
% 
% MM/DD/YY -- CHANGELOG
% 06/13/18 -- File initialized. Cloned to supercomputer. 
% 02/14/20 -- Cloned for YA_FST. Making scripts universal. MJH
% 02/18/20 -- Dropping physio and rerunning
% 02/20/20 -- New strategy to drop physio, rerunning. 
% 02/25/20 -- Model with all betas, wu- data, but only investigating TRs 2
%   through 4. 
% 03/05/20 -- New experiment structure. 
% 04/24/20 -- Copied for isss_multi
% 04/30/20 -- Updating for isss_multi. Using slightly smoothed data. 

clearvars; clc; 
do_first  = 0; 
do_second = 1; 

%% Pathing and parameters
dir_batch = pwd;
dir_MVPA  = fullfile(dir_batch, 'MVPA'); 
dir_process = fullfile(dir_batch, 'process'); 
dir_preprocess = fullfile(dir_batch, 'preprocess'); 
isss_multi_params

tic

if do_first
    try parpool(28); catch; delete(gcp('nocreate')); parpool(28); end
    
    for ii = 2:length(subj) % subject 1 disqualified from MVPA, only 1 run!
        %% Subject-specific parameters
        thissubj = subj(ii);
        disp(['MVPA on subject ' thissubj.name '...'])

        %% Make CONN regressors from lightly smoothed data
        % Done!

        %% Specify and estimate GLMs using FIR
%         disp(['Specifying 1st level GLM for subject ' thissubj.name '.'])
%         cd(dir_process); clear_spm_mat(thissubj, study, 8, 1)
%         cd(dir_process); spec_est_GLM(thissubj, study, 8, 1, 0, 0)
%         cd(dir_process); spec_est_GLM(thissubj, study, 8, 2, 0, 0)
%         cd(dir_process); spec_est_GLM(thissubj, study, 8, 3, 0, 0)
%         disp('Done!')

        %% Extract time course information
        disp(['Extracting time course information for ' thissubj.name '...'])
        cd(dir_MVPA); extract_betas_v2(thissubj, study, 8, 1) % HYBRID
        cd(dir_MVPA); extract_betas_v2(thissubj, study, 8, 2) % ISSS
        cd(dir_MVPA); extract_betas_v2(thissubj, study, 8, 3) % MULTIBAND
        disp('Done!')

        %% Create coordinates spheres for all subjects
        disp(['Creating spheres for ' thissubj.name '...'])
        cd(dir_MVPA); coord_sphere(thissubj, study, 8, 1) % HYBRID
        cd(dir_MVPA); coord_sphere(thissubj, study, 8, 2) % ISSS
        cd(dir_MVPA); coord_sphere(thissubj, study, 8, 3) % MULTIBAND
        disp('Done!')

        %% Run the searchlight analysis
        disp(['Running searchlights on ' thissubj.name '...'])
        cd(dir_MVPA); searchlight_LNG_NOI(thissubj, study, 8, 1); % HYBRID
        cd(dir_MVPA); searchlight_LNG_NOI(thissubj, study, 8, 2); % ISSS
        cd(dir_MVPA); searchlight_LNG_NOI(thissubj, study, 8, 3); % MULTIBAND
        disp('Done!')

        %% Make histograms for each subject!
    %     disp(['Making histograms for ' thissubj.name '...'])
    %     cd(dir_MVPA); plot_acc(thissubj, study, 5, 1, 'GNB_LNG_NOI') % HYBRID
    %     cd(dir_MVPA); plot_acc(thissubj, study, 5, 2, 'GNB_LNG_NOI') % ISSS
    %     cd(dir_MVPA); plot_acc(thissubj, study, 5, 3, 'GNB_LNG_NOI') % MULTIBAND
    %     disp('Done!')

    end
    
end

if do_second
    subj = subj(2:end); % drop first subject!
    
    %% Second-level stats
%     cd(dir_MVPA); second_level_MVPA(subj, study, 5, 1, 'GNB_LNG_NOI'); % HYBRID
%     cd(dir_MVPA); second_level_MVPA(subj, study, 5, 2, 'GNB_LNG_NOI'); % ISSS
%     cd(dir_MVPA); second_level_MVPA(subj, study, 5, 3, 'GNB_LNG_NOI'); % MULTI

%     cd(dir_MVPA); second_level_MVPA_noacc(subj, study, 5, 1, 'GNB_LNG_NOI'); % HYBRID
%     cd(dir_MVPA); second_level_MVPA_noacc(subj, study, 5, 2, 'GNB_LNG_NOI'); % ISSS
%     cd(dir_MVPA); second_level_MVPA_noacc(subj, study, 5, 3, 'GNB_LNG_NOI'); % MULTI

    %% Whole brain plots
%     cd(dir_MVPA); get_wb_accuracy(subj, study, 5, 'GNB_LNG_NOI', {'hybrid', 'isss'});
%     cd(dir_MVPA); get_wb_accuracy(subj, study, 5, 'GNB_LNG_NOI', {'hybrid', 'multiband'});
    
    cd(dir_MVPA); get_sig_accuracy(subj, study, 5, 'GNB_LNG_NOI', {'hybrid', 'isss'});
    cd(dir_MVPA); get_sig_accuracy(subj, study, 5, 'GNB_LNG_NOI', {'hybrid', 'multiband'});
end

disp('Completed batch!'); toc
