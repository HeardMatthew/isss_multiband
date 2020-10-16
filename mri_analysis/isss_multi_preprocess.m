%% preprocess.m
% Preprocesses images for analysis. 
% 
% MM/DD/YY: Changelog
% 02/10/20: Forked from isss_multi. 
%   + Fieldmaps are created using FSL (batch/unwarp.sh) and can be run from
%     here or from bash. 
% 03/09/20: Subjects 5976YL and 5977YL ran without issues. 
% 03/26/20: Cloned for hybrid_isss! 
% 03/30/20: Realign and unwarp done. Renamed any extra scans!
% 04/28/20: Re-preprocessed data, including first dummy scans
% 04/29/20: Back to old design. Making MVPA dataset with 3mm smoothing. 

clc; clearvars

%% Pathing
dir_batch = pwd;
dir_preprocess = fullfile(dir_batch, 'preprocess'); 
dir_snr = fullfile(dir_batch, 'Noise_script'); 
isss_multi_params

cd ..
dir_data = fullfile(pwd, 'data'); 

tic
for ii = 3:length(subj)
    %% Setup for this subject
    thissubj = subj(ii);
    dir_subj = fullfile(dir_data, thissubj.name); 
    disp(['Preprocessing subj ' thissubj.name '...'])
    
    %% Purge and archive the old stuff (BE CAREFUL WITH THIS ONE!)
% %     disp('Archiving old files!')
% %     cd(dir_preprocess); purge_and_archive(thissubj, study)
% %     disp('Done!')
    
    %% Rename existing FUNC_GLM to allow for a new set of images
%     disp('Renaming old FUNC_GLM...')
%     cd(dir_preprocess); rename_FUNC_GLM(thissubj, study)
%     disp('Done!')

    %% Delete existing preprocessed data to prevent confusion
    disp('Purging the old stuff...')
    cd(dir_preprocess); purge_the_old_stuff(thissubj, study)
    disp('Done!')

    %% Unwarp and realignment 
    % Using previously created fieldmaps
    disp('Unwarping and realigning data')
    cd(dir_preprocess); realign_unwarp_v4(thissubj, study, 1) % HYBRID
    cd(dir_preprocess); realign_unwarp_v4(thissubj, study, 2) % ISSS
    cd(dir_preprocess); realign_unwarp_v4(thissubj, study, 3) % MULTIBAND
    disp('Done!')

    %% Coregistration
    disp('Coregistration begins')
    cd(dir_preprocess); coregister(thissubj, study, 1) % HYBRID
    cd(dir_preprocess); coregister(thissubj, study, 2) % ISSS
    cd(dir_preprocess); coregister(thissubj, study, 3) % MULTIBAND
    disp('Done coregistering!')
     
    %% Normalization 
    disp('Normalizing to MNI-space')
    cd(dir_preprocess); normalize(thissubj, study, 1) % HYBRID
    cd(dir_preprocess); normalize(thissubj, study, 2) % ISSS
    cd(dir_preprocess); normalize(thissubj, study, 3) % MULTIBAND
    disp('Done normalizing!')

    %% Smoothing
    disp('Smoothing data now')
    cd(dir_preprocess); smooth(thissubj, study, 1, 'MVPA') % HYBRID
    cd(dir_preprocess); smooth(thissubj, study, 2, 'MVPA') % ISSS
    cd(dir_preprocess); smooth(thissubj, study, 3, 'MVPA') % MULTIBAND
    disp('Done smoothing!')

    %% SNR
%     disp('Running SNR scripts')
%     cd(dir_snr); snr_sd_v4_isss_multi(thissubj, study, 1) % HYBRID
%     cd(dir_snr); snr_sd_v4_isss_multi(thissubj, study, 2) % ISSS
%     cd(dir_snr); snr_sd_v4_isss_multi(thissubj, study, 3) % MULTIBAND
%     disp('Done with SNR!')

end

disp('Batch complete!')
toc
