%% get_anatomy_spa
% Calculate mean, SD, and SER across whole brain MVPA results
% 
% MM/DD/YY -- CHANGELOG
% 05/08/20 -- Initialized. 
% 09/14/20 -- Forcing this to call manuscript dataset because I'm at my
%   wit's end. Added "scans" to input
% 09/24/20 -- Forked to check whole brain SPA values across scans. Manually
%   coding the contrast of interest!
% 09/25/20 -- New version looking at anatomical compartments. 
% 1. Identify good atlases for cortex, basal ganglia, thalamus, and
%   brainstem
% 2. Make brain mask
% 3. Load atlases and screen through brain mask
% 4. Get accuracy for language and noise
% 5. Make figure checking for double dissociation

function get_anatomy_spa(subj, study, dd)
%% Flags
do_create  = 0; 
do_analyze = 0; 
do_plot    = 1;

%% check inputs
if length(subj) < 2
    error('Submit ALL subjects!')
end

if ~isstruct(subj) || ~isstruct(study)
    error('subj and study are both struct!')
end

if ~isnumeric(dd)
    error('dd specify number of design!')
end

%% Pathing
dir_process = pwd; 
cd ..
dir_batch = pwd; 

% Shame on me for using direct paths. 
dir_masks = 'C:\Users\heard\Documents\MATLAB\mri_data\atlases\HarvardOxford'; 
dir_roi   = 'C:\Users\heard\Documents\MATLAB\mri_data\isss_multi\roi'; % lazy
dir_docs  = 'C:\Users\heard\Documents\MATLAB\mri_data\isss_multi\docs\092820_anatomical_plots'; 

%% Parameters
scans       = {'hybrid', 'isss', 'multiband'}; 
design      = study.design(dd);
conditions  = {'LNG', 'NOI'}; 
regions     = {'cortex', 'subcortex', 'thalamus', 'brainstem'}; 
regions_idx = {[2 13], [5, 6, 7, 9, 10, 11, 16, 17, 18, 19, 20, 21], ...
    [4, 15], 8}; 
    % Cortex - 2 (LH) and 13 (RH) 
    % Subcortex - 5, 16 (Caudate), 6, 17 (putamen), 7, 18 (GP), 9, 19
    %   (hippocampus), 11, 21 (NAcc), 10, 20 (amygdala)
    % Thalamus - 4 (L) and 15 (R)
    % Brainstem - 8 
    
numscans   = length(scans); 
numsubjs   = length(subj); 
numcons    = length(conditions);
numregions = length(regions);

%% Load second-level whole brain masks
if do_create
    disp('Creating brain masks...')

    % Preallocate
    dir_designs = cell(numscans, 1); 

    % Scan specific template
    disp('Making brain mask across scans')
    for nscans = 1:numscans
        dir_designs{nscans} = fullfile(study.path, 'data', '050520_manuscript', [scans{nscans} '_' design.name], 'LNG_NOI'); 
        fname = fullfile(dir_designs{nscans}, 'mask.nii'); % WHOLE-BRAIN MASK
        Vmask = spm_vol(fname); [mask, xyz] = spm_read_vols(Vmask); 
        mask = logical(mask);
        if nscans == 1
            brain_mask = mask; 
        else
            brain_mask = brain_mask & mask; 
        end

    end

    voxels_in_mask = xyz(:, brain_mask(:)); % useful for debugging
    
    % Combine with Harvard Oxford atlas
    disp('Combining with Harvard Oxford atlas')
    %%% HARD CODING PATH BECAUSE LAZINESS
    fname = fullfile(dir_masks, 'reslice_HarvardOxford-sub-maxprob-thr25-1mm.nii'); 
    V = spm_vol(fname);
    HO = spm_read_vols(V); 
    
    for rr = 1:numregions
        theseidx = regions_idx{rr}; 
        disp(regions{rr})
        HO_regions = reshape(any(HO(:) == theseidx, 2), size(brain_mask)); 
        HO_regions = HO_regions & brain_mask; 
        
        % Save new mask
        V.fname = fullfile(dir_masks, ['HO_' regions{rr} '_isss_multi.nii']); 
        spm_write_vol(V, HO_regions); 
    end
    
end

%% Get accuracy per subject
if do_analyze
    
    for rr = 1:numregions
        roifile = fullfile(dir_roi, ['HO_' regions{rr} '_isss_multi.nii']);
        V = spm_vol(roifile); 
        brain_mask = logical(spm_read_vols(V)); 
        
        % Preallocate matrix for average and standard error
        SPA_means = nan(numscans, numcons);  
        SPA_serrs = nan(numscans, numcons); 
        
        for ncons = 1:numcons
            % Preallocate
%             SPA_avg_all = []; 
%             SPA_std_all = []; 
%             SPA_table = nan(numsubjs, numscans); 
            
            for nscans = 1:numscans
                SPA = zeros(numsubjs, 1); 
                
                for nsubjs = 1:numsubjs
                    % Subject-specific path and parameters
                    thissubj = subj(nsubjs); 
                    disp(thissubj.name)
                    dir_design = fullfile(study.path, 'data', thissubj.name, 'design', [scans{nscans} '_' design.name]); 
                
                    % Load SPA for subjects 
                    fname = fullfile(dir_design, ['AUE_' conditions{ncons} '.nii']); 
                    V = spm_vol(fname); 
                    data = spm_read_vols(V); 
                    SPA(nsubjs) = mean(data(brain_mask)); 
                end

%             SPA_table(:, nscan) = SPA; 

                %% Average across subjects and save
                SPA_means(nscans, ncons) = mean(SPA);
                SPA_serrs(nscans, ncons) = std(SPA)/sqrt(numsubjs); 

            end
            
        end

        %% Save results
        scanname = scans'; 
        LNG = SPA_means(:, 1); 
        NOI = SPA_means(:, 2); 
        data_means = table(scanname, LNG, NOI); 
        
        LNG = SPA_serrs(:, 1); 
        NOI = SPA_serrs(:, 2); 
        data_serrs = table(scanname, LNG, NOI); 
        
        % scannames = strjoin(scans, '_'); 
        fname = fullfile(dir_batch, ['SPA_LNG_NOI_' regions{rr} '_' date '.mat']); 
        save(fname, 'data', 'data_means', 'data_serrs')
    end
    
end

if do_plot
    close all
    for rr = 1:numregions
        % Load data
        fname = fullfile(dir_batch, ['SPA_LNG_NOI_' regions{rr} '_' date '.mat']); 
        load(fname)
        
        means = [data_means.LNG([1, 3]), data_means.NOI([1, 3])]'; 
        serrs = [data_serrs.LNG([1, 3]), data_serrs.NOI([1, 3])]'; 
        % Make figure
        figure(rr)
        H = bar(means); 
        title(regions{rr})  
        hold on
        H(1).FaceColor = [1 0 0];
        H(2).FaceColor = [0 0 0];
        xlabel('Condition'); 
        xticklabels({'LNG', 'NOI'})
        ylabel('SPA value')
        
        %%% MATHWORKS
        % Find the number of groups and the number of bars in each group
        ngroups = size(means, 1);
        nbars = size(means, 2);
        % Calculate the width for each bar group
        groupwidth = min(0.8, nbars/(nbars + 1.5));
        % Set the position of each error bar in the centre of the main bar
        % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
        for i = 1:nbars
            % Calculate center of each bar
            x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
            E = errorbar(x, means(:,i), serrs(:,i), 'k', 'linestyle', 'none');
            E.Color = [0.5 0.5 0.5];
        end
        %%%
        
        legend({'Hybrid', 'Multiband'}, 'Location', 'northeast')
        
        % Save data
        fname = fullfile(dir_docs, ['SPA_plot_' regions{rr}]); 
        saveas(gcf, fname, 'png')
    end
    
end
    
end