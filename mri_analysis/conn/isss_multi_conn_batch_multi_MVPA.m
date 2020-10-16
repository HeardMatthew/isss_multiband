%% hybrid_isss_conn_batch
% Uses CONN to perform aCompCorr, creates regressors for WM and CSF signal
% fluctuations. These regressors are used to remove noise from GLMs. 
% Author: Matthew Heard
% 
% MM/DD/YY -- CHANGELOG
% 03/16/20 -- Changelog started. 
% 04/09/20 -- Cloned for hybrid_isss. Going to run one script for each
%   design. "Basic" is now done! Loaded in structural and functional
%   images, after much consternation. 
% 04/09/20 -- Finished preprocessing anatomical images. Now ready to run
%   the regressor code. 

% Sample CONN Script for the OSU Workshop
% Created by Andrew Jahn, University of Michigan, 02.27.2020
% Adapted from Alfonso Nieto-Castanon's script, www.alfnie.com 
clearvars; clc; 
home = pwd; cd ..; 
isss_multi_params
cd(home)

load_batch = 0; 
create_reg = 1; 

%% PREALLOCATED TEMPLATE AND BATCH INFO (CHANGE THIS STUFF)
ss = 3; % MULTIBAND
conn_batch = 'conn_multi_physio.mat'; % just borrow this to load ROI info
conn_name = 'conn_multi_physio_MVPA'; 
scan = study.scan(ss);

if load_batch
    load(conn_batch) % a preallocated template
    fname = fullfile(pwd, [conn_name '.mat']); 

    %% Parameters
    NSUBJECTS = length(subj);
    scan = study.scan(ss); 

    %% Set up batch path
    CONN_x.filename = fname; 

    folders = fields(CONN_x.folders); 
    folders = folders(2:end); % not adjusting ROIs
    for ii = 1:length(folders)
        temp = CONN_x.folders.(folders{ii}); 
        temp = strsplit(temp, '/'); 
        temp{9} = conn_name; 
        temp = strjoin(temp, '/'); 
        CONN_x.folders.(folders{ii}) = temp; 
    end

    % CONN_x.Setup.isnew = 1; % not actually a field?
    CONN_x.Setup.acquisitiontype = 2; % we are modeling as discontinuous data

    CONN_x.Setup.nsubjects = NSUBJECTS;
    CONN_x.Setup.RT = repelem(scan.TR, NSUBJECTS); 
    % CONN_x.Setup.functional = repmat({{}}, [1, NSUBJECTS]);  

    for nsub = 1:length(subj)
        thissubj = subj(nsub); 
        disp(thissubj.name)
        %% Paths for each subject
        dir_subj = fullfile(study.path, 'data', thissubj.name); 
        dir_func_MVPA = fullfile(dir_subj, 'FUNCTIONAL'); 
        dir_anat = fullfile(dir_subj, 'ANATOMICAL'); 
        dir_func_GLM = fullfile(dir_subj, 'FUNC_GLM'); 

        %% Parameters
        nsessions = thissubj.runs(ss); 
        CONN_x.Setup.nsessions(nsub) = nsessions; 

        %% Anatomical files (already loaded in template)
%         % May need to preprocess anatomical images in this script...
%         target = fullfile(dir_anat, study.anat); 
%         file_anat = dir(target); 
%         file_anat = fullfile(dir_anat, file_anat(:).name); 
% 
%         % The .structural field has 3 cells inside another cell: 
%         % 1) Full path to the file. 
%         % 2) Another nested cell:
%         % 2.1) String that says how many files x size of file. I just copied
%         %      and pasted this from the first preallocated file since MPRAGE 
%         %      has same dimensions across all subjects. 
%         % 2.2) Some shorthand version of the original file path? Likely gets
%         %      overwritten once you use the GUI. 
%         % 3) Header info for the file, use spm_vol to load. 
%         % 
%         % I have a hunch that we can skip cells {1}{2} and {1}{3}, but I do not
%         % know for sure. This code works so I ain't touching it. 
%         CONN_x.Setup.structural{nsub}{1}{1} = file_anat;  
%         V = spm_vol(file_anat); 
%         CONN_x.Setup.structural{nsub}{1}{3} = V; 
%         CONN_x.Setup.structural{nsub}{1}{2}{1} = ... 
%             ['[1 file] x [size ' num2str(V.dim) ' ]']; 

        %% Functional data
        % Probably finding the same pattern as anatomical...

        % CONN_x.Setup.functional{nsub}{nses}{those three cells again}
        for nses = 1:nsessions
            target = fullfile(dir_func_MVPA, ... 
                [study.mvpa.prefix scan.runname num2str(nses) '*.nii']); 
            files_func = dir(target); 
            files_func = fullfile(dir_func_MVPA, {files_func(:).name})';
            files_func = char(files_func); 
            CONN_x.Setup.functional{nsub}{nses}{1} = files_func; 

            V = spm_vol(files_func(1, :)); 
            CONN_x.Setup.functional{nsub}{nses}{3} = V; 
            CONN_x.Setup.functional{nsub}{nses}{3}(2) = spm_vol(files_func(end, :)); 

            CONN_x.Setup.functional{nsub}{nses}{2}{1} = ...
                ['[' num2str(size(files_func, 1)) ' files] x [size ' num2str(V.dim) ' ]']; 

            %% Load condition/session info (needs work)
            CONN_x.Setup.conditions.values{nsub}{1}{nses}{1} = 0;
            CONN_x.Setup.conditions.values{nsub}{1}{nses}{2} = Inf;

        end 

    end

    save(fname, 'CONN_x')
end

%% MANUALLY RUN PREPROCESSING HERE!
% Just run structural segmentation. It takes ~1.5 hours for hybrid_isss
% Not neccessary for the MVPA code, all ROI info already in batch. 

%% Create regressors
if create_reg
    dir_data = fullfile(study.path, 'data'); 

    %% This code creates the regressors, saves them as dp_sw*.txt
    conn_module( 'preprocessing', ...
       'steps', {'functional_regression'}, ...
       'reg_names', {'White Matter','CSF'}, ...
       'reg_dimensions',[2, 2], ...
       'reg_skip', true, ...
       'reg_deriv', [0, 0]);

    for ii = 1:length(subj) 
        %% Setup for this subject
        thissubj = subj(ii);
        dir_subj = fullfile(dir_data, thissubj.name); 
        dir_func_MVPA = fullfile(dir_subj, 'FUNCTIONAL'); % location may change
        dir_func_GLM = fullfile(dir_subj, 'FUNC_GLM'); % location may change
        dir_reg  = fullfile(dir_subj, 'reg'); 
        disp(['Writing regressor for subj ' thissubj.name '...'])
        blocks = 1:thissubj.runs(ss); 

        for rr = blocks
            %% Move regressors to regressor directory
            disp(['Block ' num2str(rr)])
            
            if scan.regnum < 10
                regnum = ['0000' num2str(scan.regnum)]; 
            elseif scan.regnum < 100
                regnum = ['000' num2str(scan.regnum)];
            else
                error('too large regnum!')
            end

            fname = fullfile(dir_func_MVPA, ...
                ['dp_' study.mvpa.prefix scan.runname num2str(rr) '_' regnum '.txt']); 
            target = fullfile(dir_reg, ...
                ['dp_MVPA_' study.mvpa.prefix scan.runname num2str(rr) '_' regnum '.txt']); 
            if exist(fname, 'file')
                movefile(fname, target)
            end

            %% Load our new regressors
            T = readtable(target); 

            %% Trim off the first two columns of the data 
            % First two columns appear to be constant and linear terms
            newname = fullfile(dir_reg, ...
                ['dp_trimmed_pc2_MVPA_' study.mvpa.prefix scan.runname num2str(rr) '_' regnum '.txt']); 

            T = T(:, 3:end); 
            writetable(T, newname, 'Delimiter', '\t', 'WriteVariableNames', false)
        end

    end

end

