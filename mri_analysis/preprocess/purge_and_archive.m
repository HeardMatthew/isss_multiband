%% purge_and_archive
% Deletes old files and archives important ones

function purge_and_archive(subj, study)
%% Check input
if ~isstruct(subj) || length(subj) ~= 1
    error('Input ("subj") where subj is a SINGLE struct')
end

if ~isstruct(study)
    error('Input ("study") which has experiment info')
end

%% Pathing
dir_subj = fullfile(study.path, 'data', subj.name);

%% Parameters
folders_archive = { ...
    'batch', 'design', 'dicm2nii_out_112018', 'ps', 'reg', 'SNR', ...
    }; 
folders_purge = { ...
    'FUNC_GLM', 'FUNC_GLM_DOWNSAMPLE', 'FUNC_MVPA', 'FUNCTIONAL_DOWNSAMPLE', 'FUNC_MVPA_DOWNSAMPLE', 'FUNC_HRF', ...
    }; 

scans_purge = { ...
    'mean', 'u', ...
    }; 

%% Archive requested folders
% disp('Archiving folders...')
% for ff = 1:length(folders_archive)
%     target = fullfile(dir_subj, folders_archive{ff}); 
%     newname = fullfile(dir_subj, [folders_archive{ff}, '_old']); 
%     if exist(newname, 'file')
%         warning(['Already exists: ' newname])
%     elseif ~exist(target, 'file')
%         warning(['Cannot find: ' target])
%     else
%         zip(newname, target)
%         disp(newname)
%         rmdir(target, 's')
%         mkdir(target)
%     end
%     
% end
% 
% disp('Done archiving!')

%% Purge requested folders
% disp('Purging folders...')
% for ff = 1:length(folders_purge)
%     target = fullfile(dir_subj, folders_purge{ff}); 
%     if exist(target, 'file')
%         rmdir(target, 's')
%         mkdir(target)
%         disp(target)
%     else
%         warning(['Cannot find: ' target])
%     end
%     
% end
% 
% disp('Done purging!')

%% Purge requested scans
disp('Purging scans...')
files = dir(fullfile(dir_subj, 'FUNCTIONAL')); files = {files(:).name};
for ss = 1:length(scans_purge)
    these = ~cellfun(@isempty, regexp(files, ['^' scans_purge{ss}])); 
%     these = contains(files, scans_purge{ss}); 
    if ~any(these)
        warning('No scans to purge!')
    else
        targets = fullfile(dir_subj, 'FUNCTIONAL', files(these)); 
        for tt = 1:length(targets)
            delete(targets{tt});
        end
            
    end
    
end

disp('Done purging scans!')

end