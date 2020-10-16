%% rename_FUNC_GLM
% 
% 

function rename_FUNC_GLM(subj, study)

%% Pathing
dir_subj = fullfile(study.path, 'data', subj.name);

dir_func_GLM = fullfile(dir_subj, 'FUNC_GLM'); 
dir_renamed  = fullfile(dir_subj, 'FUNC_GLM_v1'); 

%% The lifting
if ~exist(dir_renamed, 'file')
    movefile(dir_func_GLM, dir_renamed)
    mkdir(dir_func_GLM)
else
    warning('Folder already exists!')
end

end