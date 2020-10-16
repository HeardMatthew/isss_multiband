%%

function purge_the_old_stuff(subj, study)
%% Check input
if ~isstruct(subj) || length(subj) ~= 1
    error('Input ("subj") where subj is a SINGLE struct')
end

if ~isstruct(study)
    error('Input ("study") which has experiment info')
end

%% Pathing
dir_subj = fullfile(study.path, 'data', subj.name);
dir_func = fullfile(dir_subj, 'FUNCTIONAL'); 

%% Find the old stuff
target = fullfile(dir_func, '*nii'); 
files = dir(target); files = {files(:).name}; 
these = ~cellfun(@isempty, regexp(files, '^[swu]')); % preprocessed files
delete_these = fullfile(dir_func, files(these)); 
nfiles = length(delete_these); 
for nf = 1:nfiles
    delete(delete_these{nf})
end


end