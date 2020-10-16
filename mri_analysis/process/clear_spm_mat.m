%% clear_spm_mat
% Deletes SPM.mat in the specified directory. 
% CHANGELOG
% 02/12/20 - File initialized. MJH
% 02/17/20 - Now looks for 1run designs. MJH

function clear_spm_mat(subj, study, dd, ss)
%% Check input
if ~isstruct(subj) || length(subj) ~= 1
    error('Input ("subj") where subj is a SINGLE struct')
end

if ~isstruct(study)
    error('Input ("study") which has experiment info')
end

if ~isnumeric(dd)
    error('Input ("dd") which specifies which design should be cleared')
end

if ~isnumeric(ss)
    error('Inout "ss" which specifies which scan!')
end

scan = study.scan(ss); 

%% Pathing
dir_start = pwd;
dir_design = fullfile(study.path, 'data', subj.name, 'design'); 

for ii = 1:length(dd)
    name = [scan.runname(1:end-4) '_' study.design(dd(ii)).name]; 
    dir_target = fullfile(dir_design, name); 
    try
        cd(dir_target)
        delete SPM.mat
        cd(dir_start)
    catch
        disp(['Folder ' dir_target ' does not exist!'])
    end    
    
    name = [scan.runname(1:end-4) '_' study.design(dd(ii)).name '_1run']; 
    dir_target = fullfile(dir_design, name); 
    try
        cd(dir_target)
        delete SPM.mat
        cd(dir_start)
    catch
        disp(['Folder ' dir_target ' does not exist!'])
    end    
    
end

end
