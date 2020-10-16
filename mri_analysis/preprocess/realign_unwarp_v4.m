%% realign_unwarp
% Realigns and unwarps functional images using SPM. 
% 
% MM/DD/YY: Changelog
% 02/10/20: Forked from isss_multi. 
% 02/21/20: Version 2, drops the first five images
% 03/01/20: Version 3, dropping the first 5 images + first TR of each
%   event.
% 03/26/20: Version 4, forked for hybrid_isss
% 03/30/20: Works for hybrid and isss, requires more tinkering for
%   multiband. 

function realign_unwarp_v4(varargin)
%% Check input
if length(varargin) < 2
    error('Requires (subj, study) input, (ss) is optional!')
end

subj = varargin{1}; 
study = varargin{2}; 
if length(varargin) > 2
    ss = varargin{3};
else
    ss = 1; 
end

if length(varargin) > 3
    dummys = varargin{4};
else
    dummys = ''; 
end

if ~isstruct(subj) || length(subj) ~= 1
    error('Input ("subj") where subj is a SINGLE struct')
end

if ~isstruct(study)
    error('Input ("study") which has experiment info')
end

%% Parameters
scan = study.scan(ss);  
numruns = subj.runs(ss); 

dont_rename_these = {'DUMMY', 'SBRef', 'DROPFIRST', 'DUMMYEND', 'mean', 'SILENT'}; 

%% Pathing
cd ..
batch = fullfile(pwd, 'matlabbatch', 'SPM_realign_unwarp.mat'); 
load(batch)

dir_subj = fullfile(study.path, 'data', subj.name);

dir_reg  = fullfile(dir_subj, 'reg');
dir_func = fullfile(dir_subj, 'FUNCTIONAL'); 
dir_realign  = fullfile(dir_subj, 'realign', [scan.runname '1']);
dir_ps = fullfile(dir_subj, 'ps');
dir_subj_batch = fullfile(dir_subj, 'batch'); 

%% Re-insert dummy scans
target = fullfile(dir_func, '*.nii'); 
funcfiles = dir(target); funcfiles = {funcfiles(:).name}';
these = contains(funcfiles, scan.runname) & ~contains(funcfiles, ... 
    dont_rename_these(2:end)); 
% Drop renamed, mean, etc. files, but grab dummy scans
these = these & cellfun(@isempty, regexp(funcfiles, '^[swu]')); % preprocessed files
funcfiles = funcfiles(these); 

renamethese = cell(numruns, scan.first+scan.epis);
for rr = 1:numruns
    tag = [scan.runname num2str(rr) '_']; 
    idx = 1; 
    for ii = 2:scan.first+scan.epis
        if ii < 10
            renamethese{rr, idx} = ['DUMMY' tag '0000' num2str(ii)]; 
        elseif ii < 100
            renamethese{rr, idx} = ['DUMMY' tag '000' num2str(ii)]; 
        else
            error('Unexpected scan number!')
        end

        idx = idx + 1; 
    end

end

these = false(length(funcfiles), 1);
for rr = 1:numruns
    for ii = 1:(length(renamethese) - 1) % minus the first one...
        these = these | contains(funcfiles, renamethese{rr, ii}); 
    end
end

if ~any(these) 
    warning('Not renaming any extra EPIs!')
else
    files = funcfiles(these);
    newname = cellfun((@(x) x(6:end)), files, 'UniformOutput', false); 
    files = fullfile(dir_func, files);
    newname = fullfile(dir_func, newname); 
    for ii = 1:length(files)
        movefile(files{ii}, newname{ii})
    end

end

disp('Done un-naming first dummy scans!')

%% Drop the first few images
% The first image is a dummy scan. The next few are not associated with
% any stimuli.
target = fullfile(dir_func, '*.nii'); 
funcfiles = dir(target); funcfiles = {funcfiles(:).name}';
these = contains(funcfiles, scan.runname) & ~contains(funcfiles, ... 
    dont_rename_these); % renamed & mean files
these = these & cellfun(@isempty, regexp(funcfiles, '^[swu]')); % preprocessed files
funcfiles = funcfiles(these); 

renamethese = cell(numruns, scan.first);
for rr = 1:numruns
    tag = [scan.runname num2str(rr) '_']; 
    for ii = 1:scan.first
        if ii < 10
            renamethese{rr, ii} = [tag '0000' num2str(ii)]; 
        elseif ii < 100
            renamethese{rr, ii} = [tag '000' num2str(ii)]; 
        else
            error('Unexpected scan number!')
        end

    end

end

these = false(length(funcfiles), 1);
for rr = 1:numruns
    for ii = 1:length(renamethese)
        these = these | contains(funcfiles, renamethese{rr, ii}); 
    end
end

if ~any(these) 
    warning('Not renaming any extra EPIs!')
else
    files = funcfiles(these);
    newname = cellfun((@(x) ['DUMMY' x]), files, 'UniformOutput', false);
    files = fullfile(dir_func, files);
    newname = fullfile(dir_func, newname); 
    for ii = 1:length(files)
        movefile(files{ii}, newname{ii})
    end

end

disp('Done renaming first dummy scans!')

%% Rename scans taken during "silent" period (multiband only)
if strcmp(scan.runname, 'multiband_run')
    target = fullfile(dir_func, '*.nii'); 
    funcfiles = dir(target); funcfiles = {funcfiles(:).name}';
    these = contains(funcfiles, scan.runname) & ~contains(funcfiles, ... 
        dont_rename_these); % renamed & mean files
    these = these & cellfun(@isempty, regexp(funcfiles, '^[swu]')); % preprocessed files
    funcfiles = funcfiles(these); 

    targets = sort([scan.first+scan.epis+1 : scan.epis + scan.silence : scan.numscans, ... 
                    scan.first+scan.epis+2 : scan.epis + scan.silence : scan.numscans, ...
                    scan.first+scan.epis+3 : scan.epis + scan.silence : scan.numscans, ...
                    scan.first+scan.epis+4 : scan.epis + scan.silence : scan.numscans, ...
                    ]); % only cheats a little...
    
    renamethese = cell(numruns, length(targets)); 
    for rr = 1:numruns
        for ii = 1:length(targets)
             if targets(ii) < 10
                 renamethese{rr, ii} = ['0000' num2str(targets(ii))];
             elseif targets(ii) < 100
                 renamethese{rr, ii} = ['000' num2str(targets(ii))];
             elseif targets(ii) < 1000
                 renamethese{rr, ii} = ['00' num2str(targets(ii))];
             else
                 error('Too large a number!')
             end

        end

    end

    these = false(length(funcfiles), 1);
    for rr = 1:numruns
        for ii = 1:length(renamethese)
            these = these | contains(funcfiles, renamethese{rr, ii}); 
        end
    end

    if ~any(these) 
        warning('Not renaming any extra EPIs!')
    else
        files = funcfiles(these);
        newname = cellfun((@(x) ['SILENT' x]), files, 'UniformOutput', false);
        files = fullfile(dir_func, files);
        newname = fullfile(dir_func, newname); 
        for ii = 1:length(files)
            movefile(files{ii}, newname{ii})
        end

    end

    disp('Done renaming silent scans!') 
end

%% Drop the first couple TRs of each event
target = fullfile(dir_func, '*.nii'); 
funcfiles = dir(target); funcfiles = {funcfiles(:).name}';
these = contains(funcfiles, scan.runname) & ~contains(funcfiles, ... 
    dont_rename_these); 
these = these & cellfun(@isempty, regexp(funcfiles, '^[swu]')); % preprocessed files
funcfiles = funcfiles(these); 

if strcmp(scan.runname, 'hybrid_run')
    targets = sort([scan.first+1 : scan.epis : scan.numscans, scan.first+2 : scan.epis : scan.numscans]); 
elseif strcmp(scan.runname, 'multiband_run')
    targets = sort([scan.first+1 : scan.epis+scan.silence : scan.numscans, ...
                    scan.first+2 : scan.epis+scan.silence : scan.numscans]); 
elseif strcmp(scan.runname, 'isss_run')
    targets = scan.first+1 : scan.epis : scan.numscans; 
end

lastevent = targets(end) + 1; % use this later to drop the last, empty block
renamethese = cell(numruns, length(targets)); 
for rr = 1:numruns
    for ii = 1:length(targets)
         if targets(ii) < 10
             renamethese{rr, ii} = ['0000' num2str(targets(ii))];
         elseif targets(ii) < 100
             renamethese{rr, ii} = ['000' num2str(targets(ii))];
         elseif targets(ii) < 1000
             renamethese{rr, ii} = ['00' num2str(targets(ii))];
         else
             error('Too large a number!')
         end

    end

end

these = false(length(funcfiles), 1);
for ii = 1:length(renamethese)
    these = these | contains(funcfiles, renamethese{1, ii}); 
end

if ~any(these) 
    warning('Not renaming any extra EPIs!')
else
    files = funcfiles(these);
    newname = cellfun((@(x) ['DROPFIRST' x]), files, 'UniformOutput', false);
    files = fullfile(dir_func, files);
    newname = fullfile(dir_func, newname); 
    for ii = 1:length(files)
        movefile(files{ii}, newname{ii})
    end
    
end

disp('Done renaming first TRs of each event!')

%% Drop the last (empty) run
target = fullfile(dir_func, '*.nii'); 
funcfiles = dir(target); funcfiles = {funcfiles(:).name}';
these = contains(funcfiles, scan.runname) & ~contains(funcfiles, ... 
    dont_rename_these); 
funcfiles = funcfiles(these); 

dropend = lastevent:scan.numscans; 

renamethese = cell(numruns, length(dropend)); 
for rr = 1:numruns
    for ii = 1:length(dropend)
         if dropend(ii) < 10
             renamethese{rr, ii} = ['0000' num2str(dropend(ii))];
         elseif dropend(ii) < 100
             renamethese{rr, ii} = ['000' num2str(dropend(ii))];
         elseif dropend(ii) < 1000
             renamethese{rr, ii} = ['00' num2str(dropend(ii))];
         else
             error('Too large a number!')
         end

    end

end

these = false(length(funcfiles), 1);
for ii = 1:length(renamethese)
    these = these | contains(funcfiles, renamethese{1, ii}); 
end

if ~any(these) 
    warning('Not renaming any extra EPIs!')
else
    files = funcfiles(these);
    newname = cellfun((@(x) ['DUMMYEND' x]), files, 'UniformOutput', false);
    files = fullfile(dir_func, files);
    newname = fullfile(dir_func, newname); 
    for ii = 1:length(files)
        movefile(files{ii}, newname{ii})
    end
    
end

disp('Done renaming dummy end scans!')

for rr = 1:numruns
    %% Load names and put into mattlabbatch
    [boldFiles, ~] = spm_select('List', dir_func, ['^' scan.runname num2str(rr) '_00*.*\.nii$']);
    boldFiles = [repmat([dir_func filesep], length(boldFiles),1), boldFiles];
    boldFiles = cellstr(boldFiles);
    matlabbatch{1}.spm.spatial.realignunwarp.data(rr).scans = boldFiles;
    
    [vdmFiles, ~] = spm_select('List', dir_realign, '^vdm5_*.*\.nii$');
    vdmFiles = cellstr([dir_realign filesep vdmFiles]);
    matlabbatch{1}.spm.spatial.realignunwarp.data(rr).pmscan = vdmFiles;
    
    if strcmp(matlabbatch{1}.spm.spatial.realignunwarp.data(rr).scans, '<UNDEFINED>')
        matlabbatch{1}.spm.spatial.realignunwarp.data(rr) = [];
    end
    
end

if length(matlabbatch{1}.spm.spatial.realignunwarp.data) > numruns
    matlabbatch{1}.spm.spatial.realignunwarp.data(numruns+1:end) = [];
end

filename = fullfile(dir_subj_batch, ['realign_unwarp_v4_' scan.runname '.mat']);
save(filename, 'matlabbatch') 

%% Prepare SPM
disp('Loading SPM...')
close all; 
spm('defaults', 'FMRI');
spm_jobman('initcfg');
spm_figure('GetWin','Graphics'); % Thanks Guillaume
disp('Done!')

%% Run job
disp('Starting realign_unwarp...')
spm_jobman('run',matlabbatch);

%% Move regressors and motion graph to their proper directory
target = fullfile(dir_func, '*.txt'); 
regs = dir(target);
regNames = fullfile(dir_func, {regs(:).name}');
for ii = 1:length(regNames)
    movefile(regNames{ii}, dir_reg);
end

target = fullfile(pwd, '*.png');
psfiles = dir(target); 
psNames = fullfile(pwd, {psfiles(:).name}'); 
for ii = 1:length(psfiles)
    if contains(dummys, 'first')
        target = fullfile(dir_ps, ['dummyfirst_realign_unwarp_v4_' scan.runname '_' num2str(ii) '.png']); 
    else
        target = fullfile(dir_ps, ['realign_unwarp_v4_' scan.runname '_' num2str(ii) '.png']); 
    end
    movefile(psNames{ii}, target)
end

end