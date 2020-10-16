%% extract_timing_nowrong
% Loads and processes timing information. Spits out a regressor of
% incorrect scans. 
% Input: dir_subj

% CHANGELOG
% 09/01/17  File inception
% 09/06/17  Worked on timing and duration extraction
% 09/07/17  Added button press condition, accuracy regressor
% 09/08/17  Changed accuracy regressor to a condition, started function
% 09/15/17  Updated to use cells to store onsets, events, and durations
% 09/26/17  Changing to use cells instead of structs. Way easier to handle
%   in SPM batch
% 09/29/17  Combined with convert to scans. 
% 11/08/17  Updated with new naming conventions
% 11/21/17  Went through to try and catch errors. Turned into a new,
%   simpler script. 
% 02/07/20  Forked for YA_FST. Loads in xlsx files. 
% 02/10/20  v1 finished, timing extracted for 3 subjects. 
% 02/13/20  Errors found, trying again. 
% 02/21/20  New timing scheme. 
% 02/29/20  Updated for new version of experiment. New strategy for
%   encoding onsets. 
% 03/02/20  We are dropping the first TR of each event. Updating to reflect
%   this. 
% 03/09/20  Added 'perfect' output to represent any blocks with perfect
%   accuracy. This will be used to select regressors. 
% 04/27/20  New condition: LNG. Combines OR and SR. 
% 04/28/20  New argument firstscan, specifies what should be considered the
%   first scan. 

function onsets = extract_timing_all(varargin)
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
    dd = varargin{4};
else
    error('Need to specify which design!')
end

if ~isstruct(subj) || length(subj) ~= 1
    error('Input ("subj") where subj is a SINGLE struct')
end

if ~isstruct(study)
    error('Input ("study") which has experiment info')
end

%% Parameters and path
onsets = struct( ... % Needs to be updated for each experiment. 
    'NOI', [], ... 
    'SIL', [], ... 
    'LNG', [], ...
    'OR',  [], ... 
    'SR',  []  ...
    ); 

dir_subj  = fullfile(study.path, 'data', subj.name);
dir_reg   = fullfile(dir_subj, 'reg'); 
dir_behav = fullfile(dir_subj, 'behav', 'scan'); 

design = study.design(dd); 

events_each = study.events.each; 
events_all  = study.events.all; 
 
scan = study.scan(ss); 

conNames = fields(onsets);
numCons  = length(conNames); 
runTypes = {'hybrid', 'isss', 'multi'}; 
thisrun  = runTypes{ss}; 

drop = 2/scan.TR; % should automatically identify if dropping 1 or 2 events?
numbetas = scan.epis - drop; 
scan.order = 0:numbetas:numbetas*(events_all + 1) - scan.first; 
% dropped dummy events (beginning and end) as well as and first (two) TR of each run. 
% design(6) now includes 10 of the first dummy events. 
if strcmp(design.dummyscans, 'first')
    scan.order = scan.order + numbetas; 
end

%% Load data for subject
disp('Now loading behavioral data...')
var = dir(fullfile(dir_behav, '*.xlsx'));
files = {var(:).name};
% if any(contains(files, 'splice'))
%     thisone = contains(files, 'splice'); 
% else
    thisone = contains(files, 'lang') & ~contains(files, 'pract') & ~contains(files, 'splice'); 
% end

var = var(thisone); 
disp(['Found ' num2str(length(var)) ' variable files.']) % There might be many from aborted runs

if isempty(var)
    error('No files found!')
elseif length(var) ~= 1
    %% Combining aborted runs 
    disp('Let us begin the splicing!')
    
    for vv = 1:length(var)
        % Load all tables
        oldname = fullfile(var(vv).folder, var(vv).name); 
        
        for rr = 1:subj.runs(ss)
            runname = thisrun; 
            if rr > 1
                runname = [thisrun num2str(rr)]; 
            end
            
            try
                thisT = readtable(oldname, 'Sheet', runname); 
                these = ~isnan(thisT.ActualJitter);
                if ~any(these)
                    disp('This run was cancelled!')
                else
                    newname = fullfile(dir_behav, [var(vv).name(1:19) '_splice.xlsx']);
                    writetable(thisT, newname, 'FileType', 'spreadsheet', ...
                        'Sheet', runname); 
                end
                
            catch err
                if strcmp(err.identifier, 'MATLAB:spreadsheet:book:openSheetName')
                    warning([runname ' not in sheet'])
                else
                    rethrow(err)
                end
                
            end
                
            
        end
        
    end
    
    fname = newname; 
    
elseif length(var) == 1
    fname = fullfile(var.folder, var.name); 
else
    error('Unknown error when finding files')
end

disp('Done!')

%% Populate onsets structure with preallocated matrixes
% Needs to be updated per experiment. 
for cc = 1:numCons
    thiscon = conNames{cc}; 
    onsets.(thiscon) = nan(events_each, subj.runs(ss)); 
end

%% Grab the correct events
blocks = 1:subj.runs(ss); 
perfect = nan(size(blocks)); 
        
for bb = blocks % run number
    runname = thisrun; 
    if bb > 1
        runname = [runname num2str(bb)]; 
    end

    try
        T = readtable(fname, 'Sheet', runname); 
    catch err
        disp(subj.name)
        rethrow(err)
    end

    stim = T.EventKey;
    resp = T.SubjResponse; 
    key  = T.AnswerKey; 

    noise   = ismember(stim, 193:196); 
    silence = ismember(stim, 197:200); 
    OR = ismember(stim, sort([1:4:192, 2:4:192])); 
    SR = ismember(stim, sort([3:4:192, 4:4:192])); 
    ns = noise | silence; 

    correct = resp == key; 
    correct_OR = correct(OR); 
    correct_SR = correct(SR); 

    accuracy = correct; 
    accuracy(ns) = []; 
    accuracy = mean(accuracy)*100; 
    disp(['Subject answered ' num2str(accuracy) '% sentences correctly!'])

    accuracy_OR = correct_OR; 
    accuracy_OR = mean(accuracy_OR)*100; 
    accuracy_SR = correct_SR; 
    accuracy_SR = mean(accuracy_SR)*100; 
    
    accuracy_OR_SR = accuracy_OR - accuracy_SR; 
    
    correct(noise)   = true; 
    correct(silence) = true; 

    %% Create incorrect regressor
    if all(correct)
        disp(['Perfect accuracy in block ' num2str(bb) '!!!'])
        perfect(bb) = true; 
    else
        perfect(bb) = false; 
        correct = repelem(correct, scan.epis-drop); 

        regname = fullfile(dir_reg, ['incorrect_' scan.runname num2str(bb) '.txt']); 
        fid = fopen(regname, 'w'); 
        incorrect = ~correct; 
        fprintf(fid, '%d\n', incorrect); 
        fclose(fid); 
    end

    %% Create onsets
    idx = 0; 
    idx_cond = zeros(numCons, 1); 
    for ev = 1:height(T)
        idx = idx + 1; 

        if noise(ev)
            idx_cond(1) = idx_cond(1) + 1; % not future proof...
            onsets.NOI(idx_cond(1), bb) = scan.order(idx); 
        elseif silence(ev)
            idx_cond(2) = idx_cond(2) + 1;
            onsets.SIL(idx_cond(2), bb) = scan.order(idx); 
        elseif OR(ev)
            idx_cond(3) = idx_cond(3) + 1;
            onsets.OR(idx_cond(3), bb) = scan.order(idx); 
        elseif SR(ev)
            idx_cond(4) = idx_cond(4) + 1;
            onsets.SR(idx_cond(4), bb) = scan.order(idx); 
        else
            error('Unknown stim type!')
        end
        
        if OR(ev) || SR(ev)
            idx_cond(5) = idx_cond(5) + 1; 
            onsets.LNG(idx_cond(5), bb) = scan.order(idx); 
        end

    end

end

%% Save
fname = fullfile(dir_subj, ['onsets_' scan.runname(1:end-4) '_' design.name '.mat']); 
if exist(fname, 'file')
    delete(fname)
end

save(fname, 'onsets', 'accuracy', 'accuracy_OR', 'accuracy_SR', 'accuracy_OR_SR', 'perfect')

end
