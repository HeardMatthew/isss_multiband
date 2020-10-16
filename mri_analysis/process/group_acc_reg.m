%% group_acc_reg
% Creates a group-level regressor that tracks task accuracy.
% 
% MM/DD/YY -- CHANGELOG
% 04/08/20 -- File initialized. 


function group_acc_reg(subj, study, ss)
%% Check input
if ~isstruct(subj) || length(subj) < 2
    error('subj has all subjects!')
end

% whatever, I'm getting lazy. 

%% Pathing
cd ..
dir_batch = pwd; 
dir_mlb   = fullfile(dir_batch, 'matlabbatch'); 
dir_data  = fullfile(study.path, 'data'); 

%% Parameters
scan = study.scan(ss); 
scanname = scan.runname(1:end-4); 

%% Preallocate regressor
group_level_acc = zeros(length(subj), 1); 
group_level_OR_SR = zeros(length(subj), 1); 

%% Load onsets
for ss = 1:length(subj)
    thissubj = subj(ss); 
    dir_subj  = fullfile(study.path, 'data', thissubj.name); 
    file = fullfile(dir_subj, ['onsets_' scanname '.mat']);
    load(file)
    
    %% Poke in accuracy
    group_level_acc(ss) = accuracy; 
    group_level_OR(ss)  = accuracy_OR; 
    group_level_SR(ss)  = accuracy_SR; 
    group_level_OR_SR(ss) = accuracy_OR_SR;
end

%% Save the regressor
filename = fullfile(dir_batch, ['group_level_acc_' scanname '.mat']); 
save(filename, 'group_level_acc', 'group_level_OR', 'group_level_SR', 'group_level_OR_SR')

end
