%% SPA_manipulate
% Combines data across runs and conditions to form AUE contrasts
% 
% MM/DD/YY -- CHANGELOG
% 04/07/20 -- Initialized changelog. Cloned from original hybrid_isss batch
%   for new analysis. Updated format to match new universal strategy. 

function SPA_manipulate(varargin)
%% Check input
subj = varargin{1}; 
study = varargin{2}; 
dd = varargin{3}; 
ss = varargin{4}; 

if ~isstruct(subj) || length(subj) ~= 1
    error('Input ("subj") where subj is a SINGLE struct')
end

if ~isstruct(study)
    error('Input ("study") which has experiment info')
end

if ~isnumeric(dd)
    error('Input ("dd") which specifies which designs (#) should be used')
end

if ~isnumeric(ss)
    error('Input ("ss") which specifies scan!')
end

%% Pathing
cd ..
dir_batch = pwd; 
dir_subj  = fullfile(study.path, 'data', subj.name); 
dir_design_root = fullfile(dir_subj, 'design');

%% Parameters
scan = study.scan(ss); 
scanname = scan.runname(1:end-4); 
thisdesign = study.design(dd); 
dir_thisdesign = fullfile(dir_design_root, [scanname '_' thisdesign.name]); 

numruns = subj.runs(ss);     
% numcond = length(thisdesign.cond); 
numcont = length(thisdesign.con.name); 
numepis = 8/scan.TR; 

%% Load data
for ii = 1:numcont
    these_conditions = thisdesign.con.vec(ii, :) ~= 0;  
    this_math = thisdesign.con.vec(ii, these_conditions); 

    conditions = thisdesign.cond(these_conditions);
    num_conditions = length(conditions); 

    V = cell(num_conditions, numruns);
    filenames = cell(num_conditions, numruns);
    
    for rr = 1:numruns
        for tt = 1:num_conditions
            filenames{tt, rr} = fullfile(dir_thisdesign, [ ...
                'AUE_' conditions{tt} '_run' num2str(rr) '.nii' ...
                ]); 
            V{tt, rr} = spm_vol(filenames{tt, rr}); 
            if rr == 1 && tt == 1
                math = zeros([V{tt, rr}.dim, num_conditions, numruns]);
            end
            
            Img = spm_read_vols(V{tt, rr}); 
            
            math(:, :, :, tt, rr) = this_math(tt) * Img; 
        end
        
    end
    
    %% Combine across runs and do the math
    % Math has size(X, Y, Z, condition, run)
    math_across_runs = mean(math, 5); 
    combined_AUE = sum(math_across_runs, 4); 
        
    %% Save image
    V_merge = V{1, 1}; 
    V_merge.fname = fullfile(dir_thisdesign, ['AUE_' thisdesign.con.name{ii} '.nii']); 
    spm_write_vol(V_merge, combined_AUE);
    
end

end


