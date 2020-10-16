%% make_physio_regressor(run)
% Does exactly what it says on the tin. Save your zip file holding the
% regressor in the \PHYSIO\run folder, and sit back to watch the magic
% happen. Make sure to update the number of images!!!

% CHANGELOG
% 09/11/17  File inception
% 09/15/17  Updated file to follow current conventions of isss_multi_params
% 02/08/20  Forked for YA_FST
% 02/29/20  Adding input for number of scans as part of subj struct, 
%   updated to use subj input. Revamped for second analysis. 
% 03/20/20  Updated to match new (drop first, etc.) design
% 04/08/20  Cloned for isss_multi

function make_phys_reg(subj, study, ss)
%% Check input
if ~isstruct(subj) || (length(subj) ~= 1)
    error('Input ("subj") where subj is ONE structure')
end

if ~isstruct(study)
    error('Input ("study") where study is structure')
end

if ~isnumeric(ss)
    error('Input ("ss") which specifies which scan!')
end

%% Pathing
dir_subj  = fullfile(study.path, 'data', subj.name); 
dir_physio = fullfile(dir_subj, 'PHYSIO'); 
dir_physio_reg = fullfile(dir_physio, 'reg');
dir_reg = fullfile(dir_subj, 'reg'); 

%% Params
scan = study.scan(ss); 
nscans = (scan.first - 1) + (scan.epis + scan.silence) * (study.events + 1); 
drop = 2/scan.TR; 
nscans_drop = (scan.epis - drop) * study.events; 

%% Read in raw data
niis = dir(fullfile(dir_physio, '*.nii'));
hdrs = dir(fullfile(dir_physio, '*.mat'));

V = cell(1, length(niis));
P = cell(1, length(niis));
for n = 1:length(niis)
    file = fullfile(hdrs(n).folder, hdrs(n).name); 
    load(file)
    P{n} = h.PulseRespiratoryRegressors.ProtocolName; 
    if contains(P{n}, scan.runname)
        file = fullfile(niis(n).folder, niis(n).name); 
        hdr = spm_vol(file); 
        V{n} = spm_read_vols(hdr); 
        V{n} = reshape(V{n}, size(V{n}, 3), size(V{n}, 2));
    end
    
end

skipRuns = []; 
for ii = 1:length(V)
    if length(V{ii}) ~= nscans
        skipRuns = [skipRuns ii];
        disp('Found a terminated run. Skipping...')
    end
    
end

if ~isempty(skipRuns)
    skipRuns = sort(skipRuns, 'descend'); 
    for rr = skipRuns
        P(rr) = [];
        V(rr) = [];
    end
    
end

disp(['Found ' num2str(length(P)) ' files total.'])

lengths = cellfun(@length, V);     
if length(unique(lengths)) ~= 1
    error('Inconsistent number of elements!')
end

for rr = 1:length(P)
    %% Prepare to make regressors
    skipToEnd = 0;
    firstScans = false(scan.first, 1); 
    % Not modeling first dummy event. And remember the sinister TR bug...
    
    template = vertcat(false(scan.silence+drop, 1), true(scan.epis-drop, 1));
    % All models skip the first TR because of artifact
    
    lastScans = false(scan.silence + scan.epis - 1, 1); 
    
    extract = vertcat(firstScans, ...
        repmat(template, [study.events, 1]), ...
        lastScans); 

    try
        reg = []; 
        for zz = 1:subj.runs
            temp = V{rr}; 
            
            for cc = 1:size(temp, 2)
                reg = [reg, temp(extract, cc)];                
            end
            
        end
        
        assert(size(reg, 1) == nscans_drop)
    catch
        disp(['Something happened during ' P{rr} ' and the regressor was not created.'])
        skipToEnd = 1;
    end

    if ~skipToEnd
        disp(['Saving ' P{rr}])

         if scan.regnum < 10
            tag = ['0000' num2str(scan.regnum)]; 
        elseif scan.regnum < 100
            tag = ['000' num2str(scan.regnum)]; 
        else
            error('Unexpected regnum!')
        end
        
        filename1 = fullfile(dir_reg, ['physio_full_' P{rr} '_' tag '.txt']);
        fid = fopen(filename1, 'w');
        for line = 1:length(reg)
            eline = '\t';
            for ii = 1:8
                if reg(line, ii) < 0
                    eline = [eline, '%e   '];
                elseif reg(line, ii) >= 0
                    eline = [eline, ' %e   '];
                end
            end
            fprintf(fid, eline, reg(line, :));
            fprintf(fid, '\n');
        end
        fclose(fid); 
            
        filename2 = fullfile(dir_reg, ['physio_1st_' P{rr} '_' tag '.txt']);
        fid = fopen(filename2, 'w');
        for line = 1:length(reg)
            eline = '\t';
            for ii = [1 2 5 6]
                if reg(line, ii) < 0
                    eline = [eline, '%e   '];
                elseif reg(line, ii) >= 0
                    eline = [eline, ' %e   '];
                end
            end
            fprintf(fid, eline, reg(line, [1 2 5 6]));
            fprintf(fid, '\n');
        end
        fclose(fid); 
    end

end

end
