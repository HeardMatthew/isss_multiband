% LoadStimuli.m
function [audiodata, samplingrate, duration, jitkey, eventkey, anskey, rhythmkey] = ... 
    LoadStimuliAndKeys_rhythm(fileloc, events, run) 
% Loads audio data from \stimuli into two cell arrays; one for the
% waveform (audiodata) and one for the sampling rate (samplingrate). Now 
% generates the answer, jitter, and event keys at the same time. 
% Author - Matt H

% Set variables while debugging
%  fileloc = StimuliLoc; 
%  events = p.events; 
%  run = p.runNum; 
%  numRhythm = NumberOfRTStimuli; 

%% Preparing to load stimuli
% Load file names
run = str2double(run); % Input is a string

cd(fileloc)
files = dir('*.wav'); 

% Preallocate anskey variable
anskey = NaN(1, events); 

% Preallocate stimuli variables
audiodata = cell(1, length(files));
samplingrate = cell(1, length(files));

%% Load ALL audio stimuli
for i = 1:length(files)
    [tempAudio, tempFs] = audioread(files(i).name);
    audiodata{i} = [tempAudio'; tempAudio'];
    samplingrate{i} = tempFs;
end

% Check samplingrate is same across files (assumption made in main code)
for i = 2:length(samplingrate)
    if samplingrate{i} ~= samplingrate{i-1}
        error('Your sampling rates are not all the same. Stimuli will not play correctly.')
        % A quick check to make sure stimuli will all play correctly
    end
end

% Load duration of files
clear info
info(length(audiodata)) = audioinfo(files(end).name); % Preallocate struct
for i = 1:length(files)
    info(i) = audioinfo(files(i).name); 
end

duration = NaN(1, length(info)); 
for i = 1:length(duration)
    duration(i) = info(i).Duration; 
end

%% Make keys

% eventkey -- In what order will stimuli be presented?
% Carefully choose which stimuli to present to ensure stimuli are
% counterbalanced. 
eventkey = Shuffle(horzcat( ... 
    1 * ones(1, 4), ... % Simple long
    2 * ones(1, 4), ... % Complex long
    3 * ones(1, 4), ... % Simple short
    4 * ones(1, 4), ... % Complex short
    5, 6, 7, 8 ...      % Oddball, oddball, silence, silence
    )); 

% jitkey -- How much was the jitter before stimulus presentation? This
% depends on which stimuli is being presented, as durations vary across one
% condition. 
shortstimduration = 2.100; 
longstimduration  = 3.600; 
jitkey = NaN(1, events); 
for i = 1:events
    if ~isempty(find(eventkey(i) == [1 2 5], 1)) % Long conditions
        jitkey(i) = (4-longstimduration)*rand(1); 
    elseif ~isempty(find(eventkey(i) == [3 4 6], 1)) % Short conditions
        jitkey(i) = (4-shortstimduration)*rand(1); 
    elseif ~isempty(find(eventkey(i) == [7 8], 1)) % Silent conditions
        jitkey(i) = 2*rand(1); 
    end
end

% rhythmkey -- Which rhythms were presented to the subject? 
rhythmkey = sort(eventkey);

% anskey -- What should have subjects responded with?
for i = 1:events
    if ~isempty(find(eventkey(i) == [5 6], 1)) % events that are odd
        anskey(i) = 1; 
    else % events that are not odd
        anskey(i) = 0; 
    end
end

end