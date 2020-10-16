% LoadStimuli.m
function [audiodata, samplingrate, jitkey, eventkey, anskey] = LoadStimuli_mono_anskey(fileloc, runs, eventsPerRun) 
% Loads audio data from stimuli.wav into two cell arrays; one for the
% waveform (audiodata) and one for the sampling rate (samplingrate). Now 
% generates the answer, jitter, and event keys at the same time. 
% Keep all stimuli ALONE in fileloc. Placing any other file into the 
% stimuli folder breaks this script, as it calls directly from dir. Author
% - Matt H

% Load file names
cd(fileloc)
files = dir; 

% Preallocate key variables
jitkey = cell(1, runs);
eventkey = cell(1, runs);
anskey = cell(1, runs);
for i = 1:runs
    anskey{i} = NaN(1, eventsPerRun); 
end

% Preallocate stimuli variables
audiodata = cell(1, length(files)-3);
samplingrate = cell(1, length(files)-3);
% Note the -4 above, due to '.' '..' 'clear.zip' '\unused'

% Load audio stimuli
for i = 3:length(files)-1
    [tempAudio, tempFs] = audioread(files(i).name);
    audiodata{i-2} = [tempAudio'; tempAudio'];
    samplingrate{i-2} = tempFs;
end
clear tempAudio
clear tempFs

for i = 4:length(files)-4
    if samplingrate{i} ~= samplingrate{i-1}
        error('Your sampling rates are not all the same. Stimuli will not play correctly.')
        % A quick check to make sure stimuli will all play correctly
    end
end

% Make keys
for i = 1:runs
    jitkey{i} = 1 + rand(1, eventsPerRun); % Add 1 for silent TR
    eventkey{i} = randperm(eventsPerRun); 
    anskey{i} = mod(eventkey{i}, 2); 
end