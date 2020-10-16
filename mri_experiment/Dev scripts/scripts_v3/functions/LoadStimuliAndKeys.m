% LoadStimuli.m
function [audiodata, samplingrate, duration, jitkey, eventkey, anskey, speechkey] = ... 
    LoadStimuliAndKeys(fileloc, events, numSpeech) 
% Loads audio data from stimuli.wav into two cell arrays; one for the
% waveform (audiodata) and one for the sampling rate (samplingrate). Now 
% generates the answer, jitter, and event keys at the same time. 
% Keep all stimuli ALONE in fileloc. Placing any other file into the 
% stimuli folder breaks this script, as it calls directly from dir. Author
% - Matt H

%% Preparing to load stimuli
% Load file names
cd(fileloc)
files = dir; 

% Preallocate anskey variable
anskey = NaN(1, events); 

% Preallocate stimuli variables
audiodata = cell(1, events);
samplingrate = cell(1, events);

%% Load ALL audio stimuli
for i = 3:length(files)
    [tempAudio, tempFs] = audioread(files(i).name);
    audiodata{i-2} = [tempAudio'; tempAudio'];
    samplingrate{i-2} = tempFs;
end
clear tempAudio
clear tempFs

% Check samplingrate is same across files (assumption made in main code)
for i = 2:length(samplingrate)
    if samplingrate{i} ~= samplingrate{i-1}
        error('Your sampling rates are not all the same. Stimuli will not play correctly.')
        % A quick check to make sure stimuli will all play correctly
    end
end

% Load duration of files
info(length(audiodata)) = audioinfo(files(end).name); % Preallocate variable
for i = 3:length(files)
    info(i-2) = audioinfo(files(i).name); 
end

duration = NaN(1, length(info)); 
for i = 1:length(duration)
    duration(i) = info(i).Duration; 
end

%% Make keys
% jitkey -- How much is the silent period jittered by?
jitkey = 1 + rand(1, events); % Add 1 for silent TR

% speechkey -- Which speech stimuli should we use this run?
% Carefully choose which stimuli to present to ensure stimuli are
% counterbalanced. 
randomstim = Shuffle(horzcat( ... 
    0 * ones(1, 4), ... 
    1 * ones(1, 4), ... 
    2 * ones(1, 4), ... 
    3 * ones(1, 4) ... 
    )); 
sentence = 1:4:64; 
speechkey = NaN(1, length(randomstim)); 

for i = 1:length(speechkey)
    speechkey(i) = sentence(i) + randomstim(i); 
end

% eventkey -- In what order will stimuli be presented?
% This is gnarly. In essence, I iterate over the randomized index of
% eventkey to shuffle the events. 
eventkey = Shuffle(horzcat(speechkey, numSpeech+1:length(audiodata)));

% anskey -- What should have subjects responded with?
for i = 1:events
    if eventkey(i) > numSpeech % events that are silence or noise
        anskey(i) = 0; 
    elseif mod(eventkey(i), 2) == 0 % That is, if actor is male
        anskey(i) = 1; 
    else % That is, if actor is female
        anskey(i) = 2; 
    end
end

end