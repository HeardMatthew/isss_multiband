% LoadStimuli.m
function [audiodata, samplingrate, duration, jitkey, eventkey, anskey, speechkey] = ... 
    LoadStimuliAndKeys(fileloc, events, run, numSpeech) 
% Loads audio data from stimuli.wav into two cell arrays; one for the
% waveform (audiodata) and one for the sampling rate (samplingrate). Now 
% generates the answer, jitter, and event keys at the same time. 
% Keep all stimuli in fileloc or training.Author - Matt H

% CHANGELOG
% 7/12/17 Started keeping record of changelog on file. -- MH

% Set variables while debugging
%  fileloc = StimuliLoc; 
%  events = p.events; 
%  run = p.runNum; 
%  numSpeech = NumberOfSpeechStimuli; 

%% Preparing to load stimuli
% Load file names
run = str2double(run); % Input is a string

if run == 0
    cd([fileloc filesep 'training'])
else
    cd(fileloc)
end
files = dir('*.wav'); 

% Preallocate anskey variable
anskey = NaN(1, events); 

% Preallocate stimuli variables
audiodata = cell(1, length(files));
samplingrate = cell(1, length(files));

%% Load ALL audio stimuli
for i = 1:length(files)
    [tempAudio, tempFs] = audioread(files(i).name);
    audiodata{i} = [tempAudio'; tempAudio']; % Convert mono to stereo
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
% jitkey -- How much is the silent period jittered by?
jitkey = 1 + rand(1, events); % Add 1 because stimuli are short-ish

% speechkey -- Which speech stimuli should we use this run?
% Carefully choose which stimuli to present to ensure stimuli are
% counterbalanced. 
randomstim = Shuffle(horzcat( ... % Half of events are speech stim
    0 * ones(1, events/8), ...    % Four conditions of speech means...
    1 * ones(1, events/8), ...    % events/(2*4) = events/8
    2 * ones(1, events/8), ... 
    3 * ones(1, events/8) ... 
    )); 
sentence = ((run-1)*32)+1:4:run*events*2; % just figure this out in adv.
speechkey = sentence + randomstim; 

% eventkey -- In what order will stimuli be presented?
eventkey = Shuffle(horzcat(speechkey, numSpeech+1:length(audiodata)));

if run == 0
    speechkey = 1:8;
    eventkey  = Shuffle(1:10); 
    jitkey = 1 + rand(1, 10); 
end

% anskey -- What should have subjects responded with?
for i = 1:events
    if eventkey(i) > numSpeech % events that are silence or noise
        anskey(i) = 0; 
    elseif mod(eventkey(i), 2) == 0 % That is, if actor is male
        anskey(i) = 2; 
    else % That is, if actor is female
        anskey(i) = 1; 
    end
end

end