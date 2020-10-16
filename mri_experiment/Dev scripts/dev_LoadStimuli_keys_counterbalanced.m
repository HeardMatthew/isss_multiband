% LoadStimuli.m
% function [audiodata, samplingrate, jitkey, eventkey, anskey, speechkey] = ... 
%     LoadStimuli_keys_counterbalanced(fileloc, events, numSpeech) 
% Loads audio data from stimuli.wav into two cell arrays; one for the
% waveform (audiodata) and one for the sampling rate (samplingrate). Now 
% generates the answer, jitter, and event keys at the same time. 
% Keep all stimuli ALONE in fileloc. Placing any other file into the 
% stimuli folder breaks this script, as it calls directly from dir. Author
% - Matt H

fileloc = 'C:\Users\heard.49\Documents\GitHub\isss_multiband\experiment\stimuli';
events = p.events;
numSpeech = NumberOfSpeechStimuli;

%% Preparing to load stimuli
% Load file names
cd(fileloc)
files = dir; 

% Preallocate speechkey variable
speechkey = zeros(1, numSpeech); 
anskey    = NaN(1, events); 

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

% Check duration of stimuli
info(length(audiodata)) = audioinfo(files(end).name); % Preallocate variable
for i = 3:length(files)
    info(i-2) = audioinfo(files(i).name); 
end

duration = NaN(1, 40); 
for i = 1:40
    duration(i) = info(i).Duration; 
end

durationdiff = NaN(1, 40); 
for i = 2:40
    durationdiff(i) = abs(info(i).Duration - info(i-1).Duration); 
%     if info(i).Duration ~= info(i-1).Duration
%         disp(i)
%         error('Duration of stimuli are not all the same.')
%     end
end

%% Make keys
% jitkey -- How much is the silent period jittered by?
jitkey = 1 + rand(1, events); % Add 1 for silent TR

% speechkey -- Which speech stimuli should we use this run?
% Carefully choose which stimuli to present to ensure stimuli are
% counterbalanced. 
randomstim = Shuffle([0 0 0 0, 1 1 1 1, 2 2 2 2, 3 3 3 3]); 
sentence = 1:4:64; 
speechfile = NaN(1, length(randomstim)); 

for i = 1:length(speechfile)
    speechfile(i) = sentence(i) + randomstim(i); 
end

% eventkey -- In what order will stimuli be presented?
% This is gnarly. In essence, I iterate over the randomized index of
% eventkey to shuffle the events. 
eventkey = Shuffle(horzcat(speechfile, numSpeech+1:length(audiodata)));

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

%% A quick test, for debugging. Run the next few lines to test eventkey. 
% for i = eventkey
%     i
%     sound(audiodata{i}, samplingrate{i})
%     pause(2); 
% end
% 
%     for i = 1:events
%         if eventkey(i) <= 16
%             
%         
%     anskey(i) = mod(eventkey{i}, 2); 