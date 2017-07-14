% LoadStimuli.m
function [audiodata, samplingrate, jitkey, eventkey, anskey, speechkey] = ... 
    LoadStimuli_mono_anskey_v2(fileloc, events, numSpeech) 
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

% Preallocate speechkey variable
speechkey = NaN(1, numSpeech); 
anskey    = NaN(1, events); 

% Preallocate stimuli variables
audiodata = cell(1, events);
samplingrate = cell(1, events);

%% Load ALL audio stimuli
for i = 3:length(files)-1
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

%% Make keys
% jitkey -- How much is the silent period jittered by?
jitkey = 1 + rand(1, events); % Add 1 for silent TR

% speechkey -- Which speech stimuli should we use this run?
% Randomly choose speech stimuli to use, M or F, for each O/S sentence
for i = 1:2:numSpeech
    speechkey(i:i+1) = mod(randperm(2), 2);
end
speechfile = find(speechkey);

% eventkey -- In what order will stimuli be presented?
% This is gnarly. In essence, I iterate over the randomized index of
% eventkey to shuffle the events. 
eventkey = horzcat(speechfile, 33:48);
temp2 = NaN(1, events); 
temp = randperm(events); 
for i = 1:events
    temp2(i) = eventkey(temp(i)); 
end
eventkey = temp2; 

% anskey -- What should have subjects responded with?
% eventkey entries of 33 or higher should have no response. 
for i = 1:events
    if eventkey(i) > 32
        anskey(i) = 0; 
    elseif mod(eventkey(i), 2) == 0 % That is, if actor is male
        anskey(i) = 1; 
    else
        anskey(i) = 2; 
    end
end

% A quick test, for debugging. Run the next few lines to test eventkey. 
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
end