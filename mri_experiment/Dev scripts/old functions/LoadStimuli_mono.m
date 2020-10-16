% LoadStimuli.m
function [audiodata, samplingrate] = LoadStimuli_mono(fileloc) 
% Loads audio data from stimuli.wav into two cell arrays; one for the
% waveform (audiodata) and one for the sampling rate (samplingrate). Keep
% all stimuli ALONE in fileloc. Placing any other file into the stimuli
% folder breaks this script, as it calls elements directly from dir. Author
% - Matt H
cd(fileloc)
files = dir; 
audiodata = cell(1, length(files)-4);   % Preallocate audio data cell
samplingrate = cell(1, length(files)-4);% Note the -4 above, due to 
for i = 3:length(files)-2               % '.' '..' 'clear.zip' '\unused'
    [tempAudio, tempFs] = audioread(files(i).name);
    audiodata{i-2} = [tempAudio'; tempAudio'];
    samplingrate{i-2} = tempFs;
end
clear tempAudio
clear tempFs
