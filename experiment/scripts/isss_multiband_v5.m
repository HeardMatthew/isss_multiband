%% isss_multiband.m
% Script to run for pilot scans comparing ISSS, multiband, and the new
% hybrid scanning protocol. Ported from my previous ISSS_test script. 
% Author - Matt Heard
%% Startup
sca; DisableKeysForKbCheck([]); KbQueueStop; 
Screen('Preference','VisualDebugLevel', 0); 

try 
    PsychPortAudio('Close'); 
catch
    disp('PsychPortAudio is already closed.')
end

InitializePsychSound

clearvars; clc; 
codeStart = GetSecs(); 

%% Parameters
% Many of these parameters need to be determined while testing. 
AudioDevice = PsychPortAudio('GetDevices', 3); 
    
prompt = {...
    'Subject number (###):', ...
    'Subject initials (XX):', ...
    'Scan protocol (train/isss/multi/hybrid):', ...
    'Run number (0-6)', ...
    'Scanner connected (0/1):', ...
    'RTBox connected (0/1):', ...
    }; 
dlg_ans = inputdlg(prompt); 

% Convert dlg_ans into my p.arameters struct
p.subjNum  = dlg_ans{1};
p.subjInit = dlg_ans{2}; 
p.scanType = dlg_ans{3}; 
p.runNum   = dlg_ans{4}; 
ConnectedToScanner = str2double(dlg_ans{5});
ConnectedToRTBox   = str2double(dlg_ans{6}); 

ShowInstructions = 0; 

p.TR     = 1.000; 
p.epiNum = 10; 
if strcmp(p.scanType, 'train')
    if ~strcmp(p.runNum, '0')
        error('Set Run number to 0')
    end
    ShowInstructions = 1;
elseif strcmp(p.scanType, 'isss')
    p.TR     = 2.000;
    p.epiNum = 5;
elseif isempty(find(strcmp(p.scanType, {'hybrid', 'multi'}), 1))
    error('Invalid scan protocol')
end

NumberOfSpeechStimuli = 192; % 192 different speech clips
NumberOfStimuli       = 200; % 200 .wav files in stimuli folder
% Of these 200 .wav files, 4 are silent, 4 of noise, and 192 are sentences.
% Of these 192 speech sounds, this script chooses 8 per run (2 MO, 2 FO, 2 
% MS, 2 FS) for presentation. Subjects will not hear a
% repeated "sentence structure" (i.e. one stimuli of of 001, one stimuli of
% 002) in the entire experiment. 

p.events      = 16; 
if strcmp(p.scanType, 'train')
    p.events  = 10; 
    NumberOfSpeechStimuli = 8; 
    NumberOfStimuli       = 10; 
end
p.presTime    = 4.000;  % 4 seconds
p.epiTime     = 10.000; % 10 seconds
p.eventTime   = p.presTime + p.epiTime;
p.runDuration = p.epiTime + ...   % After first pulse
    p.eventTime * p.events + ...  % Each event
    p.eventTime;                  % After last acquisition
p.rxnWindow = 3.000;  % 3 seconds
p.jitWindow = 1.000;  % 1 second, see notes below
    % For this experiment, the first second of the silent window will not
    % have stimuli presented. To code for this, I add an additional 1 s
    % to the jitterKey. So, the jitter window ranges from 1 s to 2 s.

% Paths
cd ..
direc = pwd; 

StimuliLoc   = [direc, '\stimuli_lang'];
ScriptsLoc   = [direc, '\scripts'];
ResultsLoc   = [direc, '\results']; 
FuncsLoc     = [ScriptsLoc, '\functions'];
Instructions = 'instructions.txt';

cd ..
RTBoxLoc = [pwd, '\RTBox']; 

filetag    = [p.subjNum '_' p.subjInit '_' p.runNum '_' p.scanType]; 
ResultsTxt = [filetag '_results.txt']; 
ResultsXls = [filetag '_results.xlsx']; 
Variables  = [filetag '_variables.mat']; 

% Preallocating timing variables
eventStart = NaN(1, p.events);
stimStart  = NaN(1, p.events); 
stimEnd    = NaN(1, p.events); 
eventEnd   = NaN(1, p.events); 

respTime = cell(1, p.events); 
respKey  = cell(1, p.events); 
    % I use cells here so that responses can be empty when subjects time 
    % out, and because responses from RTBox come back as strings. 

%% Prepare test
% Load stimuli and check counterbalance
cd(FuncsLoc) 
[audio, fs, stimDuration, jitterKey, eventKey, answerKey, speechKey] = ...
    LoadStimuliAndKeys(StimuliLoc, p.events, p.runNum, NumberOfSpeechStimuli);
fs = fs{1}; % Above func checks that all fs are the same. 

if ~strcmp(p.scanType, 'train')
    cd(FuncsLoc)
    stimulicheck(NumberOfSpeechStimuli, eventKey); 
end
cd(direc)

% Prepare timing keys
eventStartKey = p.epiTime + [0:p.eventTime:((p.events-1)*p.eventTime)]; %#ok<NBRAK>
stimStartKey  = eventStartKey + jitterKey; 
stimEndKey    = stimStartKey + stimDuration(eventKey);
rxnEndKey     = stimEndKey + p.rxnWindow; 
eventEndKey   = eventStartKey + p.eventTime;

% Open PTB screen on scanner, prepare fixation cross coords
[wPtr, rect] = Screen('OpenWindow', 0, 185);
frameDuration = Screen('GetFlipInterval', wPtr);
centerX = rect(3)/2;
centerY = rect(4)/2;
crossCoords = [-30, 30, 0, 0; 0, 0, -30, 30]; 
HideCursor(); 

% Open audio connection
InitializePsychSound
pahandle = PsychPortAudio('Open', [], [], [], fs);
    % Stimuli are presented on the scanner computer, which shares a screen
    % and audio output with the scanner projector and headphones. 

% Check if using RTBox or Keyboard
if ConnectedToRTBox == 0
    cd(RTBoxLoc)
    RTBox('fake', 1)
    cd(direc)
end

% Display instructions
if ShowInstructions
    cd(FuncsLoc)
    DisplayInstructions_bkfw_rtbox(Instructions, wPtr, RTBoxLoc); 
    cd(direc)
end
DrawFormattedText(wPtr, 'Waiting for experimenters...'); 
Screen('Flip', wPtr); 

% Wait for first pulse
cd(RTBoxLoc)
RTBox('Clear'); 
RTBox('UntilTimeout', 1);
firstPulse = RTBox('WaitTR'); 

% Draw onto screen after recieving first pulse
Screen('DrawLines', wPtr, crossCoords, 2, 0, [centerX, centerY]);
Screen('Flip', wPtr); 

% Generate absolute time keys
AbsEvStart   = firstPulse + eventStartKey; 
AbsStimStart = firstPulse + stimStartKey; 
AbsStimEnd   = firstPulse + stimEndKey; 
AbsRxnEnd    = firstPulse + rxnEndKey; 
AbsEvEnd     = firstPulse + eventEndKey; 

WaitTill(firstPulse + p.epiTime); 

%% Present audio stimuli
try
    for j = 1:p.events
        eventStart(j) = GetSecs(); 
        PsychPortAudio('FillBuffer', pahandle, audio{eventKey(j)});

        WaitTill(AbsStimStart(j)); 
        
        stimStart(j) = GetSecs; 
        PsychPortAudio('Start', pahandle, 1);
        WaitTill(AbsStimEnd(j)); 
        stimEnd(j) = GetSecs; 
        RTBox('Clear'); 

        [respTime{j}, respKey{j}] = RTBox(AbsRxnEnd(j)); 

        [~, eventEnd(j)] = WaitTill(AbsEvEnd(j));    
    end
    WaitSecs(p.eventTime); 
    runEnd = GetSecs(); 
    sca;  
catch err
    sca; 
    rethrow(err)
end

%% Save data
cd(FuncsLoc)
OutputData
cd(ScriptsLoc)

%% Closing down
Screen('CloseAll');
PsychPortAudio('Close'); 
cd(ScriptsLoc)
DisableKeysForKbCheck([]); 