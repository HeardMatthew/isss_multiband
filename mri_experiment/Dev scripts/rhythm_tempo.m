%% rhythm_tempo.m
% Script to run for pilot scans testing rhythm and tempo discrimination. 
% Ported from isss_multiband script. 
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
    'First run (7-8)', ... 
    'Last run (7-8)', ... 
    'Scanner connected (0/1):', ...
    'RTBox connected (0/1):', ...
    }; 
dlg_ans = inputdlg(prompt); 

% Convert dlg_ans into my p.arameters struct
subj.Num  = dlg_ans{1};
subj.Init = dlg_ans{2}; 
subj.firstRun = dlg_ans{3}; 
subj.lastRun  = dlg_ans{4}; 
ConnectedToScanner = str2double(dlg_ans{5});
ConnectedToRTBox   = str2double(dlg_ans{6}); 

NumberOfRTStimuli = 4;  % 4 different rhythm/tempo clips. Needs to be 4. 
NumberOfStimuli   = 8;  % 8 .wav files in folder (norm + oddball + silent)
% The script presents four sets of the 4 RT stimuli per run in a random 
% order, and includes two oddball and two silent stimuli. 

scan.TR = 1.000; 
scan.epiNum = 10; 

% Timing
p.runs        = subj.lastRun - subj.firstRun + 1; 
p.events      = 20;     % 20 events (16 stim + 2 silent + 2 oddball)
p.presTime    = 4.000;  % 4 seconds
p.epiTime     = 10.000; % 10 seconds
p.eventTime   = p.presTime + p.epiTime;
p.runDuration = p.epiTime + ...   % After first pulse
    p.eventTime * p.events + ...  % Each event
    p.eventTime;                  % After last acquisition
p.rxnWindow   = 3.000;  % 3 seconds
p.jitWindow   = [1.000, 2.300];  % Varies, see notes below
    % For this experiment, each stimuli onset will be jittered. However,
    % stimuli may vary in duration. "Long" stimuli last 3.0 seconds, while
    % "short" stimuli are 1.7 seconds. So, jitter can be either up to 1.000
    % or 2.300 seconds long. This is taken care of in LoadStimuli. 

%% Paths
cd ..
direc = pwd; 
StimuliTag = [direc, '\stimuli_rhythm']; 
if mod(str2double(p.subjNum), 2) == 1 % Odd numbered subject
    StimuliLoc   = [StimuliTag, '\oddsubj'];
else                                  % Even numbered subject
    StimuliLoc   = [StimuliTag, '\evensubj'];
end

ScriptsLoc   = [direc, '\scripts'];
ResultsLoc   = [direc, '\results']; 
FuncsLoc     = [ScriptsLoc, '\functions'];
Instructions = 'instructions.txt';

cd ..
RTBoxLoc = [pwd, '\RTBox']; 

filetag    = [p.subjNum '_' p.subjInit '_' p.runNum '_rhythm']; 
ResultsTxt = [filetag '_results.txt']; 
ResultsXls = [filetag '_results.xlsx']; 
Variables  = [filetag '_variables.mat']; 

%% Preallocating timing variables (hard-coded with n = 2)
n = 2; 
eventStart   = NaN(p.events, n);
stimStart    = NaN(p.events, n); 
stimEnd      = NaN(p.events, n); 
eventEnd     = NaN(p.events, n); 
stimEndKey   = NaN(p.events, n); 
stimDuration = NaN(p.events, n); 

eventStartKey = NaN(p.events, n); 

respTime = cell(1, p.events); 
respKey  = cell(1, p.events); 
    % I use cells here so that responses can be empty when subjects time 
    % out, and because responses from RTBox come back as strings. 
    
firstPulse = NaN(1, n);
runEnd     = NaN(1, n); 

%% Prepare test
% Load stimuli and check counterbalance
cd(FuncsLoc) 
[audio, fs, stimDuration, jitterKey, eventKey, answerKey, rhythmKey] = ...
    LoadStimuliAndKeys_rhythm(StimuliLoc, p.events, p.runNum);
fs = fs{1}; % Above func checks that all fs are the same. 

% if ~strcmp(p.scanType, 'train')
%     cd(FuncsLoc)
%     stimulicheck(NumberOfRTStimuli, eventKey); 
% end
cd(direc)

% Prepare timing keys
%eventStartKey = p.epiTime + [0:p.eventTime:((p.events-1)*p.eventTime)]; %#ok<NBRAK>
%stimStartKey  = eventStartKey + jitterKey; 
%stimEndKey    = stimStartKey + stimDuration(eventKey);
%rxnEndKey     = stimEndKey + p.rxnWindow; 
%eventEndKey   = eventStartKey + p.eventTime;

% Open PTB screen on scanner, prepare fixation cross coords
[wPtr, rect] = Screen('OpenWindow', 0, 185);
frameDuration = Screen('GetFlipInterval', wPtr);
centerX = rect(3)/2;
centerY = rect(4)/2;
crossCoords = [-30, 30, 0, 0; 0, 0, -30, 30]; 
HideCursor(); 

% Open audio connection
pahandle = PsychPortAudio('Open', [], [], [], fs);
% Play a silent audio clip just to eliminate clipping on first presentation
PsychPortAudio('FillBuffer', pahandle, audio{rhythmKey(end)});
PsychPortAudio('Start', pahandle, 1);
    % Stimuli are presented on the scanner computer, which shares a screen
    % and audio output with the scanner projector and headphones. 

% Check if using RTBox or Keyboard
if ConnectedToRTBox == 0
    cd(RTBoxLoc)
    RTBox('fake', 1)
    cd(direc)
end

% Display instructions
% if ShowInstructions
%     cd(FuncsLoc)
%     DisplayInstructions_bkfw_rtbox(Instructions, wPtr, RTBoxLoc); 
%     cd(direc)
% end
DrawFormattedText(wPtr, 'Waiting for experimenters...'); 
Screen('Flip', wPtr); 

%% Prepare test

for run = subj.firstRun:subj.lastRun

    DrawFormattedText(wPtr, 'Please wait, preparing run...');
    Screen('Flip', wPtr); 
    
    filetag = [subj.Num '_' subj.Init]; 
    
    ResultsTxt  = [filetag 'rhythm_results.txt']; 
    ResultsXls  = [filetag 'rhythm_results.xlsx']; 
    Variables   = [filetag 'rhythm_variables.mat']; 

    
    
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
        PsychPortAudio('Start', pahandle);
        WaitTill(AbsStimEnd(j)); 
        PsychPortAudio('Stop', pahandle); 
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
OutputData_rhythm
cd(ScriptsLoc)

%% Closing down
Screen('CloseAll');
PsychPortAudio('Close'); 
cd(ScriptsLoc)
DisableKeysForKbCheck([]); 