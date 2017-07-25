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
    'First run (0-6, input 0 for training block)', ... 
    'Last run (0-6, input 0 for training block)', ... 
    'Scanner connected (0/1):', ...
    'RTBox connected (0/1):', ...
    }; 
dlg_ans = inputdlg(prompt); 

% Convert dlg_ans into subj struct
subj.Num  = dlg_ans{1};
subj.Init = dlg_ans{2}; 
subj.firstRun = str2double(dlg_ans{3}); 
subj.lastRun  = str2double(dlg_ans{4}); 
ConnectedToScanner = str2double(dlg_ans{5});
ConnectedToRTBox   = str2double(dlg_ans{6}); 

% Training exception
if subj.firstRun == 0
    Training = 1; 
    subj.firstRun = 1; 
    subj.lastRun = 1; 
else
    Training = 0; 
end

% Structures per scan type
scanType.hybrid.TR     = 1.000; 
scanType.hybrid.epiNum = 10; 
scanType.multi.TR      = 1.000; 
scanType.multi.epiNum  = 14; 
scanType.isss.TR       = 2.000; 
scanType.isss.epiNum   = 5; 

NumSpeechStimuli = 192; % 192 different speech clips
NumStimuli       = 200; % 200 .wav files in stimuli folder
% Of these 200 .wav files, 4 are silent, 4 of noise, and 192 are sentences.
% Of these 192 speech sounds, this script chooses 8 per run (2 MO, 2 FO, 2 
% MS, 2 FS) for presentation. Subjects will not hear a
% repeated "sentence structure" (i.e. one stimuli of of 001, one stimuli of
% 002) in the entire experiment.

% Timing
p.runs        = subj.lastRun - subj.firstRun + 1; 
p.events      = 16; 
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
    
%% Paths
cd ..
direc = pwd; 
StimuliLoc   = [direc, '\stimuli_lang'];
ScriptsLoc   = [direc, '\scripts'];
ResultsLoc   = [direc, '\results'];
FuncsLoc     = [ScriptsLoc, '\functions'];
Instructions = 'instructions_lang.txt';
cd ..
RTBoxLoc = [pwd, '\RTBox']; 

% Training exceptions    
if Training
    NumSpeechStimuli = 4; 
    NumStimuli       = 5; 
    p.events = 5; 
end

%% Preallocating timing variables (hard-coded width so if runs = (not 1):n, not broken...
if Training
    n = 1;
else
    n = 6; 
end

eventStart   = NaN(p.events, n);
stimStart    = NaN(p.events, n); 
stimEnd      = NaN(p.events, n); 
eventEnd     = NaN(p.events, n); 
stimEndKey   = NaN(p.events, n); 
stimDuration = NaN(p.events, n); 

eventStartKey = NaN(p.events, n); 

respTime = cell(p.events, n); 
respKey  = cell(p.events, n); 
    % I use cells here so that responses can be empty when subjects time 
    % out, and because responses from RTBox come back as strings. 

firstPulse = NaN(1, n); 
runEnd     = NaN(1, n); 
    
%% Load protocol master list, PTB, and stimuli
cd(ScriptsLoc)
load('master_scan_order.mat')
protocol_order = cell(1, 6); 
for event = 1:10
    if master{1, event} == subj.Num
        for k = 1:6
            protocol_order{k} = master{k+1, event}; 
        end
    end
end
cd(direc)

% Open PTB screen on scanner, prepare fixation cross coords
[wPtr, rect] = Screen('OpenWindow', 0, 185);
DrawFormattedText(wPtr, 'Please wait, preparing experiment...');
Screen('Flip', wPtr); 
centerX = rect(3)/2;
centerY = rect(4)/2;
crossCoords = [-30, 30, 0, 0; 0, 0, -30, 30]; 
HideCursor(); 

% Load stimuli and check counterbalance
cd(FuncsLoc) 
[audio, fs, rawStimDuration, jitterKey, eventKey, answerKey, speechKey] = ...
    LoadStimuliAndKeys(StimuliLoc, p.events, subj.firstRun, subj.lastRun, NumSpeechStimuli, Training);
fs = fs{1}; % Above func checks that all fs are the same.  

if Training
    protocol_order = {'multi'}; 
else
    cd(FuncsLoc)
    stimulicheck(NumSpeechStimuli, eventKey); 
end
cd(direc)

% Open audio connection
pahandle = PsychPortAudio('Open', [], [], [], fs);
    % Stimuli are presented on the scanner computer, which shares a screen
    % and audio output with the scanner projector and headphones. 
    
% Check if using RTBox or Keyboard
if ConnectedToRTBox == 0
    cd(RTBoxLoc)
    RTBox('fake', 1); 
    cd(direc)
end
    
%% Prepare test

for run = subj.firstRun:subj.lastRun
    
    DrawFormattedText(wPtr, 'Please wait, preparing run...');
    Screen('Flip', wPtr); 
    
    % File names
    if Training
        filetag = [subj.Num '_' subj.Init '_practice']; 
    else
        filetag = [subj.Num '_' subj.Init]; 
    end
    ResultsTxt  = [filetag '_results.txt']; 
    ResultsXls  = [filetag '_results.xlsx']; 
    Variables   = [filetag '_variables.mat']; 

    % Prepare timing keys
    eventStartKey(:,run) = p.epiTime + [0:p.eventTime:((p.events-1)*p.eventTime)]'; %#ok<NBRAK>
    stimStartKey = eventStartKey + jitterKey; 
    
    if Training
        stimEndKey = stimStartKey + rawStimDuration(eventKey)';
    else
        stimEndKey(:, subj.firstRun:subj.lastRun) = ... 
            stimStartKey(:, subj.firstRun:subj.lastRun) + rawStimDuration(eventKey);
    end
    
    rxnEndKey   = stimEndKey + p.rxnWindow; 
    eventEndKey = eventStartKey + p.eventTime;

    % Display instructions
    if Training
        cd(FuncsLoc)
        DisplayInstructions_bkfw_rtbox(Instructions, wPtr, RTBoxLoc); 
        cd(direc)
    end
    DrawFormattedText(wPtr, 'Waiting for first pulse...'); 
    Screen('Flip', wPtr); 

    % Wait for first pulse
    cd(RTBoxLoc)
    RTBox('Clear'); 
    RTBox('UntilTimeout', 1);
    firstPulse(run) = RTBox('WaitTR'); 

    % Draw onto screen after recieving first pulse
    Screen('DrawLines', wPtr, crossCoords, 2, 0, [centerX, centerY]);
    Screen('Flip', wPtr); 

    % Generate absolute time keys
    AbsEvStart   = firstPulse + eventStartKey; 
    AbsStimStart = firstPulse + stimStartKey; 
    AbsStimEnd   = firstPulse + stimEndKey; 
    AbsRxnEnd    = firstPulse + rxnEndKey; 
    AbsEvEnd     = firstPulse + eventEndKey; 

    WaitTill(firstPulse(run) + p.epiTime); 

    %% Present audio stimuli
    try
        for event = 1:p.events
            eventStart(event, run) = GetSecs(); 
            
            PsychPortAudio('FillBuffer', pahandle, audio{eventKey(event)});
            WaitTill(AbsStimStart(event, run)); 

            stimStart(event, run) = GetSecs; 
            PsychPortAudio('Start', pahandle, 1);
            WaitTill(AbsStimEnd(event, run)); 
            stimEnd(event, run) = GetSecs; 
            RTBox('Clear'); 

            [respTime{event, run}, respKey{event, run}] = RTBox(AbsRxnEnd(event, run)); 

            [~, eventEnd(event, run)] = WaitTill(AbsEvEnd(event, run));    
        end
        
        WaitSecs(p.eventTime); 
        runEnd(run) = GetSecs(); 
        

        
    catch err
        sca; 
        runEnd(run) = GetSecs();  %#ok<NASGU>
        cd(FuncsLoc)
        OutputData
        cd(ScriptsLoc)
        rethrow(err)
    end
    
    if run ~= subj.lastRun
    endstring = ['End of run. Nice job! Experimenters, prepare for '...
    	protocol_order{run+1}, num2str(run+1)]; 
	DrawFormattedText(wPtr, endstring); 
	Screen('Flip', wPtr)
    RTBox('Clear'); 
    RTBox(inf); 
    end 

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