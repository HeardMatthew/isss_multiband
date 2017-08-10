%% rhythm_tempo.m
% Script to run for pilot scans testing rhythm and tempo discrimination. 
% Ported from isss_multiband script. 
% Author - Matt Heard

% CHANGELOG
% 07/31/17  Latency issues discovered. Stimuli duration is hard-coded in
%    LoadStimuliAndKeys_rhythm and was incorrect. -- MH
% 08/08/17  Started changelog. -- MH
% 08/08/17  Completed v2. Beginning testing. Work on OutputData. -- MH
% 08/09/17  Continued testing and fixing bugs. Bugs are now all squashed?
%    Ready for experimenting! -- MH

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
AudioDevice = PsychPortAudio('GetDevices', 3); 

% Subject info
prompt = {...
    'Subject number (###):', ...
    'Subject initials (XX):', ...
    'First run (1-2)', ... 
    'Last run (1-2)', ... 
    'Scanner connected (0/1):', ...
    'RTBox connected (0/1):', ...
    }; 
dlg_ans = inputdlg(prompt); 

subj.Num  = dlg_ans{1};
subj.Init = dlg_ans{2}; 
subj.firstRun = str2double(dlg_ans{3}); 
subj.lastRun  = str2double(dlg_ans{4}); 
ConnectedToScanner = str2double(dlg_ans{5});
ConnectedToRTBox   = str2double(dlg_ans{6}); 

% Stimuli
NumberOfRTStimuli = 4;  % 4 different rhythm/tempo clips. Needs to be 4. 
NumberOfStimuli   = 8;  % 8 .wav files in folder (norm + oddball + silent)
% The script presents four sets of the 4 RT stimuli per run in a random 
% order, and includes two oddball and two silent stimuli. 

% Scan
scan.TR = 1.000; 
scan.epiNum = 10; 

% Timing
p.runs        = length(subj.firstRun:subj.lastRun); 
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
expDir = pwd; 

StimuliTag     = [expDir, '\stimuli_rhythm']; 
if mod(str2double(subj.Num), 2) == 1 % Odd numbered subject
    StimuliLoc = [StimuliTag, '\oddsubj'];
else                                  % Even numbered subject
    StimuliLoc = [StimuliTag, '\evensubj'];
end

ScriptsLoc = [expDir, '\scripts'];
ResultsLoc = [expDir, '\results']; 
FuncsLoc   = [ScriptsLoc, '\functions'];

cd ..
RTBoxLoc = [pwd, '\RTBox']; 

%% Preallocating timing variables 
maxNumRuns = 2; 
eventStart    = NaN(p.events, maxNumRuns);
stimStart     = NaN(p.events, maxNumRuns); 
stimEnd       = NaN(p.events, maxNumRuns); 
eventEnd      = NaN(p.events, maxNumRuns); 
stimStartKey  = NaN(p.events, maxNumRuns); 
stimEndKey    = NaN(p.events, maxNumRuns); 
eventStartKey = NaN(p.events, maxNumRuns); 

respTime = cell(p.events, maxNumRuns); 
respKey  = cell(p.events, maxNumRuns); 
    
firstPulse = NaN(1, maxNumRuns);
runEnd     = NaN(1, maxNumRuns); 

%% Prepare experiment
% File names
filetag    = [subj.Num '_' subj.Init '_']; 
ResultsTxt = [filetag 'rhythm_results.txt']; 
ResultsXls = [filetag 'rhythm_results.xlsx']; 
Variables  = [filetag 'rhythm_variables.mat']; 

% Load stimuli and check counterbalance
cd(FuncsLoc) 
[audio, fs, rawStimDur, jitterKey, eventKey, answerKey, rhythmKey] = ...
    LoadStimAndKeys_rhythm(StimuliLoc, p.events, subj.firstRun, subj.lastRun, maxNumRuns);
fs = fs{1}; % Above func checks that all fs are the same. 
cd(expDir)

% Open PTB screen on scanner, prepare fixation cross coords
[wPtr, rect]  = Screen('OpenWindow', 0, 185);
frameDuration = Screen('GetFlipInterval', wPtr);

font = 'Cambria'; % The best font. 
fontsize = 72;

Screen('TextFont',  wPtr, font); 
Screen('TextSize',  wPtr, fontsize);
Screen('TextStyle', wPtr, 2);

centerX = rect(3)/2;
centerY = rect(4)/2;
crossCoords = [-30, 30, 0, 0; 0, 0, -30, 30]; 
HideCursor(); 

% Open audio connection
pahandle = PsychPortAudio('Open', [], [], [], fs);
% Play a silent audio clip just to eliminate clipping on first presentation
PsychPortAudio('FillBuffer', pahandle, audio{8});
PsychPortAudio('Start', pahandle, 1);
    % Stimuli are presented on the scanner computer, which shares a screen
    % and audio output with the scanner projector and headphones. 

% Check if using RTBox or Keyboard
if ConnectedToRTBox == 0
    cd(RTBoxLoc)
    RTBox('fake', 1)
    cd(expDir)
end

%% Prepare each run

for run = subj.firstRun:subj.lastRun
    
    DrawFormattedText(wPtr, 'Please wait, preparing run...');
    Screen('Flip', wPtr); 

	% Prepare timing keys
    eventStartKey(:, run) = p.epiTime + [0:p.eventTime:((p.events-1)*p.eventTime)]'; %#ok<NBRAK>
    stimStartKey(:, run)  = eventStartKey(:, run) + jitterKey(:, run); 

    stimEndKey(:, run) = stimStartKey(:, run) + rawStimDur(eventKey(:, run))';
    
    rxnEndKey   = stimEndKey + p.rxnWindow; 
    eventEndKey = eventStartKey + p.eventTime;
    
    % Wait for first pulse
    DrawFormattedText(wPtr, 'Waiting for first pulse...'); 
    Screen('Flip', wPtr); 

    cd(RTBoxLoc)
    RTBox('Clear'); 
    RTBox('UntilTimeout', 1);
    firstPulse(run) = RTBox('WaitTR'); 

    % Draw onto screen after recieving first pulse
    Screen('DrawLines', wPtr, crossCoords, 2, 0, [centerX, centerY]);
    Screen('Flip', wPtr); 

    % Generate absolute time keys
    AbsEvStart   = firstPulse + eventStartKey(:,run); 
    AbsStimStart = firstPulse + stimStartKey(:,run); 
    AbsStimEnd   = firstPulse + stimEndKey(:,run); 
    AbsRxnEnd    = firstPulse + rxnEndKey(:,run); 
    AbsEvEnd     = firstPulse + eventEndKey(:,run); 

    WaitTill(firstPulse(run) + p.epiTime); 

    %% Present audio stimuli

    try
        for event = 1:p.events
            eventStart(event, run) = GetSecs(); 
            
            PsychPortAudio('FillBuffer', pahandle, audio{eventKey(event, run)});
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
        OutputData_rhythm
        cd(ScriptsLoc)
        PsychPortAudio('Close'); 
        rethrow(err)
    end
    
	if run ~= subj.lastRun %#ok<ALIGN>
        endstring = 'End of run. Nice job! Press any button when ready to continue.'; 
        DrawFormattedText(wPtr, endstring); 
        Screen('Flip', wPtr)
        RTBox('Clear'); 
        RTBox(inf); 
    end 
end

%% Closing down
Screen('CloseAll');
PsychPortAudio('Close'); 
cd(ScriptsLoc)
DisableKeysForKbCheck([]); 

%% Save data
cd(FuncsLoc)
disp('Saving data. Please wait...')
OutputData_rhythm
disp('All done!')
cd(ScriptsLoc)