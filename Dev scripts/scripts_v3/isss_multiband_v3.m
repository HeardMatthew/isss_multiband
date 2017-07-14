%% isss_multiband.m
% Script to run for pilot scans comparing ISSS, multiband, and the new
% hybrid scanning protocol. Ported from my previous ISSS_test script. 
% Author - Matt Heard
%% Startup
sca; DisableKeysForKbCheck([]); 
Screen('Preference','VisualDebugLevel', 0); 

try 
    PsychPortAudio('Close'); 
catch
    disp('PsychPortAudio is already closed.')
end

clearvars; clc; 
codeStart = GetSecs(); 
InitializePsychSound

%% Parameters
% Many of these parameters need to be determined while testing. 
AudioDevice = PsychPortAudio('GetDevices', 3); 
    
prompt = {...
    'Subject number (###):', ...
    'Subject initials (XX):', ...
    'Scan protocol (isss/multi/hybrid):', ...
    'Show instructions (0/1):', ...
    'Scanner connected (0/1):', ...
    'RTBox connected (0/1):', ...
    }; 
dlg_ans = inputdlg(prompt); 

% Convert dlg_ans into my p.arameters struct
p.subjNum  = dlg_ans{1};
p.subjInit = dlg_ans{2}; 
p.scanType = dlg_ans{3}; 
ShowInstructions   = str2double(dlg_ans{4});
ConnectedToScanner = str2double(dlg_ans{5});
ConnectedToRTBox   = str2double(dlg_ans{6}); 

if strcmp(p.scanType, 'isss')
    p.TR = 2.000; % Either 2.000 or 2.500 s (7/3)
    p.epiNum = 5; % Depends on TR
elseif strcmp(p.scanType, 'multi')
    p.TR = 1.000; 
    p.epiNum = 10; 
elseif strcmp(p.scanType, 'hybrid')
    p.TR = 1.000; 
    p.epiNum = 10;
else
    error('Invalid scan protocol')
end

NumberOfSpeechStimuli = 64; % 64 different speech clips
NumberOfStimuli       = 80; % 80 .wav files in stimuli folder
% Of these 80 .wav files, 8 of them are silent, 8 of them are noise, and 64
% of them are speech sounds. Of these 32 speech sounds, this script chooses
% 16 (4 MO, 4 FO, 4 MS, 4 FS) for presentation. Subjects never hear a
% repeated "sentence structure" i.e. one stimuli of of 001, one stimuli of
% 002... 

p.events     = 32;     % 32 events
p.stimTime   = 4.000;  % 4 seconds
p.epiTime    = 10.000; % 10 seconds
p.jitWindow  = 1.000;  % 1 second, see notes below
p.respWindow = 3.000;  % 3 seconds
    % For this experiment, the first second of the silent window will not
    % have stimuli presented. To code for this, I add 1 second to the 
    % jitterKey. This variable is just here to remind me that the addition 
    % occurs later, within the LoadStimuli function. 

% Buttons
triggerCode = KbName('5%'); % This is the trigger recieved from the MRI. 
button1     = KbName('1!'); % Advance backwards through instructions/Male.  
button2     = KbName('2@'); % Advance forwards through instructions/Female.
escape      = KbName('esc'); 

% Paths
cd ..
direc = pwd; 

StimuliLoc   = [direc, '\stimuli'];
ScriptsLoc   = [direc, '\scripts'];
ResultsLoc   = [direc, '\results']; 
FuncsLoc     = [ScriptsLoc, '\functions'];
Instructions = 'instructions.txt';

cd ..
RTBoxLoc     = [pwd, '\RTBox']; 

Results    = [p.subjNum '_' p.subjInit '_' p.scanType '_results.txt']; 
Variables  = [p.subjNum '_' p.subjInit '_' p.scanType '_variables.mat']; 

% Preallocating timing variables
presentStart = NaN(1, p.events); 
presentEnd   = NaN(1, p.events); 
epiEnd    = NaN(1, p.events); 
RTEnd     = NaN(1, p.events); 

respTimeKey = cell(1, p.events); 
responseKey = cell(1, p.events); 
    % I use cells here so that responses can be empty when subjects time out,
    % and because responses from RTBox come back as strings. 

%% Prepare test
% Load Stimuli
cd(FuncsLoc) 
[audio, fs, stimDuration, jitterKey, eventKey, answerKey, speechKey] = ...
    LoadStimuliAndKeys(StimuliLoc, p.events, NumberOfSpeechStimuli);
fs = fs{1}; % Above func checks that all fs are the same. 

cd(FuncsLoc)
stimulicheck(NumberOfSpeechStimuli, eventKey); 
cd(direc)

% Open PTB screen on scanner, prepare fixation cross coords
[wPtr, rect] = Screen('OpenWindow', 0, 185);
HideCursor(); 
frameDuration = Screen('GetFlipInterval', wPtr);
centerX = rect(3)/2;
centerY = rect(4)/2;
crossCoords = [-20, 20, 0, 0; 0, 0, -20, 20]; 

% Open audio connection
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
if ShowInstructions == 1
    cd(FuncsLoc)
    DisplayInstructions_bkfw_rtbox(Instructions, wPtr, RTBoxLoc); 
    cd(direc)
end
DrawFormattedText(wPtr, 'Waiting for experimenters...'); 
Screen('Flip', wPtr); 

% Wait for first pulse
cd(RTBoxLoc)
RTBox('Clear'); 
firstpulse = RTBox('WaitTR'); 

% Draw onto screen after recieving first pulse
Screen('DrawLines', wPtr, crossCoords, 2, 0, [centerX, centerY]);
[~, runStart] = Screen('Flip', wPtr); 

WaitSecs(p.epiTime); 

%% Present audio stimuli
for j = 1:p.events
    RTBox('Clear'); 
    PsychPortAudio('FillBuffer', pahandle, audio{eventKey(j)});
    presentStart(j) = GetSecs(); 

    WaitSecs(jitterKey(j)); 
    PsychPortAudio('Start', pahandle);
    WaitSecs(stimDuration(eventKey(j)));
    PsychPortAudio('Stop', pahandle, 1);
    presentEnd(j) = GetSecs();

    [respTimeKey{j}, responseKey{j}] = RTBox(p.respWindow); 
    if ~strcmp(responseKey{j}, '')
        RTEnd(j) = WaitSecs('UntilTime', presentEnd(j) + p.respWindow); 
    end
    
    epiEnd(j) = WaitSecs('UntilTime', presentStart(j) + p.stimTime + p.epiTime);    
end
WaitSecs(p.stimTime + p.epiTime); 

DrawFormattedText(wPtr, 'End of run. Saving results.', 'center', 'center');
[~, runEnd] = Screen('Flip', wPtr);
WaitSecs(3); 

%% Saving relevant timing information
% Creates comparisons
p.eventTime    = p.stimTime + p.epiTime;
p.runDuration  = p.epiTime + ... % After first pulse
    p.eventTime * p.events + ... % Each event
    p.stimTime  + p.epiTime;     % After last acquisition

% Convert to relative time, instead of system
presentDuration = presentEnd - presentStart;
epiStart        = presentStart + p.stimTime; % Approximated
epiDuration     = epiEnd - epiStart;         % Approximated 
totalDuration   = epiEnd - presentStart; 
runDuration     = runEnd - runStart; 

reaction = NaN(1, p.events); 
response = NaN(1, p.events); 

% Convert keys from cells to vectors
for i = 1:p.events
    if isempty(respTimeKey{i})
        respTimeKey{i} = NaN; 
    elseif isempty(responseKey{i})
        responseKey{i} = 'NaN'; 
    end
    reaction(i) = respTimeKey{i} - presentEnd(i); 
    response(i) = str2double(responseKey{i});
end

% Path
cd(ResultsLoc)
mkdir(p.subjNum); cd(p.subjNum); 

% Prepare to print to file
fid = fopen(Results, 'w');    
dstring = '';
fstring = '';
smalldstring = ''; 

for i = 1:p.events
    dstring = strcat(dstring, ' %d '); 
    fstring = strcat(fstring, ' %f ');
end

for i = 1:NumberOfSpeechStimuli/4
    smalldstring = strcat(smalldstring, ' %d ');
end

%% Begin printing
fprintf(fid, 'DATA FOR RUN ---------- \n');

% Timing data
fprintf(fid, '# Timing data \n'); 
fprintf(fid, 'Run started %6.2f after code started \n', ...
    runStart - codeStart); 
fprintf(fid, 'Run duration: %6.2f \n', runDuration);
fprintf(fid, 'Expected run duration: %6.2f \n', p.runDuration); 

jitterstring = ['Jitter key (msec): ', fstring, '\n']; 
fprintf(fid, jitterstring, (jitterKey * 1000)); 

presentdurationstring = ['Presentation durations: ', fstring, '\n'];
fprintf(fid, presentdurationstring, presentDuration);

epidurationstring = ['Approximated EPI durations: ', fstring, '\n']; 
fprintf(fid, epidurationstring, epiDuration); 

totaldurationstring = ['Total durations: ', fstring, '\n'];
fprintf(fid, totaldurationstring, totalDuration);

fprintf(fid, 'Expected total duration: %f \n', p.eventTime); 

% Stimuli data
fprintf(fid, '# Stimuli data \n'); 
speechstring = ['Speech samples used: ', smalldstring, '\n']; 
fprintf(fid, speechstring, speechKey); 

keystring = ['Event key: ', dstring, '\n'];
fprintf(fid, keystring, eventKey);

ansstring = ['Answer key: ', dstring, '\n'];
fprintf(fid, ansstring, answerKey); 

% Subject data
fprintf(fid, '# Response data \n'); 

respstring = ['Subject responses: ', dstring, '\n']; 
fprintf(fid, respstring, response); 

reactionstring = ['Reaction time (msec): ', fstring, '\n'];
fprintf(fid, reactionstring, (reaction * 1000)); 

% Done printing!
fprintf(fid, '\n'); 

fclose(fid); 
save(Variables); 

%% Closing down
Screen('CloseAll');
PsychPortAudio('Close'); 
cd(ScriptsLoc)
DisableKeysForKbCheck([]); 