%% isss_multiband.m
% Script to run for pilot scans comparing ISSS, multiband, and the new
% hybrid scanning protocol. Ported from my previous ISSS_test script. 
% Author - Matt Heard
%% Startup
sca; PsychPortAudio('Close'); clearvars; clc;
codeStart = GetSecs(); 
cd ..
direc = pwd; 
InitializePsychSound
%% Parameters
% Many of these parameters need to be determined while testing. 
AudioDevice = PsychPortAudio('GetDevices', 3); 
    % Is this PsychPortAudio call correct? I want to call all ASIO devices
    % on the Windows computer. In the lab and on my laptop, there is only
    % one. Are there multiple in CCBBI? If so, then I'll have to change
    % which device is called. 
    
prompt = {'Subject number (###):', 'Subject initials (XX):', ...
    'Scan protocol (isss/multi/hybrid):', 'Show instructions (0/1):'}; 
dlg_ans = inputdlg(prompt); 

% Convert dlg_ans into my p.arameters struct
p.subjNum  = dlg_ans{1};
p.subjInit = dlg_ans{2}; 
p.scanType = dlg_ans{3}; 
ShowInstructions = str2double(dlg_ans{4});

    % Go over these numbers with Dr. Lee
if strcmp(p.scanType, 'isss')
    p.TR = 1.000;
    p.epiNum = 10;
elseif strcmp(p.scanType, 'multi')
    p.TR = 1.000;
    p.epiNum = 10; 
elseif strcmp(p.scanType, 'hybrid')
    p.TR = 1.000; 
    p.epiNum = 10;
else
    error('Invalid scan protocol')
end
    
p.runs = 1;                 % 1
p.eventsPerRun = 40;        % 40 events/run
NumberOfStimuli = 40;       % 40 .wav files
p.silentTime = 4.000;       % 4 seconds
p.epiTime = 1.000;          % 10 seconds
p.jitterWindow = 1.000;     % 1 second, see notes below
% For this experiment, the first second of the silent window will not have
% stimuli presented. To code for this, I add 1 second to the jitterKey. 

% Estimates to compare later
p.eventTime = p.silentTime + p.epiTime;
p.runDuration = p.eventTime * p.eventsPerRun; 

triggerCode = KbName('5%'); % This is the trigger recieved from the MRI. 
button1 = KbName('1!'); % Advance backwards through instructions/Male.  
button2 = KbName('2@'); % Advance forwards through instructions/Female.

StimuliLoc = [direc, '\stimuli'];
ScriptsLoc = [direc, '\scripts'];
FuncsLoc = [ScriptsLoc, '\functions'];
Instructions = 'instructions.txt';
Results = [p.subjNum '_' p.subjInit '_' p.scanType '_results.txt']; 
Variables = [p.subjNum '_' p.subjInit '_' p.scanType '_variables.mat']; 

    % Determine these params when testing. 
ScreenNumber = 0; % Which screen does the subject see?
HeadphonesID = 2; % What is the ID of the MRI-safe headphones?

% Debugging
ConnectedToScanner = 0;

% Preallocating certain variables, and preparing the event and jitter keys. 
runStart = cell(1, p.runs);
runEnd= cell(1, p.runs);
runDuration = cell(1, p.runs);
eventStart = NaN(p.runs, p.eventsPerRun); 
eventEnd = NaN(p.runs, p.eventsPerRun); 
epiEnd = NaN(p.runs, p.eventsPerRun); 

%% Prepare test
% Load Stimuli
cd(FuncsLoc)
[audio, fs, jitterKey, eventKey, answerKey] = ...
    LoadStimuli_mono_anskey(StimuliLoc, p.runs, p.eventsPerRun);
fs = fs{1};
cd(direc)

% Open PTB screen on scanner, prepare fixation cross
[wPtr, rect] = Screen('OpenWindow', ScreenNumber, 185);
frameDuration = Screen('GetFlipInterval', wPtr);
centerX = rect(3)/2;
centerY = rect(4)/2;
crossCoords = [-50, 50, 0, 0; 0, 0, -50, 50]; 

% Open audio connection
pahandle = PsychPortAudio('Open', HeadphonesID, [], [], fs);

% Instructions
if ShowInstructions == 1
    cd(FuncsLoc)
    DisplayInstructions_bkfw(Instructions, wPtr, button1, button2); 
    cd(direc)
end

%% Present audio stimuli
% I have coded for the program to present the stimuli and wait. 
for i = 1:p.runs

    % Wait for first pulse
    cd(FuncsLoc)
	WaitForScannerTrigger(ConnectedToScanner, wPtr, triggerCode);
	% WaitSecs(p.epiTime) 
	% Again, uncomment if first pulse is multiple EPIs
    
    % Draw onto screen
    % DrawFormattedText(wPtr, 'Presenting stimuli...', 'center', 'center');
    Screen('DrawLines', wPtr, crossCoords, 2, 0, [centerX, centerY]);
    [~, runStart{i}] = Screen('Flip', wPtr); 
    
    for j = 1:p.eventsPerRun
        PsychPortAudio('FillBuffer', pahandle, audio{eventKey{i}(j)});
        eventStart(i, j) = GetSecs(); 
        
        WaitSecs(jitterKey{i}(j)); 
        PsychPortAudio('Start', pahandle);
        WaitSecs(p.silentTime - jitterKey{i}(j));
        PsychPortAudio('Stop', pahandle, 1);
        eventEnd(i, j) = GetSecs();
        
        WaitSecs(p.epiTime); % Here's where the serious changes need to begin.
        % Maybe a while loop will work, with a break clause if the subject
        % takes too long?
        % Make sure to recover RT and response! Later, compare to correct
        % answer. 
        % Pretty sure RTBox scripts will be the best for this. 
        epiEnd(i, j) = GetSecs();     
    end
    
    if i~= p.runs
        DrawFormattedText(wPtr, 'End of run. Press button2 to continue.',...
            'center', 'center'); % EoR prompt. Gives break to subjects. 
        [~, runEnd{i}] = Screen('Flip', wPtr);

        while 1 % Wait for subject to press button to advance. 
            [keyIsDown, ~, keyCode] = KbCheck(-1);
            if keyIsDown
                if find(keyCode) == button2
                    break
                end
            end
        end

    else
        DrawFormattedText(wPtr, 'End of experiment.', 'center', 'center');
        [~, runEnd{i}] = Screen('Flip', wPtr);
        WaitSecs(3); 
    end

end

%% Saving relevant timing information
eventDuration = eventEnd - eventStart;
epiStart = eventEnd;
epiDuration = epiEnd - epiStart; 
totalDuration = eventDuration + epiDuration; 

cd(direc)
fid = fopen(Results, 'w');    
dstring = ''; 
fstring = '';

for i = 1:p.eventsPerRun
    dstring = strcat(dstring, ' %d '); 
    fstring = strcat(fstring, ' %f ');
end

for i = 1:p.runs
    runDuration{i} = runEnd{i} - runStart{i};     
    fprintf(fid, 'DATA FOR RUN %d ---------- \n', i);
    
    fprintf(fid, 'Run started %6.2f after code started \n', ...
        runStart{i} - codeStart); 
    
    fprintf(fid, 'Run duration: %6.2f \n', runDuration{i});
    fprintf(fid, 'Expected run duration: %6.2f \n', p.runDuration); 

    keystring = ['Event key: ', dstring, '\n'];
    fprintf(fid, keystring, eventKey{i});
    
    jitterstring = ['Jitter key (msec): ', fstring, '\n']; 
    fprintf(fid, jitterstring, (jitterKey{i} * 1000)); 
    
    eventdurationstring = ['Silence durations: ', fstring, '\n'];
    fprintf(fid, eventdurationstring, eventDuration(i,:));
    
    epidurationstring = ['EPI durations: ', fstring, '\n'];
    fprintf(fid, epidurationstring, epiDuration(i, :));
    
    totaldurationstring = ['Total durations: ', fstring, '\n'];
    fprintf(fid, totaldurationstring, totalDuration(i, :));
    
    fprintf(fid, 'Expected event duration: %f \n', p.eventTime); 
    
    fprintf(fid, '\n'); 
end

fclose(fid); 
save(Variables); 

% Closing down
Screen('CloseAll');
PsychPortAudio('Close'); 
cd(ScriptsLoc)