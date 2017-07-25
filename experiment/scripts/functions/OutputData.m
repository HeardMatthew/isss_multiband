%% OutputData.m
% Saves data into a txt and xlsx file. Run as part of the main experiment 
% script. Author -- Matt H

%% Saving relevant timing information
% Convert to relative time, instead of system
actualJitter        = stimStart - eventStart; 
actualStimDuration  = stimEnd - stimStart; 
actualEventDuration = eventEnd - eventStart; 
stimDurationPrint   = NaN(p.events, n); 
if Training 
    speechKeyPrint  = NaN(5, 1); % This might be wrong, check later
else
    speechKeyPrint  = NaN(8, n); 
end
eventKeyPrint       = NaN(p.events, n); 
answerKeyPrint      = NaN(p.events, n); 

% Print variables have the required size, even if the run nums are weird
stimDurationPrint(:, subj.firstRun:subj.lastRun) = stimDuration(eventKey); 
speechKeyPrint(:, subj.firstRun:subj.lastRun) = speechKey; 
eventKeyPrint(:, subj.firstRun:subj.lastRun) = eventKey; 
answerKeyPrint(:, subj.firstRun:subj.lastRun) = answerKey; 

if Training
    stimDurationTemp = stimDuration(eventKey)'; 
else
    stimDurationTemp = stimDuration(eventKey); 
end
    
headers = {'Jitter key', 'Actual jitter', ... 
        'Stim duration key', 'Actual stim duration', ...
        'Event duration', 'Event key', 'Answer key', ... 
        'Subj response', 'RT'}; 

runDuration = runEnd - firstPulse; 

reaction = NaN(p.events, 6); 
response = NaN(p.events, 6); 

% Convert keys from cells to vectors
for i = 1:p.events
    for j = 1:size(respKey, 2)
        if isempty(respKey{i, j})
            respKey{i, j}  = '0'; 
            respTime{i, j} = NaN; 
        end
        reaction(i, j) = respTime{i, j} - stimEnd(i, j); 
        response(i, j) = str2double(respKey{i, j});
    end
end

% Path
cd(ResultsLoc)
if ~exist(subj.Num, 'file')
    mkdir(subj.Num); 
end
cd(subj.Num)

%% Checks if files already exists to prevent overwrite
while exist(ResultsXls, 'file') == 2
	ResultsXls = [ResultsXls(1:end-5), '_new', ResultsXls(end-4:end)]; 
end

while exist(ResultsTxt, 'file') == 2
	ResultsTxt = [ResultsTxt(1:end-4), '_new', ResultsTxt(end-3:end)]; 
end

%% Prepare to print to txt file
fid = fopen(ResultsTxt, 'w');    
dstring = '';
fstring = '';
smalldstring = ''; 

for i = 1:p.events
    dstring = strcat(dstring, ' %d '); 
    fstring = strcat(fstring, ' %f ');
end

for i = 1:NumSpeechStimuli/4
    smalldstring = strcat(smalldstring, ' %d ');
end


%% Begin printing to txt file
for run = subj.firstRun:subj.lastRun
    fprintf(fid, 'DATA FOR RUN %d ---------- \n', run);
    fprintf(fid, 'SCAN TYPE: %6s \n', protocol_order{run}); 
    
    % Timing data
    fprintf(fid, '# Timing data \n'); 
    fprintf(fid, 'Run started %6.2f after code started \n', ...
        firstPulse(run) - codeStart); 
    fprintf(fid, 'Run duration: %6.2f \n', runDuration(run));
    fprintf(fid, 'Expected run duration: %6.2f \n', p.runDuration); 

    jitterstring = ['Jitter key (msec): ', fstring, '\n']; 
    fprintf(fid, jitterstring, (jitterKey(:, run) * 1000)); 

    actualjitterstring = ['Actual jitter (msec): ', fstring, '\n']; 
    fprintf(fid, actualjitterstring, (actualJitter(:, run) * 1000)); 

    stimulistring = ['Stimuli duration key: ', fstring, '\n'];
    fprintf(fid, stimulistring, stimDurationPrint(:, run));

    actualstimulistring = ['Actual stimuli duration: ', fstring, '\n']; 
    fprintf(fid, actualstimulistring, actualStimDuration(:, run)); 

    eventdurationstring = ['Event duration: ', fstring, '\n']; 
    fprintf(fid, eventdurationstring, actualEventDuration(:,run)); 

    fprintf(fid, 'Expected total duration: %f \n', p.eventTime); 

    % Stimuli data
    fprintf(fid, '# Stimuli data \n'); 
    speechstring = ['Speech samples used: ', smalldstring, '\n']; 
    fprintf(fid, speechstring, speechKeyPrint(:,run)); 

    keystring = ['Event key: ', dstring, '\n'];
    fprintf(fid, keystring, eventKeyPrint(:,run));

    ansstring = ['Answer key: ', dstring, '\n'];
    fprintf(fid, ansstring, answerKeyPrint(:,run)); 

    % Subject data
    fprintf(fid, '# Response data \n'); 

    respstring = ['Subject responses: ', dstring, '\n']; 
    fprintf(fid, respstring, response(:, run)); 

    reactionstring = ['Reaction time (msec): ', fstring, '\n'];
    fprintf(fid, reactionstring, (reaction(:,run) * 1000)); 
    
    fprintf(fid, '\n'); 
    
    %% Prepare to print to xlsx file
    data = cell(p.events + 1, 9); 
    
    M    = horzcat(jitterKey(:,run), actualJitter(:,run), ...
        stimDurationPrint(:,run), actualStimDuration(:,run), ... 
        actualEventDuration(:,run), ...
        eventKeyPrint(:,run), answerKeyPrint(:,run), ...
        response(:,run), reaction(:,run)); 

    data(1,:) = headers; 
    for i = 1:p.events
        for j = 1:9
            data{i+1, j} = M(i, j); 
        end
    end
    
    %% Print to xlsx file
    warning off
    if run ~= 1
        for i = 1:run - 1
            if strcmp(protocol_order{run}, protocol_order{i})
                protocol_order{run} = [protocol_order{run}, '2']; 
            end
        end
    end
    
    xlswrite(ResultsXls, data, protocol_order{run})
    warning on
    
end

%% All done!
fclose(fid); 
save(Variables); 