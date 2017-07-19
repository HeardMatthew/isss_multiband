%% OutputData.m
% Saves data into a txt and xlsx file. Run as part of the main experiment 
% script. Author -- Matt H

%% Saving relevant timing information
% Convert to relative time, instead of system
actualJitter        = stimStart - eventStart; 
actualStimDuration  = stimEnd - stimStart; 
actualEventDuration = eventEnd - eventStart; 

runDuration = runEnd - firstPulse; 

reaction = NaN(1, p.events); 
response = NaN(1, p.events); 

% Convert keys from cells to vectors
for i = 1:p.events
    if isempty(respKey{i})
        respKey{i}  = '0'; 
        respTime{i} = NaN; 
    end
    reaction(i) = respTime{i} - stimEnd(i); 
    response(i) = str2double(respKey{i});
end

% Path
cd(ResultsLoc)
if ~exist(p.subjNum, 'file')
    mkdir(p.subjNum); 
end
cd(p.subjNum)

%% Prepare to print to txt file
fid = fopen(ResultsTxt, 'w');    
dstring = '';
fstring = '';

for i = 1:p.events
    dstring = strcat(dstring, ' %d '); 
    fstring = strcat(fstring, ' %f ');
end

%% Begin printing to txt file
fprintf(fid, 'DATA FOR RUN ---------- \n');

% Timing data
fprintf(fid, '# Timing data \n'); 
fprintf(fid, 'Run started %6.2f after code started \n', ...
    firstPulse - codeStart); 
fprintf(fid, 'Run duration: %6.2f \n', runDuration);
fprintf(fid, 'Expected run duration: %6.2f \n', p.runDuration); 

jitterstring = ['Jitter key (msec): ', fstring, '\n']; 
fprintf(fid, jitterstring, (jitterKey * 1000)); 

actualjitterstring = ['Actual jitter (msec): ', fstring, '\n']; 
fprintf(fid, actualjitterstring, (actualJitter* 1000)); 

stimulistring = ['Stimuli duration key: ', fstring, '\n'];
fprintf(fid, stimulistring, stimDuration(eventKey));

actualstimulistring = ['Actual stimuli duration: ', fstring, '\n']; 
fprintf(fid, actualstimulistring, actualStimDuration); 

eventdurationstring = ['Event duration: ', fstring, '\n']; 
fprintf(fid, eventdurationstring, actualEventDuration); 

fprintf(fid, 'Expected total duration: %f \n', p.eventTime); 

% Stimuli data
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

%% Prepare to print to xlsx file
data = cell(p.events + 1, 9); 
M    = horzcat(jitterKey', actualJitter', ...
    stimDuration(eventKey)', actualStimDuration', ... 
    actualEventDuration', ...
    eventKey', answerKey', ...
    response', reaction'); 

headers = {'Jitter key', 'Actual jitter', ... 
        'Stim duration key', 'Actual stim duration', ...
        'Event duration', 'Event key', 'Answer key', ... 
        'Subj response', 'RT'}; 
    
data(1,:) = headers; 
for i = 1:p.events
    for j = 1:9
        data{i+1, j} = M(i, j); 
    end
end

%% Print to xlsx file
xlswrite(ResultsXls, data)

%% All done!
save(Variables); 