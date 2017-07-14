% WaitForScannerTrigger.m
function WaitForScannerTrigger(ScannerConnected, windowPointer, fmriCode)
% First image is taken when scanner sends a '5%' to the computer. This 
% script pauses the rest of the experiment from starting until it recieves 
% the first EPI. Adapted from Jonas Kaplan fMRI 
% tutorial by Matt H. 
if ScannerConnected == 0
    DrawFormattedText(windowPointer, 'Press 5% to continue', 'center', 'center');
    Screen('Flip', windowPointer);
end

% Find the key code for the fMRI scanner
keyIsDown = 0; %#ok<NASGU>

% Make sure no keys are disabled for key check
DisableKeysForKbCheck([]); 

% Get the time when program starts to wait. Not needed?
% startToWait = GetSecs();  %#ok<NASGU>

% Wait for the first trigger

while 1 % Wait for a certain # of EPIs to be taken. 
    [keyIsDown, ~, keyCode] = KbCheck(-1);
    if keyIsDown
        if find(keyCode) == fmriCode
            break
        end
    end
end
% Uncomment the following line to ignore 5s for the rest of the experiment.
DisableKeysForKbCheck([fmriCode]); %#ok<NBRAK>
end