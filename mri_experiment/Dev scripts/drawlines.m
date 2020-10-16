% drawlines.m
[wPtr, rect] = Screen('OpenWindow', 0, 185);
frameDuration = Screen('GetFlipInterval', wPtr);
centerX = rect(3)/2;
centerY = rect(4)/2;
crossCoords = [-25, 25, 0, 0; 0, 0, -25, 25]; 
Screen('DrawLines', wPtr, crossCoords, 2, 0, [centerX, centerY]);
Screen('Flip', wPtr); 
WaitSecs(2); 
Screen('CloseAll'); 