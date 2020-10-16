%% sigtest.,
% Test if the durations of the 1ch noise stimuli are significantly
% different from the durations of the intelligible stimuli. Run this after
% loading all stimuli. 

noiseDuration = stimDuration(end-15:end-8); 

voiceDuration = stimDuration(1:64); 
avgVD = mean(voiceDuration);  

%%

[h, p] = ttest(noiseDuration, avgVD); 