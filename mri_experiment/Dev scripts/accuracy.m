%% Accuracy
% Subject accuracy on LDT

%% Params
numRuns = length(subj.firstRun:subj.lastRun); 

%% Convert to vectors
answerVec   = []; 
responseVec = [];

for run = subj.firstRun:subj.lastRun
    answerVec   = vertcat(answerKey(:, run), answerVec); 
    responseVec = vertcat(response(:, run), responseVec); 
end

%% Pop out noise trials
for answer = 1:length(answerVec)
    if answerVec(answer) == 3
        answerVec(answer)   = NaN; 
        responseVec(answer) = NaN; 
    end
end

answerVec(isnan(answerVec))     = []; 
responseVec(isnan(responseVec)) = []; 

%% Calculate accuracy
correct = answerVec == responseVec; 
acc = 100 * sum(correct) / length(correct); 