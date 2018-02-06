function [onsets] = PTB_Protocol_OnsetExtract(Data, tpb)
%[onsets] = PTB_ProtocolFile_v1(subjectNum, gainsloss)
%Extract onset times for each trial from PTB Data files, structured as 'Start - Response - Feedback - ITI - End'
%All trials included, did NOT delete the first trial

%INPUTS:
%Data - a struct loaded from a subject file
%tpb - trials per block

%OUTPUT:
%onsets - A structure that includes a matrix of various onsets times
%   .Description - Headings for block onset time matrices
%   .b1 - Onset times for block 1 by trial
%   .b2 - Onset times for block 2 by trial
%   .b3 - Onset times for block 3 by trial
%   .b4 - Onset times for block 4 by trial

% Initialize fields

trialNum = length(Data.trialTime); % possible that some blocks are not done

% column 1-3: hour, min, sec; column 4: convert 1-3 into sec; column 5:
% relative time in sec to the first trial in the block 
start = zeros(trialNum, 5);
ends = zeros(trialNum, 5);
resp = zeros(trialNum, 5);
feedback = zeros(trialNum, 5);
ITIs = zeros(trialNum, 5);

% NOTE: for-loop necessary due to odd nested struct issues
for i = 1:trialNum
  % Extract onset times and convert from hr:min:sec into seconds
  % NOTE: To avoid repetition of indices: `var(i, 4) = var(i, 1:3) * [3600; 60; 1]`
  start(i, 1:3) = Data.trialTime(i).trialStartTime(4:6);
  start(i, 4) = start(i, 1)*3600 + start(i, 2)*60 + start(i, 3);

  % NOTE: Rarely, `ITIStartTime` and `trialEndTime` are not recorded 
  % -> get next trial's trialStartTime and settle for NaN
  try
    ends(i, 1:3) = Data.trialTime(i).trialEndTime(4:6);
    ITIs(i, 1:3) = Data.trialTime(i).ITIStartTime(4:6);
  catch ME
    fprintf('Error caught: %s in %s at choice %d\n', ME.identifier, Data.filename, i) % MATLAB.badsubscript
    if (floor(i / tpb) * tpb + 1) == (i + 1) % if next choice is a different block
      ends(i, 1:3) = NaN;
    else
      ends(i, 1:3) = Data.trialTime(i + 1).trialStartTime(4:6);
    end
    ITIs(i, 1:3) = NaN;
  end
  ends(i, 4) = ends(i, 1)*3600 + ends(i, 2)*60 + ends(i, 3);
  ITIs(i, 4) = ITIs(i, 1)*3600 + ITIs(i, 2)*60 + ITIs(i, 3);

  resp(i, 1:3) = Data.trialTime(i).respStartTime(4:6);
  resp(i, 4) = resp(i, 1)*3600 + resp(i, 2)*60 + resp(i, 3);

  feedback(i, 1:3) = Data.trialTime(i).feedbackStartTime(4:6);
  feedback(i, 4) = feedback(i, 1)*3600 + feedback(i, 2)*60 + feedback(i, 3);


  % Calculate time in relation to the appropriate block's first trial onset
  if i < tpb+1
    first_trial = 1;
  elseif i < 2*tpb+1
    first_trial = tpb+1;
  elseif i < 3*tpb+1
    first_trial = 2*tpb+1;
  else
    first_trial = 3*tpb+1;
  end
  % NOTE: `floor(i / tpb) * tpb + 1` designates each block's first trial for
  % an arbitrary number of blocks, but this is more straightforward.

  start(i, 5) = start(i, 4) - start(first_trial, 4);
  ends(i, 5) = ends(i, 4) - start(first_trial, 4);
  resp(i, 5) = resp(i, 4) - start(first_trial, 4);
  feedback(i, 5) = feedback(i, 4) - start(first_trial, 4);
  ITIs(i, 5) = ITIs(i, 4) - start(first_trial, 4);
end

%% Take the time in seconds and divide into blocks
% each row is a trial in the block, each column is the time point in sec
% structured as 'Start - Response - Feedback - ITI - End'
block1 = [start(1:tpb, 5) resp(1:tpb, 5) feedback(1:tpb, 5) ITIs(1:tpb, 5) ends(1:tpb, 5)];
block2 = [start(tpb+1:2*tpb, 5) resp(tpb+1:2*tpb, 5) feedback(tpb+1:2*tpb, 5) ITIs(tpb+1:2*tpb, 5) ends(tpb+1:2*tpb, 5)];
block3 = [start(2*tpb+1:3*tpb, 5) resp(2*tpb+1:3*tpb, 5) feedback(2*tpb+1:3*tpb, 5) ITIs(2*tpb+1:3*tpb, 5) ends(2*tpb+1:3*tpb, 5)];
block4 = [start(3*tpb+1:4*tpb, 5) resp(3*tpb+1:4*tpb, 5) feedback(3*tpb+1:4*tpb, 5) ITIs(3*tpb+1:4*tpb, 5) ends(3*tpb+1:4*tpb, 5)];

%% Return
% FirstTrial1 = Data.trialTime(1).trialStartTime(1:6);  % First Day 1 trial
% FirstTrial2 = Data.trialTime(63).trialStartTime(1:6); % First Day 2 trial

% onsets.Day1 = [num2str(FirstTrial1(4)) ':' num2str(FirstTrial1(5)) ':' num2str(FirstTrial1(6)) ', ' num2str(FirstTrial1(2)) '/' num2str(FirstTrial1(3)) '/' num2str(FirstTrial1(1))];
% onsets.Day2 = [num2str(FirstTrial2(4)) ':' num2str(FirstTrial2(5)) ':' num2str(FirstTrial2(6)) ', ' num2str(FirstTrial2(2)) '/' num2str(FirstTrial2(3)) '/' num2str(FirstTrial2(1))];
onsets.Description = 'Start - Response - Feedback - ITI - End';
onsets.b1 = block1;
onsets.b2 = block2;
onsets.b3 = block3;
onsets.b4 = block4;
end
