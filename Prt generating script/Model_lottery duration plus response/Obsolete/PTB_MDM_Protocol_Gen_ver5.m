function PTB_MDM_Protocol_Gen_ver5(subjectNum, domain, tr, trialduration, DiscAcq, ParametricMod, path_in, path_out, PRT, tpb)
%PTB_PROTOCOL_GEN Generates PRT Files from Psychtoolbox Data Files

% TODo: 
% Conditionally use raw value or rating value for RewardValue
% How to calculate Subjective value 

%INPUTS:
%       subjectNum - Specifies subject file to be loaded
%       domain - A string of either 'monetary' or 'medical' to indicate the domain to extract
%       tr - Temporal resolution of the visual data, in seconds
%       trialduration - Number of seconds to analyze after trial onset
%       DiscAcq - Number of seconds to discard at the beginning of each block
%       ParametricMod - which value to use as parameter? Values: 'RewardValue', 'RiskLevel', 'AmbiguityLevel', 'SV', or ''
%       path_in - the dominant folder that both auxiliary function files and data files are stored in
%       path_out - the folder to save the generated PRT files into
%       PRT - a struct of settings for the PRT file
%       tpb - trials per block
%
%OUTPUT: PRT files for given domain in the folder specified by `path_out`
%
% NOTE: In order to create a protocol file with onsets measured in seconds,
% simply switch tr value to 1


%% Process arguments
% Store in logical value whether protocol should be generated for monetary or medical
is_mon = strcmp(domain, 'mon');
if (~is_mon && ~strcmp(domain, 'med'))
  error('Aborting: `domain` must be set to either "mon" or "med"!')
  return
end

% Set NumParametricWeights on the basis of ParametricMod
permissible_parameters = {'RewardValue', 'RiskLevel', 'AmbiguityLevel', 'SubjRating','SV'};
if sum(strcmp(ParametricMod, permissible_parameters)) 
  % using sum because strcmp returns a matrix of logicals if string is compared 
  % with an array of strings
  NumParametricWeights = 1;
elseif strcmp(ParametricMod, '')
  NumParametricWeights = 0;
else
  error('Aborting: `ParametricMod` must be set to a valid parameter value!')
  return
end
PRT.ParametricWeights = num2str(NumParametricWeights);

% For PRT properties, fix colors into correct order
if ~is_mon
  Colors = {PRT.ColorMed_5, PRT.ColorMed_8, PRT.ColorMed_12, PRT.ColorMed_25};
  alt_colors = {PRT.ColorMon_5, PRT.ColorMon_8, PRT.ColorMon_12, PRT.ColorMon_25};
else
  Colors = {PRT.ColorMon_5, PRT.ColorMon_8, PRT.ColorMon_12, PRT.ColorMon_25};
  alt_colors = {PRT.ColorMed_5, PRT.ColorMed_8, PRT.ColorMed_12, PRT.ColorMed_25};
end

%% Load Mon and Med Data files
% Add the directory & all subdirs to path
addpath(genpath(path_in)); 

load(['MDM_MON_' num2str(subjectNum) '_fitpar.mat']);
load(['MDM_MED_' num2str(subjectNum) '_fitpar.mat']);


% Pick the domain to analyze
if (is_mon)
  data = Datamon;
else
  data = Datamed;
end

%% Get correct block order
block_order = getBlockOrder(is_mon, Datamon, Datamed, tpb);

%% Compute subjective value of each choice only for monetary, use constrained fitting

% ToDo: If logic, to use raw or rating value 

if is_mon
    for reps = 1:length(data.choice)
      sv(reps, 1) = ambig_utility(0, ...
          data.vals(reps), ...
          data.probs(reps), ...
          data.ambigs(reps), ...
          data.mfCstr.alpha, ...
          data.mfCstr.beta, ...
          'ambigNrisk');
    end
end
% % Flip sign, since the data files store only value magnitudes 
% if ~is_gains
%   sv(:, 1) = -1 * sv(:, 1);
% end

%% Load onset times
monon = PTB_Protocol_OnsetExtract(Datamon,tpb);
medon = PTB_Protocol_OnsetExtract(Datamed,tpb);

% Extract per-block time info from returned argument
% Each element in the cell is a matrix, each row is a trial, columns structured as 'Start - Response - Feedback - ITI - End'
mononsets = {monon.b1, monon.b2, monon.b3, monon.b4};
medonsets = {medon.b1, medon.b2, medon.b3, medon.b4};

%% Iterate over blocks in domain
% NOTE: 4 is magic number, since currently each domain has exactly 4 blocks
for blocknum = 1:4
  %% Select onset/offset time block to use
  if is_mon
    prtblock = mononsets{blocknum};
  else
    prtblock = medonsets{blocknum};
  end

  %% Process onset/offset times
  % first trial (DiscAcq) is discarded after this
  [ onsets, offsets ] = getBlockOnsetOffset(prtblock, DiscAcq, tr, trialduration,tpb);

  %% Compute indexes from current block
  % Get indices for risk trials and ambiguity trials (exclude first trial in each block)
  current_block_range = (2:tpb) + (blocknum - 1) * tpb;

  % Divide blocks by type of lottery (4 outcome magnitudes)
  % These _indexes are in terms of the first trial deleted!!
  index_5 = data.vals(current_block_range) == 5;
  index_8 = data.vals(current_block_range) == 8;
  index_12 = data.vals(current_block_range) == 12;
  index_25 = data.vals(current_block_range) == 25;

  % Parametric weights for current block
  if NumParametricWeights > 0 & is_mon
    block_amt = data.vals(current_block_range);
    block_rlevel = data.probs(current_block_range);
    block_alevel = data.ambigs(current_block_range);
    block_sv = sv(current_block_range);
    block_rtg = data.subjRatings(current_block_range);
  elseif NumParametricWeights > 0 & is_mon == 0
    block_amt = data.vals(current_block_range);
    block_rlevel = data.probs(current_block_range);
    block_alevel = data.ambigs(current_block_range);
    block_rtg = data.subjRatings(current_block_range);
  end

  %% Select what parametric value to write (or not write) into PRT file
  % Store the basic onset/offset, computed earlier
  block_5 = [onsets(index_5,1) offsets(index_5,1)];
  block_5_resp = [onsets(index_5,2) offsets(index_5,2)];
  
  block_8 = [onsets(index_8,1) offsets(index_8,1)];
  block_8_resp = [onsets(index_8,2) offsets(index_8,2)];

  block_12 = [onsets(index_12,1) offsets(index_12,1)];
  block_12_resp = [onsets(index_12,2) offsets(index_12,2)];

  block_25 = [onsets(index_25,1) offsets(index_25,1)];
  block_25_resp = [onsets(index_25,2) offsets(index_25,2)];


  % Add the selected parametric value if required
  if NumParametricWeights > 0
    if strcmp(ParametricMod, 'RewardValue')& is_mon
      block_amb = [block_amb block_amt(amb_index)]; % onset-offset-paramatricValue
      block_risk = [block_risk block_amt(risk_index)];
      block_amb_resp = [block_amb_resp block_amt(amb_index)];
      block_risk_resp = [block_risk_resp block_amt(risk_index)];
    elseif strcmp(ParametricMod, 'RiskLevel')
      block_amb = [block_amb block_rlevel(amb_index)];
      block_risk = [block_risk block_rlevel(risk_index)];
      block_amb_resp = [block_amb_resp block_rlevel(amb_index)];
      block_risk_resp = [block_risk_resp block_rlevel(risk_index)];
    elseif strcmp(ParametricMod, 'AmbiguityLevel')
      block_amb = [block_amb block_alevel(amb_index)];
      block_risk = [block_risk block_alevel(risk_index)];
      block_amb_resp = [block_amb_resp block_alevel(amb_index)];
      block_risk_resp = [block_risk_resp block_alevel(risk_index)];
    elseif strcmp(ParametricMod, 'SubjRating')
      block_amb = [block_amb block_rtg(amb_index)];
      block_risk = [block_risk block_rtg(risk_index)];
      block_amb_resp = [block_amb_resp block_rtg(amb_index)];
      block_risk_resp = [block_risk_resp block_rtg(risk_index)];
    elseif strcmp(ParametricMod, 'SV') & is_mon
      block_amb = [block_amb block_sv(amb_index)];
      block_risk = [block_risk block_sv(risk_index)];
      block_amb_resp = [block_amb_resp block_sv(amb_index)];
      block_risk_resp = [block_risk_resp block_sv(risk_index)];
    end
  end

  %% Write file to txt
  % TODO: Called once per block - abstract into a subfunction / nested function?

  % Open file for writing
  if strcmp(ParametricMod, '')
    ParametricMod = 'none';
  end
  fname = [path_out num2str(subjectNum) '_block' block_order{blocknum} ...
      '_' domain num2str(blocknum) '_model_dispresp' '_type_' 'Magnitude_NoRiskAmb_' ParametricMod '.prt'];
  fileID = fopen(fname, 'w');
  
  % TODO: Figure out how %Ns works -- it's highly unlikely that both values
  % actually do a useful thing
  fprintf(fileID, '%12s %10s\r\n \r\n', 'FileVersion:', PRT.FileVersion);
  fprintf(fileID, '%17s %11s\r\n \r\n', 'ResolutionOfTime:', PRT.ResolutionOfTime);
  fprintf(fileID, '%11s %22s\r\n \r\n', 'Experiment:', PRT.Experiment);
  fprintf(fileID, '%16s %16s\r\n', 'BackgroundColor:', PRT.BackgroundColor);
  fprintf(fileID, '%10s %16s\r\n', 'TextColor:', PRT.TextColor);
  fprintf(fileID, '%16s %10s\r\n', 'TimeCourseColor:', PRT.TimeCourseColor);
  fprintf(fileID, '%16s %6s\r\n', 'TimeCourseThick:', PRT.TimeCourseThick);
  fprintf(fileID, '%19s %11s\r\n', 'ReferenceFuncColor:', PRT.ReferenceFuncColor);
  fprintf(fileID, '%19s %3s\r\n \r\n', 'ReferenceFuncThick:', PRT.ReferenceFuncThick);

%   if NumParametricWeights > 0
%       fprintf(fileID, '%18s %4s\r\n', 'ParametricWeights:', PRT.ParametricWeights);
%   end
% 
  fprintf(fileID, '\r\n%15s %7s\r\n\r\n', 'NrOfConditions:', PRT.NrOfConditions);

  % NOTE: Empty conditions -- set here for purposes of consistent order. 
  % Equivalent but opposite block further below
%   if ~is_mon & NumParametricWeights ~= 0 % medical block with parametric
%     fprintf(fileID, '%15s\r\n%3s\r\n%4s %6s\r\n\r\n', 'Amb_mon_Display', '0', 'Color:', alt_colors{1});
%     fprintf(fileID, '%20s\r\n%3s\r\n%4s %6s\r\n\r\n', 'Amb_mon_Display x p1', '0', 'Color:', '0 0 0');
%     fprintf(fileID, '%12s\r\n%3s\r\n%4s %6s\r\n\r\n', 'Amb_mon_Resp', '0', 'Color:', alt_colors{1});
%     fprintf(fileID, '%17s\r\n%3s\r\n%4s %6s\r\n\r\n', 'Amb_mon_Resp x p1', '0', 'Color:', '0 0 0');
%     fprintf(fileID, '%16s\r\n%3s\r\n%4s %6s\r\n\r\n', 'Risk_mon_Display', '0', 'Color:', alt_colors{2});
%     fprintf(fileID, '%21s\r\n%3s\r\n%4s %6s\r\n\r\n', 'Risk_mon_Display x p1', '0', 'Color:', '0 0 0');
%     fprintf(fileID, '%13s\r\n%3s\r\n%4s %6s\r\n\r\n', 'Risk_mon_Resp', '0', 'Color:', alt_colors{2});
%     fprintf(fileID, '%18s\r\n%3s\r\n%4s %6s\r\n\r\n', 'Risk_mon_Resp x p1', '0', 'Color:', '0 0 0');
%   elseif ~is_mon & NumParametricWeights == 0
  if ~is_mon
    fprintf(fileID, '%13s\r\n%3s\r\n%4s %6s\r\n\r\n', 'mon_5_Display', '0', 'Color:', alt_colors{1}); %\r\n is starting a new line
    fprintf(fileID, '%10s\r\n%3s\r\n%4s %6s\r\n\r\n', 'mon_5_Resp', '0', 'Color:', alt_colors{1});
    fprintf(fileID, '%13s\r\n%3s\r\n%4s %6s\r\n\r\n', 'mon_8_Display', '0', 'Color:', alt_colors{2});
    fprintf(fileID, '%10s\r\n%3s\r\n%4s %6s\r\n\r\n', 'mon_8_Resp', '0', 'Color:', alt_colors{2});
    fprintf(fileID, '%14s\r\n%3s\r\n%4s %6s\r\n\r\n', 'mon_12_Display', '0', 'Color:', alt_colors{3});
    fprintf(fileID, '%11s\r\n%3s\r\n%4s %6s\r\n\r\n', 'mon_12_Resp', '0', 'Color:', alt_colors{3});
    fprintf(fileID, '%14s\r\n%3s\r\n%4s %6s\r\n\r\n', 'mon_25_Display', '0', 'Color:', alt_colors{4});
    fprintf(fileID, '%11s\r\n%3s\r\n%4s %6s\r\n\r\n', 'mon_25_Resp', '0', 'Color:', alt_colors{4});
  end

  % Print the $5 Display block for given domain
  fprintf(fileID, '%13s\r\n', [domain '_5' '_Display']); % even with parametric, 'x p1' is not included in the name
  fprintf(fileID, '%4s\r\n', num2str(size(block_5,1))); % Use size instead of length function, to preven when there is only 1 trial in for condition. Length funtion will give '2' when there is only 1 trial.
  if ~isempty(block_5)
   fprintf(fileID, '%4.0f\t %3.0f\r\n', block_5');
  end
  fprintf(fileID, '%6s %7s\r\n\r\n', 'Color:', Colors{1}); 
  % Print $5 block RESP for given domain
  fprintf(fileID, '%10s\r\n', [domain '_5' '_Resp']);
  fprintf(fileID, '%4s\r\n', num2str(size(block_5,1)));
  if ~isempty(block_5)  
    fprintf(fileID, '%4.0f\t %3.0f\r\n', block_5_resp');
  end
  fprintf(fileID, '%6s %7s\r\n\r\n', 'Color:', Colors{1});
  

  % Print the $8 Display block for given domain
  fprintf(fileID, '%13s\r\n', [domain '_8' '_Display']); 
  fprintf(fileID, '%4s\r\n', num2str(size(block_8,1)));
  if ~isempty(block_8)
    fprintf(fileID, '%4.0f\t %3.0f\r\n', block_8');
  end
  fprintf(fileID, '%6s %7s\r\n\r\n', 'Color:', Colors{2}); 
  % Print $8  block RESP for given domain
  fprintf(fileID, '%10s\r\n', [domain '_8' '_Resp']);
  fprintf(fileID, '%4s\r\n', num2str(size(block_8,1)));
  if ~isempty(block_8)
    fprintf(fileID, '%4.0f\t %3.0f\r\n', block_8_resp');
  end
  fprintf(fileID, '%6s %7s\r\n\r\n', 'Color:', Colors{2});
  

  % Print the $12 Display block for given domain
  fprintf(fileID, '%14s\r\n', [domain '_12' '_Display']); 
  fprintf(fileID, '%4s\r\n', num2str(size(block_12,1)));
  if ~isempty(block_12)
    fprintf(fileID, '%4.0f\t %3.0f\r\n', block_12');
  end
  fprintf(fileID, '%6s %7s\r\n\r\n', 'Color:', Colors{3}); 
  % Print $12 block RESP for given domain
  fprintf(fileID, '%11s\r\n', [domain '_12' '_Resp']);
  fprintf(fileID, '%4s\r\n', num2str(size(block_12,1)));
  if ~isempty(block_12)
    fprintf(fileID, '%4.0f\t %3.0f\r\n', block_12_resp');
  end
  fprintf(fileID, '%6s %7s\r\n\r\n', 'Color:', Colors{3});
  
  % Print the $25 Display block for given domain
  fprintf(fileID, '%14s\r\n', [domain '_25' '_Display']); 
  fprintf(fileID, '%4s\r\n', num2str(size(block_25,1)));
  if ~isempty(block_25)
    fprintf(fileID, '%4.0f\t %3.0f\r\n', block_25');
  end
  fprintf(fileID, '%6s %7s\r\n\r\n', 'Color:', Colors{4}); 
  % Print $25 block RESP for given domain
  fprintf(fileID, '%11s\r\n', [domain '_25' '_Resp']);
  fprintf(fileID, '%4s\r\n', num2str(size(block_25,1)));
  if ~isempty(block_25)
    fprintf(fileID, '%4.0f\t %3.0f\r\n', block_25_resp');
  end
  fprintf(fileID, '%6s %7s\r\n\r\n', 'Color:', Colors{4});
  
  % Equivalent empty-condition block from above
%   if is_mon & NumParametricWeights ~= 0
%     fprintf(fileID, '%15s\r\n%3s\r\n%4s %6s\r\n\r\n', 'Amb_med_Display', '0', 'Color:', alt_colors{1});
%     fprintf(fileID, '%20s\r\n%3s\r\n%4s %6s\r\n\r\n', 'Amb_med_Display x p1', '0', 'Color:', '0 0 0');
%     fprintf(fileID, '%12s\r\n%3s\r\n%4s %6s\r\n\r\n', 'Amb_med_Resp', '0', 'Color:', alt_colors{1});
%     fprintf(fileID, '%17s\r\n%3s\r\n%4s %6s\r\n\r\n', 'Amb_med_Resp x p1', '0', 'Color:', '0 0 0');
%     fprintf(fileID, '%16s\r\n%3s\r\n%4s %6s\r\n\r\n', 'Risk_med_Display', '0', 'Color:', alt_colors{2});
%     fprintf(fileID, '%21s\r\n%3s\r\n%4s %6s\r\n\r\n', 'Risk_med_Display x p1', '0', 'Color:', '0 0 0');
%     fprintf(fileID, '%13s\r\n%3s\r\n%4s %6s\r\n\r\n', 'Risk_med_Resp', '0', 'Color:', alt_colors{2});
%     fprintf(fileID, '%18s\r\n%3s\r\n%4s %6s\r\n\r\n', 'Risk_med_Resp x p1', '0', 'Color:', '0 0 0');
%   elseif is_mon & NumParametricWeights == 0
  if is_mon
    fprintf(fileID, '%13s\r\n%3s\r\n%4s %6s\r\n\r\n', 'med_5_Display', '0', 'Color:', alt_colors{1}); %\r\n is starting a new line
    fprintf(fileID, '%10s\r\n%3s\r\n%4s %6s\r\n\r\n', 'med_5_Resp', '0', 'Color:', alt_colors{1});
    fprintf(fileID, '%13s\r\n%3s\r\n%4s %6s\r\n\r\n', 'med_8_Display', '0', 'Color:', alt_colors{2});
    fprintf(fileID, '%10s\r\n%3s\r\n%4s %6s\r\n\r\n', 'med_8_Resp', '0', 'Color:', alt_colors{2});
    fprintf(fileID, '%14s\r\n%3s\r\n%4s %6s\r\n\r\n', 'med_12_Display', '0', 'Color:', alt_colors{3});
    fprintf(fileID, '%11s\r\n%3s\r\n%4s %6s\r\n\r\n', 'med_12_Resp', '0', 'Color:', alt_colors{3});
    fprintf(fileID, '%14s\r\n%3s\r\n%4s %6s\r\n\r\n', 'med_25_Display', '0', 'Color:', alt_colors{4});
    fprintf(fileID, '%11s\r\n%3s\r\n%4s %6s\r\n\r\n', 'med_25_Resp', '0', 'Color:', alt_colors{4});
  end
  
  fclose(fileID);
end
end
