function PTB_MDM_Protocol_Gen_ver8(subjectNum, domain, tr, trialduration, duration2use, DiscAcq, ParametricMod, path_in, path_out, PRT, tpb)
%PTB_PROTOCOL_GEN Generates PRT Files from Psychtoolbox Data Files
% This version does not create paramatric modulators for response predictor

%INPUTS:
%       subjectNum - Specifies subject file to be loaded
%       domain - A string of either 'monetary' or 'medical' to indicate the domain to extract
%       tr - Temporal resolution of the visual data, in seconds
%       trialduration - TRs fo trial display in the experimental design
%       duration2use - TRs of trial display to model in the GLM
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
% Both_RiskLevel-Risk Level not separated by domain, Both_AmbiguityLevel-Ambiguity level not separated by domain
% CV-Chosen Subjective value, CR-Chosen Reward magnitude, CRating-chosen subjective rating
permissible_parameters = {'RewardValue', 'RiskLevel', 'AmbiguityLevel', 'Both_RiskLevel', 'Both_AmbiguityLevel', 'SubjRating','SV','CV','CR','CRating'}; 


if sum(strcmp(ParametricMod, permissible_parameters))
  % using sum because strcmp returns a matrix of logicals if string is compared 
  % with an array of strings
  NumParametricWeights = 1;
elseif strcmp(ParametricMod, 'none')
  NumParametricWeights = 0;
else
  error('Aborting: `ParametricMod` must be set to a valid parameter value!')
  return
end
PRT.ParametricWeights = num2str(NumParametricWeights);

% For PRT properties, fix colors into correct order
if ~is_mon
  Colors = {PRT.ColorAmb_Med, PRT.ColorRisk_Med};
  alt_colors = {PRT.ColorAmb_Mon, PRT.ColorRisk_Mon};
else
  Colors = {PRT.ColorAmb_Mon, PRT.ColorRisk_Mon};
  alt_colors = {PRT.ColorAmb_Med, PRT.ColorRisk_Med};
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

%% Compute subjective value of each choice only for monetary
% the fir_parameters_mdm.m is updated so that sv is already calculated, can
% decide whether to use constrained or unconstrained

% ToDo: If logic, to use raw or rating value 

if is_mon
    sv = data.svUncstr; % subjective value of lottery
    cv = data.chosenSV; % subjective value of chosen side
end


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
  [ onsets, offsets ] = getBlockOnsetOffset(prtblock, DiscAcq, tr, trialduration, duration2use, tpb);

  %% Compute values from current block
  % Get indices for risk trials and ambiguity trials (exclude first trial in each block)
  current_block_range = (2:tpb) + (blocknum - 1) * tpb;

  % Divide blocks by type of lottery (store indices with risk-only vs. risk-and-ambiguity)
  % These _indexes are in terms of the first trial deleted!!
  amb_index = data.ambigs(current_block_range) > 0;
  risk_index = data.ambigs(current_block_range) == 0;
 

  % Parametric weights for current block
  if NumParametricWeights > 0 & is_mon
    block_amt = data.vals(current_block_range);
    block_rlevel = data.probs(current_block_range);
    block_alevel = data.ambigs(current_block_range);
    block_rtg = data.subjRatings(current_block_range); % subjective rating of lottery
    block_crating = data.chosenRating(current_block_range); % subjective rating of chosen one
    block_sv = sv(current_block_range); %subjective value of lottery
    block_cv = cv(current_block_range); % subjective value of the chosen one
    block_cr = data.chosenVal(current_block_range); % reward magnitude of chosen one
  elseif NumParametricWeights > 0 & is_mon == 0
    block_amt = data.vals(current_block_range);
    block_rlevel = data.probs(current_block_range);
    block_alevel = data.ambigs(current_block_range);
    block_rtg = data.subjRatings(current_block_range); % subjective rating of lottery
    block_crating = data.chosenRating(current_block_range); % subjective rating of chosen one
  end

  %% Select what parametric value to write (or not write) into PRT file
  % Store the basic onset/offset, computed earlier
  block_amb = [onsets(amb_index,1) offsets(amb_index,1)];
  block_risk = [onsets(risk_index,1) offsets(risk_index,1)];
  resp = [onsets(:,2) offsets(:,2)]; % response for all trials
%   block_amb_resp = [onsets(amb_index,2) offsets(amb_index,2)];
%   block_risk_resp = [onsets(risk_index,2) offsets(risk_index,2)];

  % Add the selected parametric value if required
  if NumParametricWeights > 0
    if strcmp(ParametricMod, 'RewardValue')& is_mon
      block_amb = [block_amb block_amt(amb_index)]; % onset-offset-paramatricValue
      block_risk = [block_risk block_amt(risk_index)];
    elseif strcmp(ParametricMod, 'RiskLevel')|strcmp(ParametricMod,'Both_RiskLevel')
      block_amb = [block_amb block_rlevel(amb_index)];
      block_risk = [block_risk block_rlevel(risk_index)];
    elseif strcmp(ParametricMod, 'AmbiguityLevel')|strcmp(ParametricMod,'Both_AmbiguityLevel')
      block_amb = [block_amb block_alevel(amb_index)];
      block_risk = [block_risk block_alevel(risk_index)];
    elseif strcmp(ParametricMod, 'SubjRating')
      block_amb = [block_amb block_rtg(amb_index)];
      block_risk = [block_risk block_rtg(risk_index)];
    elseif strcmp(ParametricMod, 'SV') & is_mon
      block_amb = [block_amb block_sv(amb_index)];
      block_risk = [block_risk block_sv(risk_index)];
    elseif strcmp(ParametricMod, 'CV') & is_mon
      block_amb = [block_amb block_cv(amb_index)];
      block_risk = [block_risk block_cv(risk_index)];
    elseif strcmp(ParametricMod, 'CR') & is_mon
      block_amb = [block_amb block_cr(amb_index)];
      block_risk = [block_risk block_cr(risk_index)];
    elseif strcmp(ParametricMod, 'CRating')
      block_amb = [block_amb block_crating(amb_index)];
      block_risk = [block_risk block_crating(risk_index)];
    end
  end

  %% Write file to txt
  % TODO: Called once per block - abstract into a subfunction / nested function?

  % Open file for writing
%   if strcmp(ParametricMod, 'none')
%     ParametricMod = 'none';
%   end
  fname = [path_out num2str(subjectNum) '_block' block_order{blocknum} ...
      '_' domain num2str(blocknum) '_model_disp' num2str(duration2use) 'TR' 'resp' '_type_' ParametricMod '.prt'];
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

  % TODO: Shouldn't ParametricWeights be always included, even if they *are* 0?
  if NumParametricWeights > 0
      fprintf(fileID, '%18s %4s\r\n', 'ParametricWeights:', PRT.ParametricWeights);
  end

  fprintf(fileID, '\r\n%15s %7s\r\n\r\n', 'NrOfConditions:', PRT.NrOfConditions);

  % NOTE: Empty conditions -- set here for purposes of consistent order. 
  % Equivalent but opposite block further below
  if ~is_mon & NumParametricWeights ~= 0 % medical block with parametric
    fprintf(fileID, '%15s\r\n%3s\r\n%4s %6s\r\n\r\n', 'Amb_mon_Display', '0', 'Color:', alt_colors{1});
    fprintf(fileID, '%20s\r\n%3s\r\n%4s %6s\r\n\r\n', 'Amb_mon_Display x p1', '0', 'Color:', '0 0 0');
    fprintf(fileID, '%16s\r\n%3s\r\n%4s %6s\r\n\r\n', 'Risk_mon_Display', '0', 'Color:', alt_colors{2});
    fprintf(fileID, '%21s\r\n%3s\r\n%4s %6s\r\n\r\n', 'Risk_mon_Display x p1', '0', 'Color:', '0 0 0');
  elseif ~is_mon & NumParametricWeights == 0
    fprintf(fileID, '%15s\r\n%3s\r\n%4s %6s\r\n\r\n', 'Amb_mon_Display', '0', 'Color:', alt_colors{1});
    fprintf(fileID, '%16s\r\n%3s\r\n%4s %6s\r\n\r\n', 'Risk_mon_Display', '0', 'Color:', alt_colors{2});
  end

  % Print the ambiguity Display block for given domain
  % The other domain will be added after SDM is created
  fprintf(fileID, '%15s\r\n', ['Amb_' domain '_Display']); % even with parametric, 'x p1' is not included in the name
  fprintf(fileID, '%4s\r\n', num2str(length(block_amb)));

  if NumParametricWeights == 0
      fprintf(fileID, '%4.0f\t %3.0f\r\n', block_amb');
  else
      fprintf(fileID, '%4.0f\t %3.0f\t %1.3f\r\n', block_amb');
  end

  fprintf(fileID, '%6s %7s\r\n\r\n', 'Color:', Colors{1});

  % Print the risk block Display for given domain
  fprintf(fileID, '%16s\r\n', ['Risk_' domain '_Display']);
  fprintf(fileID, '%4s\r\n', num2str(length(block_risk)));

  if NumParametricWeights == 0
    fprintf(fileID, '%4.0f\t %3.0f\r\n', block_risk');
  else
    fprintf(fileID, '%4.0f\t %3.0f\t %1.3f\r\n', block_risk');
  end

  fprintf(fileID, '%6s %7s\r\n\r\n', 'Color:', Colors{2});
  
  % Equivalent empty-condition block from above
  if is_mon & NumParametricWeights ~= 0
    fprintf(fileID, '%15s\r\n%3s\r\n%4s %6s\r\n\r\n', 'Amb_med_Display', '0', 'Color:', alt_colors{1});
    fprintf(fileID, '%20s\r\n%3s\r\n%4s %6s\r\n\r\n', 'Amb_med_Display x p1', '0', 'Color:', '0 0 0');
    fprintf(fileID, '%16s\r\n%3s\r\n%4s %6s\r\n\r\n', 'Risk_med_Display', '0', 'Color:', alt_colors{2});
    fprintf(fileID, '%21s\r\n%3s\r\n%4s %6s\r\n\r\n', 'Risk_med_Display x p1', '0', 'Color:', '0 0 0');
  elseif is_mon & NumParametricWeights == 0
    fprintf(fileID, '%15s\r\n%3s\r\n%4s %6s\r\n\r\n', 'Amb_med_Display', '0', 'Color:', alt_colors{1});
    fprintf(fileID, '%16s\r\n%3s\r\n%4s %6s\r\n\r\n', 'Risk_med_Display', '0', 'Color:', alt_colors{2});
  end
  
  % Print the response for all trials
  fprintf(fileID, '%4s\r\n', 'Resp');
  fprintf(fileID, '%4s\r\n', num2str(length(resp)));
  fprintf(fileID, '%4.0f\t %3.0f\r\n', resp');
  fprintf(fileID, '%6s %7s\r\n\r\n', 'Color:', Colors{2});
  
  fclose(fileID);
end
end
