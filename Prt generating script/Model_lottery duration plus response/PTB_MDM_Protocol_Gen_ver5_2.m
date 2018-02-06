function PTB_MDM_Protocol_Gen_ver5_2(subjectNum, domain, tr, trialduration, duration2use, DiscAcq, ParametricMod, path_in, path_out, PRT, tpb)
%PTB_PROTOCOL_GEN Generates PRT Files from Psychtoolbox Data Files

% TODo: 
% Conditionally use raw value or rating value for RewardValue
% How to calculate Subjective value 

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
permissible_parameters = {'RewardValue', 'RiskLevel', 'AmbiguityLevel', 'SubjRating','SV'};
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
  Colors = {PRT.ColorMed_r25, PRT.ColorMed_r50, PRT.ColorMed_r75, PRT.ColorMed_a24, PRT.ColorMed_a50, PRT.ColorMed_a74, PRT.ColorMed_resp};
  alt_colors = {PRT.ColorMon_r25, PRT.ColorMon_r50, PRT.ColorMon_r75, PRT.ColorMon_a24, PRT.ColorMon_a50, PRT.ColorMon_a74, PRT.ColorMon_resp};
else
  Colors = {PRT.ColorMon_r25, PRT.ColorMon_r50, PRT.ColorMon_r75, PRT.ColorMon_a24, PRT.ColorMon_a50, PRT.ColorMon_a74, PRT.ColorMon_resp};
  alt_colors = {PRT.ColorMed_r25, PRT.ColorMed_r50, PRT.ColorMed_r75, PRT.ColorMed_a24, PRT.ColorMed_a50, PRT.ColorMed_a74, PRT.ColorMed_resp};
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
  [ onsets, offsets ] = getBlockOnsetOffset(prtblock, DiscAcq, tr, trialduration, duration2use, tpb);

  %% Compute indexes from current block
  % Get indices for risk trials and ambiguity trials (exclude first trial in each block)
  current_block_range = (2:tpb) + (blocknum - 1) * tpb;

  % Divide blocks by type of lottery (4 outcome magnitudes)
  
  % These _indexes are in terms of the first trial deleted!!
  
  % index for risk levels
  index_r25 = data.probs(current_block_range) == 0.25;
  index_r50 = data.probs(current_block_range) == 0.5 & data.ambigs(current_block_range) == 0;
  index_r75 = data.probs(current_block_range) == 0.75;
   
  % index for ambiguity levels
  index_a24 = data.ambigs(current_block_range) == 0.24;
  index_a50 = data.ambigs(current_block_range) == 0.5;
  index_a74 = data.ambigs(current_block_range) == 0.74;
  

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
  % if lottery value or chosen value
  if strcmp(ParametricMod, 'none')
      block_r25 = [onsets(index_r25,1) offsets(index_r25,1)];
      block_r50 = [onsets(index_r50,1) offsets(index_r50,1)];
      block_r75 = [onsets(index_r75,1) offsets(index_r75,1)];
      block_a24 = [onsets(index_a24,1) offsets(index_a24,1)];
      block_a50 = [onsets(index_a50,1) offsets(index_a50,1)];
      block_a74 = [onsets(index_a74,1) offsets(index_a74,1)];
  end
  
  resp = [onsets(:,2) offsets(:,2)];

  % Add the selected parametric value if required
%   if NumParametricWeights > 0
%     if strcmp(ParametricMod, 'RewardValue')& is_mon
%       block_amb = [block_amb block_amt(amb_index)]; % onset-offset-paramatricValue
%       block_risk = [block_risk block_amt(risk_index)];
%     elseif strcmp(ParametricMod, 'RiskLevel')
%       block_amb = [block_amb block_rlevel(amb_index)];
%       block_risk = [block_risk block_rlevel(risk_index)];
%     elseif strcmp(ParametricMod, 'AmbiguityLevel')
%       block_amb = [block_amb block_alevel(amb_index)];
%       block_risk = [block_risk block_alevel(risk_index)];
%     elseif strcmp(ParametricMod, 'SubjRating')
%       block_amb = [block_amb block_rtg(amb_index)];
%       block_risk = [block_risk block_rtg(risk_index)];
%     elseif strcmp(ParametricMod, 'SV') & is_mon
%       block_amb = [block_amb block_sv(amb_index)];
%       block_risk = [block_risk block_sv(risk_index)];
%     end
%   end

  %% Write file to txt
  % TODO: Called once per block - abstract into a subfunction / nested function?

  % Open file for writing
%   if strcmp(ParametricMod, '')
%     ParametricMod = 'none';
%   end
  fname = [path_out num2str(subjectNum) '_block' block_order{blocknum} ...
      '_' domain num2str(blocknum) '_model_disp' num2str(duration2use) 'TR' 'resp' '_type_' 'RiskAmbLevel_' ParametricMod '.prt'];
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
        fprintf(fileID, '%15s\r\n%3s\r\n%4s %6s\r\n\r\n', 'mon_r25_Display', '0', 'Color:', alt_colors{1}); %\r\n is starting a new line
        fprintf(fileID, '%15s\r\n%3s\r\n%4s %6s\r\n\r\n', 'mon_r50_Display', '0', 'Color:', alt_colors{2});
        fprintf(fileID, '%15s\r\n%3s\r\n%4s %6s\r\n\r\n', 'mon_r75_Display', '0', 'Color:', alt_colors{3});
        fprintf(fileID, '%15s\r\n%3s\r\n%4s %6s\r\n\r\n', 'mon_a24_Display', '0', 'Color:', alt_colors{4});
        fprintf(fileID, '%15s\r\n%3s\r\n%4s %6s\r\n\r\n', 'mon_a50_Display', '0', 'Color:', alt_colors{5});
        fprintf(fileID, '%15s\r\n%3s\r\n%4s %6s\r\n\r\n', 'mon_a74_Display', '0', 'Color:', alt_colors{6});        
      end

      % Print the r25 Display block for given domain
      fprintf(fileID, '%15s\r\n', [domain '_r25' '_Display']); % even with parametric, 'x p1' is not included in the name
      fprintf(fileID, '%4s\r\n', num2str(size(block_r25,1))); % Use size instead of length function, to preven when there is only 1 trial in for condition. Length funtion will give '2' when there is only 1 trial.
      if ~isempty(block_r25)
       fprintf(fileID, '%4.0f\t %3.0f\r\n', block_r25');
      end
      fprintf(fileID, '%6s %7s\r\n\r\n', 'Color:', Colors{1});   

      % Print the r50 Display block for given domain
      fprintf(fileID, '%15s\r\n', [domain '_r50' '_Display']); 
      fprintf(fileID, '%4s\r\n', num2str(size(block_r50,1)));
      if ~isempty(block_r50)
        fprintf(fileID, '%4.0f\t %3.0f\r\n', block_r50');
      end
      fprintf(fileID, '%6s %7s\r\n\r\n', 'Color:', Colors{2}); 


      % Print the r75 Display block for given domain
      fprintf(fileID, '%15s\r\n', [domain '_r75' '_Display']); 
      fprintf(fileID, '%4s\r\n', num2str(size(block_r75,1)));
      if ~isempty(block_r75)
        fprintf(fileID, '%4.0f\t %3.0f\r\n', block_r75');
      end
      fprintf(fileID, '%6s %7s\r\n\r\n', 'Color:', Colors{3}); 

      % Print the a24 Display block for given domain
      fprintf(fileID, '%15s\r\n', [domain '_a24' '_Display']); 
      fprintf(fileID, '%4s\r\n', num2str(size(block_a24,1)));
      if ~isempty(block_a24)
        fprintf(fileID, '%4.0f\t %3.0f\r\n', block_a24');
      end
      fprintf(fileID, '%6s %7s\r\n\r\n', 'Color:', Colors{4}); 

      % Print the a50 Display block for given domain
      fprintf(fileID, '%15s\r\n', [domain '_a50' '_Display']); 
      fprintf(fileID, '%4s\r\n', num2str(size(block_a50,1)));
      if ~isempty(block_a50)
        fprintf(fileID, '%4.0f\t %3.0f\r\n', block_a50');
      end
      fprintf(fileID, '%6s %7s\r\n\r\n', 'Color:', Colors{5}); 

       % Print the a74 Display block for given domain
      fprintf(fileID, '%15s\r\n', [domain '_a74' '_Display']); 
      fprintf(fileID, '%4s\r\n', num2str(size(block_a74,1)));
      if ~isempty(block_a74)
        fprintf(fileID, '%4.0f\t %3.0f\r\n', block_a74');
      end
      fprintf(fileID, '%6s %7s\r\n\r\n', 'Color:', Colors{6}); 


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
        fprintf(fileID, '%15s\r\n%3s\r\n%4s %6s\r\n\r\n', 'med_r25_Display', '0', 'Color:', alt_colors{1}); %\r\n is starting a new line
        fprintf(fileID, '%15s\r\n%3s\r\n%4s %6s\r\n\r\n', 'med_r50_Display', '0', 'Color:', alt_colors{2});
        fprintf(fileID, '%15s\r\n%3s\r\n%4s %6s\r\n\r\n', 'med_r75_Display', '0', 'Color:', alt_colors{3});
        fprintf(fileID, '%15s\r\n%3s\r\n%4s %6s\r\n\r\n', 'med_a24_Display', '0', 'Color:', alt_colors{4});
        fprintf(fileID, '%15s\r\n%3s\r\n%4s %6s\r\n\r\n', 'med_a50_Display', '0', 'Color:', alt_colors{5});
        fprintf(fileID, '%15s\r\n%3s\r\n%4s %6s\r\n\r\n', 'med_a74_Display', '0', 'Color:', alt_colors{6});
      end

  
  % Print the response for all trials
  fprintf(fileID, '%4s\r\n', 'Resp');
  fprintf(fileID, '%4s\r\n', num2str(length(resp)));
  fprintf(fileID, '%4.0f\t %3.0f\r\n', resp');
  fprintf(fileID, '%6s %7s\r\n\r\n', 'Color:', Colors{7});
  
  fclose(fileID);
end
end
