%%  This script adds column of zeros to the SDM files with parametric modulator of AmbiguityLevel or RiskLevel.
%   The reason is that, e.g,  in AmbiguityLevel prt files, parametric modulators are all 0 for risky trials, thus the two columns in sdm is missing 
%   After saving the new SDM, the number of RTCMatrix column will automatically match the new SDMMatrix
%   PredictorColors will automatically match the new SDMMatrix ONLY IF the new SDM file has more column than the former SDM files.
%   Need to change PredictorColors too.

clear all
close all
%% 
root = 'D:\Ruonan\Projects in the lab\MDM Project\Medical Decision Making Imaging\MDM_imaging\Behavioral Analysis';
sdmwave = 'SDM_Model_disp_resp_07172017';
path = [root filesep 'SDM files' filesep sdmwave filesep];

% make an AddingZeros directory if necessary.
if ~isdir(fullfile(path,'AddingZeros'))
    disp('Making AddingZeros directory');
    mkdir(path,'AddingZeros');
end

pathout = [root filesep 'SDM files' filesep sdmwave filesep 'AddingZeros' filesep];


% find the files to add zeros, _AmbiguityLevel and _RiskLevel, in a loop
ParametricModType = {'SubjRating_AL_RL'};

for i = 1:length(ParametricModType)
    modulator = ParametricModType{i};
    files2add = dir([path '*_' modulator '.sdm']); % change for _RiskLevel or _AmbiguityLevel

      for s = 1:length(files2add)

    %     for run = 1:8
    %         
    %         subject = 3;
    %         run = 1;
    %         
    %         session = fix(run/4)+1; %session (day) number for each run
    %         
    %         if ismember(run,[1 2 7 8]) %check if this run is gain or loss block
    %             gainloss = 'gains';
    %             block = mod(run, 4); % relative block number for each domain
    %         else
    %             gainloss = 'loss';
    %             block = run - 2;
    %         end
    %         
    %         sdm_name = ['Prt files' num2str(subject) '_S' num2str(session) '_block' num2str(run) '_' gainloss num2str(block) '_type_AmbiguityLevel.sdm']; % change for RiskLevel

            sdm_name = files2add(s).name
            sdm = xff([path sdm_name]);
            PredictorNames_in = sdm.PredictorNames; % same as PredictorNames_full, should be 16 in total

            PredictorNames_out = {'mon_Display', ...
                                  'mon_Display x SubjRating x p1', 'mon_Display x RiskLevel x p1', 'mon_Display x AmbLevel x p1', ...
                                  'med_Display',...
                                  'med_Display x SubjRating x p1', 'med_Display x RiskLevel x p1', 'med_Display x AmbLevel x p1' ...
                                  'Resp', 'Constant'}; % 10 in total
                              
            sdm.PredictorNames = PredictorNames_out;            
            sdm.NrOfPredictors = length(PredictorNames_out);
            sdm.FirstConfoundPredictor = sdm.NrOfPredictors;
                              
%             PredictorNames_full = {'mon_Display',...
%                                    'mon_Display x SubjRating', 'mon_Display x SubjRating x p1',...
%                                    'mon_Display x RiskLevel', 'mon_Display x RiskLevel x p1',...
%                                    'mon_Display x AmbLevel', 'mon_Display x AmbLevel x p1',...
%                                    'med_Display',...
%                                    'med_Display x SubjRating', 'med_Display x SubjRating x p1',...
%                                    'med_Display x RiskLevel', 'med_Display x RiskLevel x p1',...
%                                    'med_Display x AmbLevel', 'med_Display x AmbLevel x p1',...
%                                    'Resp', 'Constant'}; % 16



            % find the missing columns in _out relative to the _in, delete do not include them in _out
            member = 0; %How many columns in _in is also member of _out
            for m = 1:length(PredictorNames_in)
                exist = 0;
                for n = 1:length(PredictorNames_out)
                     if strcmp(PredictorNames_out(n), PredictorNames_in(m))
                        exist = 1;
                     end
                end
                if exist == 1; %the m-th column in PredictorNames_in is in PredictorNames_out
                    member = member+1;
                    include(member)= m; % store which columns of PredictorNames_in will be included
                end
            end

            %  sdm, include only the columns with exist in both _in and _out
            a = zeros(300,length(PredictorNames_out)); % volume number x condition numbers
            j=1; %column number for out SDM 
            for i=1:length(PredictorNames_in) %column number for in SDM
                if ismember(i,include)    %if the column is a match
                    a(:,j) = sdm.SDMMatrix(:,i);
                    j = j+1;
                end
            end
            
            sdm.SDMMatrix = a;


            %  Predictors color, include only the columns with exist in both _in and _out
            c = zeros(length(PredictorNames_out),3); % ondition numbers x rgb
            j=1; %column number for out SDM 
            for i=1:length(PredictorNames_in)
                if ismember(i,include)    %if the column is a match
                    c(j,:) = sdm.PredictorColors(i,:);
                    j = j+1;
                end
            end
            
            sdm.PredictorColors = c;

            
            % save it in a new folder
            sdm.SaveAs([pathout sdm_name]);
      end
end

% move the old sdm to the old folder
desti_old = fullfile(path,'sdm_missingColumns\');
source_old = fullfile(path, '*_SubjRating_AL_RL.sdm');
   
movefile(source_old,desti_old);

% move the new sdm out of the AddingZeros folder
source_new = fullfile(pathout, '*_SubjRating_AL_RL.sdm');

movefile(source_new,path);

rmdir(pathout) % remove AddingZeros folder
  
        
        
       