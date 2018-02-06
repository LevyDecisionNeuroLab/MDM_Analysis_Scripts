%%  This script adds column of zeros to the SDM files with parametric modulator of AmbiguityLevel or RiskLevel.
%   The reason is that, e.g,  in AmbiguityLevel prt files, parametric modulators are all 0 for risky trials, thus the two columns in sdm is missing 
%   After saving the new SDM, the number of RTCMatrix column and the number
%   of predictor colors will automatically match the new SDMMatrix.

clear all
close all
%% 
root = 'D:\Ruonan\Projects in the lab\MDM Project\Medical Decision Making Imaging\MDM_imaging\Behavioral Analysis';
sdmwave = 'SDM_Model_disp_resp_04242017';
path = [root filesep 'SDM files' filesep sdmwave filesep];

% make an AddingZeros directory if necessary.
if ~isdir(fullfile(path,'AddingZeros'))
    disp('Making AddingZeros directory');
    mkdir(path,'AddingZeros');
end

pathout = [root filesep 'SDM files' filesep sdmwave filesep 'AddingZeros' filesep];


% find the files to add zeros, _AmbiguityLevel and _RiskLevel, in a loop
ParametricModType = {'RiskLevel','AmbiguityLevel'};

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
            sdm.NrOfPredictors = 13;
            sdm.FirstConfoundPredictor = 13;
            PredictorNames12 = sdm.PredictorNames;
            PredictorNames13 = {'Amb_mon_Display', 'Amb_mon_Display x p1', 'Amb_mon_Resp', ...
                                'Risk_mon_Display', 'Risk_mon_Display x p1', 'Risk_mon_Resp', ...
                                'Amb_med_Display', 'Amb_med_Display x p1', 'Amb_med_Resp', ...
                                'Risk_med_Display', 'Risk_med_Display x p1', 'Risk_med_Resp', ...
                                'Constant'};
            sdm.PredictorNames = PredictorNames13;

            % find the missing column
            miss = 0; %How many columns are missing
            for m = 1:length(PredictorNames13)
                exist = 0;
                for n = 1:length(PredictorNames12)
                     if strcmp(PredictorNames12(n), PredictorNames13(m))
                        exist = exist+1; %the m-th column in PredictorNames13 is found in PredictorNames12
                     end
                end
                if exist == 0; %the m-th column in PredictorNames13 is missing in PredictorNames12
                    miss = miss+1;
                    mismatch(miss)= m; % store which column of PredictorNames13 is missing
                end
            end

            % insert the missing column and fill it with zeros
            a = zeros(300,sdm.NrOfPredictors); % volume number x condition numbers
            i=1; %column number for new SDM (length = 13)
            j=1; %column number for old SDM (length = 12)
            while i < sdm.NrOfPredictors+1
                if isempty(find(mismatch(mismatch == i)))== 1    %if the column is not missting
                    a(:,i) = sdm.SDMMatrix(:,j);
                    i = i+1;
                    j = j+1;
                else
                    i=i+1;
                end
            end
            sdm.SDMMatrix = a;

            % save it in a new file name
            sdm.SaveAs([pathout sdm_name]);
      end
end

% move the missing-column  _RiskLevel and _Ambiguitylevel sdm to the old folder
desti_old = fullfile(path,'sdm_missingColumns\');
source_risk = fullfile(path, '*_RiskLevel.sdm');
source_ambig = fullfile(path, '*_AmbiguityLevel.sdm');
   
movefile(source_risk,desti_old);
movefile(source_ambig,desti_old);

% move the adding-column sdm out of the AddingZeros folder
source_risk = fullfile(pathout, '*_RiskLevel.sdm');
source_ambig = fullfile(pathout, '*_AmbiguityLevel.sdm');

movefile(source_risk,path);
movefile(source_ambig,path);

rmdir(pathout) % remove AddingZeros folder
  
        
        
       