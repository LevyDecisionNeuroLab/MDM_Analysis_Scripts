%%  This script adds column of zeros to the SDM files with parametric modulator of AmbiguityLevel or RiskLevel.
%   The reason is that, e.g,  in AmbiguityLevel prt files, parametric modulators are all 0 for risky trials, thus the two columns in sdm is missing 
%   After saving the new SDM, the number of RTCMatrix column and the number
%   of predictor colors will automatically match the new SDMMatrix.

clear all
close all
%% 
root = 'D:\Ruonan\Projects in the lab\MDM Project\Medical Decision Making Imaging\MDM_imaging\Behavioral Analysis';
sdmwave = 'SDM_Model_disp_resp_08142017';
path = [root filesep 'SDM files' filesep sdmwave filesep];

% make an AddingZeros directory if necessary.
if ~isdir(fullfile(path,'AddingZeros'))
    disp('Making AddingZeros directory');
    mkdir(path,'AddingZeros');
end

pathout = [root filesep 'SDM files' filesep sdmwave filesep 'AddingZeros' filesep];


% find the files to add zeros, _AmbiguityLevel and _RiskLevel, in a loop
% % ParametricModType = {'RiskLevel','AmbiguityLevel'};
ParametricModType = {'Both_RiskLevel', 'Both_AmbiguityLevel'};

for i = 1:length(ParametricModType)
    modulator = ParametricModType{i};
    files2add = dir([path '*_' modulator '.sdm']); % change for _RiskLevel or _AmbiguityLevel
    
      for s = 1:length(files2add)
            sdm_name = files2add(s).name
            domain = sdm_name(13:15);
            sdm = xff([path sdm_name]);
            sdm.NrOfPredictors = 10;
            sdm.FirstConfoundPredictor = 10;
            PredictorNames9 = sdm.PredictorNames;
            PredictorNames10 = {'Amb_mon_Display', 'Amb_mon_Display x p1',  ...
                                'Risk_mon_Display', 'Risk_mon_Display x p1', ...
                                'Amb_med_Display', 'Amb_med_Display x p1',  ...
                                'Risk_med_Display', 'Risk_med_Display x p1', ...
                                'Resp', 'Constant'};
            sdm.PredictorNames = PredictorNames10;

            % find the missing column
            miss = 0; %How many columns are missing
            for m = 1:length(PredictorNames10)
                exist = 0;
                for n = 1:length(PredictorNames9)
                     if strcmp(PredictorNames9(n), PredictorNames10(m))
                        exist = exist+1; %the m-th column in PredictorNames10 is found in PredictorNames9
                     end
                end
                if exist == 0; %the m-th column in PredictorNames10 is missing in PredictorNames9
                    miss = miss+1;
                    mismatch(miss)= m; % store which column of PredictorNames10 is missing
                end
            end

            % insert the missing column and fill it with zeros
            a = zeros(300,sdm.NrOfPredictors); % volume number x condition numbers
            i=1; %column number for new SDM (length = 13)
            j=1; %column number for old SDM (length = 12)
            while i < sdm.NrOfPredictors+1
                if isempty(find(mismatch(mismatch == i)))== 1    %if the column is not missing
                    a(:,i) = sdm.SDMMatrix(:,j);
                    i = i+1;
                    j = j+1;
                else
                    i=i+1;
                end
            end
            
            if strcmp(modulator,'RiskLevel')|strcmp(modulator,'AmbiguityLevel') 
                sdm.SDMMatrix = a;
            elseif strcmp(modulator, 'Both_RiskLevel')|strcmp(modulator,'Both_AmbiguityLevel')
                sdm.NrOfPredictors = 8;
                sdm.FirstConfoundPredictor = 8;
                sdm.PredictorNames = {'Amb_mon_Display', 'Risk_mon_Display',  ...
                                'Amb_med_Display', 'Risk_med_Display',  ...
                                'Amb_Display x p1', 'Risk_Display x p1', ...
                                'Resp', 'Constant'};
                b = zeros(300,8);
                b(:,1:4) = a(:,[1,3,5,7]); % display predictors
                b(:,7:8) = a(:,9:10); % response and constant predictors
                if strcmp(domain,'mon')
                    b(:,5:6) = a(:,[2,4]); % parametric modulators for monetary
                elseif strcmp(domain,'med')
                    b(:,5:6) = a(:,[6,8]); % parametric modulators for medical
                end
                sdm.SDMMatrix = b;
                sdm.PredictorColors = sdm.PredictorColors(1:8,:);
            end
            % save it in original file name
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
  
        
        
       