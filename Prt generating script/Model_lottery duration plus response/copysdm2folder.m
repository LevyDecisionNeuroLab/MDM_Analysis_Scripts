%% Copy sdm files into folders of different parametric modulators
clear all
close all
%%
root = 'D:\Ruonan\Projects in the lab\MDM Project\Medical Decision Making Imaging\MDM_imaging\Behavioral Analysis';
sdmwave = 'SDM_Model_disp_resp_08222017';
path = [root filesep 'SDM files' filesep sdmwave filesep];
% cd (path)

% ParametricModType = {'type_none','RiskLevel','AmbiguityLevel','SV','RewardValue','SubjRating','Magnitude_none','Magnitude_NoRiskAmb_none'}; % subjective to change 
% ParametricModType = {'type_none','RiskLevel','AmbiguityLevel','SubjRating','RewardValue','SV','CV','CR','CRating'};
% ParametricModType = {'disp1TRresp_type_Magnitude_NoRiskAmb_none','disp1TRresp_type_Magnitude_NoRiskAmb_chosen_none','disp1TRresp_type_Magnitude_NoRiskAmb_lowhigh_none'};
% ParametricModType = {'Both_RiskLevel','Both_AmbiguityLevel'};
ParametricModType = {'MagnitudeByRiskOrAmb_none'};


for i = 1:length(ParametricModType)
    modulator = ParametricModType{i};
    
    desti = fullfile(path, modulator);
    source = fullfile(path, ['*_' modulator '.sdm']);
    movefile(source, desti)
end