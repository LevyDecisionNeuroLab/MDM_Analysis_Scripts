clear all
close all
%% Setup and settings
%% Input:
fitparwave = 'Behavior data fitpar_07172017\';
prtwave = 'Prt_Model_disp_resp_07172017\';

duration2use = 6; % How many tr in the display is modeled in the GLM

% Exclude subjects with ineligible imaging data
% subj 2062 has all subjective rating the same value. but the behavior is not the same
exclude = [2587 2590 2599];

%% Location of data files
root = 'D:\Ruonan\Projects in the lab\MDM Project\Medical Decision Making Imaging\MDM_imaging\Behavioral Analysis';
path_in = fullfile(root, 'Behavior fitpar files', fitparwave);

% Location to save PRT files in
path_out = fullfile(root, 'Prt files',prtwave);
% make folder if does not exist
if exist(path_out)==0
    mkdir(fullfile(root,'Prt files'),prtwave)
end

%% Computational parameters
tr = 1; % Temporal resolution, in seconds
trialduration = 6; % How many volumes *including onset* we analyze, in TR
DiscardedAcquisition = 15; % How many initial volumes we discard, it is the duration of the first trial. MDM_v1 and MDM_v2 are different.
tpb = 21; % trials per block

% Permissible values: 'RewardValue', 'RiskLevel', 'AmbiguityLevel', 'SV', or '' for no parameter
% NOTE: For more parameters, PTB_Protocol_Gen must be edited to (a) accept them, (b) calculate them
% ParametricModType = {'CR','CRating'};
ParametricModType = {'SubjRating_AL_RL'}; % one design matrix contains three parametric modulators : subjective rating of lottery, ambiguity level and risk level.

% NumParametricWeights is set by the script, depending on which ParametricModType is passed

% Get all subjects
% NOTE: Assuming that all subjects with MON files also have MED files
subj_files = dir([fullfile(path_in), filesep, 'MDM_MON*fitpar.mat']);
SubjectNums = zeros(1, length(subj_files));

% Extract subject from filename
for file_idx = 1:length(subj_files)
  fname = subj_files(file_idx).name;
  matches = regexp(fname, 'MDM_(?<domain>MON|MED)_(?<subjectNum>[\d]{1,4})', 'names');
  SubjectNums(file_idx) = str2num(matches.subjectNum); 
end
% SubjectNums = [2062]; % for single subject processing and testing

SubjectNums = SubjectNums(~ismember(SubjectNums, exclude));

domain = {'mon', 'med'}; % changed from 'gainloss'

% PRT file parameters
PRT.FileVersion =         '3'; % This has to be setup as '3' in order to add the parametricWeights
PRT.ResolutionOfTime =    'Volumes';
PRT.Experiment =          'RNA_MDM_FMRI';
PRT.NrOfConditions =      '8'; % NOTE: this is re-written for each ParametircModType in the loop below. 

PRT.BackgroundColor =     '255 255 255';
PRT.TextColor =           '0 0 0';
PRT.TimeCourseColor =     '0 0 0';
PRT.TimeCourseThick =     '2';
PRT.ReferenceFuncColor =  '30 200 30';
PRT.ReferenceFuncThick =  '2';

PRT.ColorMed =       '255 0 0 ';
PRT.ColorMed_SubjRating =      '140 80 80';
PRT.ColorMed_AmbLevel =      '255 40 40';
PRT.ColorMed_RiskLevel =     '77 0 0';
PRT.ColorMon =       '0 0 255';
PRT.ColorMon_SubjRating =      '80 80 140';
PRT.ColorMon_AmbLevel =      '40 40 255';
PRT.ColorMon_RiskLevel =     '0 0 140';
PRT.ColorResp = '50 50 50';

%% Run for all of the above
% Iterate for each subject, each domain (each _fitpar data file because loss and gains are separate)
for i = 1:length(SubjectNums)
    for j = 1:length(domain)
        for k = 1:length(ParametricModType)
            % for non parametric design matrics, only 3 condistions (2 lottery duration + 1 responses)
            if strcmp(ParametricModType{k}, '')
                PRT.NrOfConditions =      '3';
            else
                PRT.NrOfConditions =      '12';
            end
            PTB_MDM_Protocol_Gen_ver9(SubjectNums(i), domain{j}, ...
                    tr, trialduration, duration2use, DiscardedAcquisition, ...
                    ParametricModType{k}, ...
                    path_in, path_out, PRT,tpb)
        end
    end
end
