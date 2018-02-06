clear all
close all
%% Setup and settings
% Location of data files
root = 'D:\Ruonan\Projects in the lab\MDM Project\Medical Decision Making Imaging\MDM_imaging\Behavioral Analysis';
fitparwave = 'Behavior data fitpar_04182017\';
prtwave = 'Prt_Model_disp_resp_04182017\';

path_in = fullfile(root, 'Behavior fitpar files', fitparwave);

% Location to save PRT files in
path_out = fullfile(root, 'Prt files',prtwave);
% make folder if does not exist
if exist(path_out)==0
    mkdir(fullfile(root,'Prt files'),prtwave)
end


% Computational parameters
tr = 1; % Temporal resolution, in seconds
trialduration = 6; % How many volumes *including onset* we analyze, in volumes
DiscardedAcquisition = 15; % How many initial volumes we discard, it is the duration of the first trial. MDM_v1 and MDM_v2 are different.
tpb = 21; % trials per block

% Permissible values: 'RewardValue', 'RiskLevel', 'AmbiguityLevel', 'SV', or '' for no parameter
% NOTE: For more parameters, PTB_Protocol_Gen must be edited to (a) accept them, (b) calculate them
% ParametricModType = {'SV', 'RewardValue', 'RiskLevel', 'AmbiguityLevel', ''};
ParametricModType = {'','RiskLevel','AmbiguityLevel','SubjRating','RewardValue','SV'};

% NumParametricWeights is set by the script, depending on which ParametricModType is passed

% Get all subjects
% NOTE: Assuming that all subjects with MON files also have MED files
subj_files = dir([fullfile(path_in), filesep, 'MDM_MON*fitpar.mat']);
SubjectNums = zeros(1, length(subj_files));
% SubjectNums = [2581]; % for single subject processing and testing

% Extract subject from filename
for file_idx = 1:length(subj_files)
  fname = subj_files(file_idx).name;
  matches = regexp(fname, 'MDM_(?<domain>MON|MED)_(?<subjectNum>[\d]{1,4})', 'names');
  SubjectNums(file_idx) = str2num(matches.subjectNum); 
end

% Exclude subjects with ineligible imaging data
exclude = [2587 2590];
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

PRT.ColorAmb_Med =       '0 0 77';
PRT.ColorRisk_Med =      '140 0 0';
PRT.ColorAmb_Mon =      '0 0 255';
PRT.ColorRisk_Mon =     '255 0 0';
PRT.ColorAmb_Med_Resp =       '40 40 77';
PRT.ColorRisk_Med_Resp =      '140 80 80';
PRT.ColorAmb_Mon_Resp =      '80 80 255';
PRT.ColorRisk_Mon_Resp =     '255 40 40';

%% Run for all of the above
% Iterate for each subject, each domain (each _fitpar data file because loss and gains are separate)
for i = 1:length(SubjectNums)
    for j = 1:length(domain)
        for k = 1:length(ParametricModType)
            % for non parametric design matrics, only 8 condistions (4 lottery duration + 4 responses)
            if strcmp(ParametricModType{k}, '')
                PRT.NrOfConditions =      '8';
            else
            % for parametic design matrics, 8+4 conditions (add a parametric modulator for each empty predictor, the other four will  be added when SDM is created)
                PRT.NrOfConditions =      '12';
            end
            if (strcmp(ParametricModType{k},'SV')|strcmp(ParametricModType{k},'RewardValue')) & strcmp(domain{j},'med') % for medical blocks, these two parameters are meanlingless
            else
                PTB_MDM_Protocol_Gen_ver3(SubjectNums(i), domain{j}, ...
                    tr, trialduration, DiscardedAcquisition, ...
                    ParametricModType{k}, ...
                    path_in, path_out, PRT,tpb)
            end
        end
    end
end
