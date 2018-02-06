clear all
close all

%% This version of protocol generating scripts create prt files for seperate binary predictors for each reward magnitude level. 
% 4 seperate predictors, for $5, $8, $12, $25
% 4 seperate predictors, for slight, moderate, major, recovery
% seperate for mon and med
% DO NOT seperate for risk and ambig
% display + resposne
% so it should be 4magnitude*2display/response=8 predicotrs for each domain

%% Setup and settings
%% Changing Input
fitparwave = 'Behavior data fitpar_07172017\';
prtwave = 'Prt_Model_disp_resp_07172017\';

duration2use = 6;

%% Location of data files
root = 'D:\Ruonan\Projects in the lab\MDM Project\Medical Decision Making Imaging\MDM_imaging\Behavioral Analysis';

path_in = fullfile(root, 'Behavior fitpar files', fitparwave);

% Location to save PRT files in
path_out = fullfile(root, 'Prt files',prtwave);
% make folder if does not exist
if exist(path_out)==0
    mkdir(fullfile(root,'Prt files\'),prtwave)
end

%% Computational parameters
tr = 1; % Temporal resolution, in seconds
trialduration = 6; % How many volumes *including onset* we analyze
DiscardedAcquisition = 15; % How many initial volumes we discard, it is the duration of the first trial. MDM_v1 and MDM_v2 are different.
tpb = 21; % trials per block

% Permissible values: 'RewardValue', 'RiskLevel', 'AmbiguityLevel', 'SV', or '' for no parameter
% NOTE: For more parameters, PTB_Protocol_Gen must be edited to (a) accept them, (b) calculate them
% ParametricModType = {'SV', 'RewardValue', 'RiskLevel', 'AmbiguityLevel', ''};

% none - grouped by lottery value, chosen_none - gourped by chosen value, lowhigh_none - 5 and 8 as low, 12 and 25 as high
ParametricModType = {'none'};


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
exclude = [2587 2590 2599];
SubjectNums = SubjectNums(~ismember(SubjectNums, exclude));

domain = {'mon', 'med'}; % changed from 'gainloss'

% PRT file parameters
PRT.FileVersion =         '3'; % This has to be setup as '3' in order to add the parametricWeights
PRT.ResolutionOfTime =    'Volumes';
PRT.Experiment =          'RNA_MDM_FMRI';
PRT.NrOfConditions =      '16'; % NOTE: this is re-written for each ParametircModType in the loop below. 

PRT.BackgroundColor =     '255 255 255'; % white
PRT.TextColor =           '0 0 0'; % black
PRT.TimeCourseColor =     '0 0 0';
PRT.TimeCourseThick =     '2';
PRT.ReferenceFuncColor =  '30 200 30'; % Green
PRT.ReferenceFuncThick =  '2';

% for ifat's grant graph
PRT.ColorMed_r25 =    '254 0 3';
PRT.ColorMed_r50 =    '0 43 209';
PRT.ColorMed_r75 =   '105 255 42';
% PRT.ColorMed_r25 =    '255 158 158';
% PRT.ColorMed_r50 =    '231 80 80';
% PRT.ColorMed_r75 =   '168 0 0';
PRT.ColorMed_a24 =   '158 158 255';
PRT.ColorMed_a50 =   '80 80 231';
PRT.ColorMed_a74 =   '0 0 168';
PRT.ColorMed_resp = '0 0 0'; % color for response

% for ifat's grant graph
PRT.ColorMon_r25 =    '254 0 3';
PRT.ColorMon_r50 =    '0 43 209';
PRT.ColorMon_r75 =   '105 255 42';
% PRT.ColorMon_r25 =    '255 152 0';
% PRT.ColorMon_r50 =    '153 91 0';
% PRT.ColorMon_r75 =   '76 45 0';
PRT.ColorMon_a24 =   '0 229 96';
PRT.ColorMon_a50 =   '0 127 53';
PRT.ColorMon_a74 =   '0 76 32';
PRT.ColorMon_resp = '0 0 0'; % color for response

%% Run for all of the above
% Iterate for each subject, each domain (each _fitpar data file because loss and gains are separate)
for i = 1:length(SubjectNums)
    for j = 1:length(domain)
        for k = 1:length(ParametricModType)
            % for non parametric design matrics((6magnitude)*2domains=12 lottery duration + 1 response)
            if strcmp(ParametricModType{k}, 'none') 
                PRT.NrOfConditions =      '13';
%           else
%             % for parametic design matrics, 8+4 conditions (add a parametric modulator for each empty predictor, the other four will  be added when SDM is created)
%                 PRT.NrOfConditions =      '12'; % not sure if this situation need parametirc modulators
            end
            if (strcmp(ParametricModType{k},'SV')|strcmp(ParametricModType{k},'RewardValue')) & strcmp(domain{j},'med') % for medical blocks, these two parameters are meanlingless
            else
                PTB_MDM_Protocol_Gen_ver5_2(SubjectNums(i), domain{j}, ...
                    tr, trialduration, duration2use, DiscardedAcquisition, ...
                    ParametricModType{k}, ...
                    path_in, path_out, PRT,tpb)
            end
        end
    end
end
