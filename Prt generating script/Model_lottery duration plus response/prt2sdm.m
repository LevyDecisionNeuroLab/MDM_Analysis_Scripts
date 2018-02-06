% PRT2SDM Take all PRT files and re-make them into SDM files
%
% NOTE: If you get an error about `checkstruct` called in `CreateSDM`,
% this is caused by Neuroelf's function being named identically to
% MATLAB's Finance Toolkit's function. Remove the Finance Toolkit from
% path to solve the problem.
clear all
close all

%% Input
prtwave = 'Prt_Model_disp_resp_08222017';
sdmwave = 'SDM_Model_disp_resp_08222017';

%% Load + save settings
root_path = 'D:\Ruonan\Projects in the lab\MDM Project\Medical Decision Making Imaging\MDM_imaging\Behavioral Analysis';
prt_loc = [root_path filesep 'Prt files' filesep prtwave filesep];
sdm_loc = [root_path filesep 'SDM files' filesep sdmwave filesep];

if exist(sdm_loc)==0
    mkdir(fullfile(root_path,'SDM files'),sdmwave)
end

n = neuroelf;
prts = n.findfiles(prt_loc, '*.prt');

for i = 1:length(prts)
    prt_filename = prts{i};
    prt = xff(prt_filename);

    sdm = prt.CreateSDM(struct('nvol', 300, ... 
            'prtr', 1000, ...
            'rcond', []));

    [~, name, ~] = fileparts(prt_filename); % Only filename, without path or extension
    sdm.SaveAs([sdm_loc name '.sdm']);

    % The operation takes time and operates on many files -> 
    fprintf('%.2f%%: Generated %s\n', 100 * i / length(prts), name);
end
