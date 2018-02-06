%% This script is meant to take parameters files and print some data to Excel
clearvars
close all

fitpar_wave = 'financial fitpar 1204022017';
fitbyrating = false;

root = 'D:\Ruonan\Projects in the lab\MDM Project\Medical Decision Making Imaging\MDM_imaging\Behavioral Analysis\';
data_path = fullfile(root, 'PTB Behavior Log/'); % Original log from PTB
subjects = getSubjectsInDir(data_path, 'subj'); %function
exclude = [2581]; % TEMPORARY: subjects incomplete data (that the script is not ready for)
subjects = subjects(~ismember(subjects, exclude));

path = fullfile(root, 'Behavior data finacial fitpar', fitpar_wave,filesep);
cd(path);

% defining monetary values
valueP = [5 8 12 25];

output_file1 = ['param_' fitpar_wave '.txt'];
% might not need
output_file2 = 'choice_data.txt'

% results file
fid1 = fopen([output_file1],'w')


if ~fitbyrating
    fprintf(fid1,'subject\tmonetary\n')
    fprintf(fid1,'\talpha\talphase\tbeta\tbetase\tgamma\tgammase\tLL\tr2\tAIC\tBIC\tmodel\texitFlag\toptimizer\n')
else
    fprintf(fid1,'subject\tmonetary\t\t\t\t\t\t\t\t\t\t\t\t\tmedical\n')
    fprintf(fid1,'\talpha\talphase\tbeta\tbetase\tgamma\tgammase\tLL\tr2\tAIC\tBIC\tmodel\texitFlag\toptimizer\talpha\talphase\tbeta\tbetase\tgamma\tgammase\tLL\tr2\tAIC\tBIC\tmodel\texitFlag\toptimizer\n')
end

% Fill in subject numbers separated by commas
% subjects = {'87','88'};
for s = 1:length(subjects)
    
    subject = subjects(s); 
    % load mon file for subject and extract params & choice data
    load(['MDM_MON_' num2str(subject) '_financ_fitpar.mat']);
    alphaP = Datamon.alpha;
    alphaseP = Datamon.MLE.se(3);
    betaP = Datamon.beta;
    betaseP = Datamon.MLE.se(2);
    gammaP = Datamon.gamma;
    gammaseP = Datamon.MLE.se(1);
    LLP = Datamon.MLE.LL;
    r2P = Datamon.MLE.r2;
    AICP = Datamon.MLE.AIC;
    BICP = Datamon.MLE.BIC;
    modelP = Datamon.MLE.model;
    exitFlagP = Datamon.MLE.exitflag;
    optimizerP = Datamon.MLE.optimizer;
    
%     riskyChoicesP = Datamon.riskyChoices;
%     ambigChoicesP = Datamon.ambigChoices;
    
    if fitbyrating
        % load med file for subject and extract params & choice data
        load(['MDM_MED_' num2str(subject) '_financ_fitpar.mat']);
        alphaN = Datamed.alpha;
        alphaseN = Datamed.MLE.se(3);
        betaN = Datamed.beta;
        betaseN = Datamed.MLE.se(2);
        gammaN = Datamed.gamma;
        gammaseN = Datamed.MLE.se(1);
        LLN = Datamed.MLE.LL;
        r2N = Datamed.MLE.r2;
        AICN = Datamed.MLE.AIC;
        BICN = Datamed.MLE.BIC;
        modelN = Datamed.MLE.model;
        exitFlagN = Datamed.MLE.exitflag;
        optimizerN = Datamed.MLE.optimizer;

        %write into param text file
        fprintf(fid1,'%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%f\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%f\t%s\n',...
                num2str(subject),alphaP,alphaseP,betaP,betaseP,gammaP,gammaseP,LLP,r2P,AICP,BICP,modelP,exitFlagP,optimizerP,...
                                 alphaN,alphaseN,betaN,betaseN,gammaN,gammaseN,LLN,r2N,AICN,BICN,modelN,exitFlagN,optimizerN)
    else
        %write into param text file
        fprintf(fid1,'%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%f\t%s\n',...
                num2str(subject),alphaP,alphaseP,betaP,betaseP,gammaP,gammaseP,LLP,r2P,AICP,BICP,modelP,exitFlagP,optimizerP)
    end
    
%     % for Excel file - choice data
%     
%     % Firt, combine choice data with and without $4
%     choices_allP = [riskyChoicesP; ambigChoicesP];
% %     choices_allP = [choices4P,choices_rnaP];
%     choices_allN = [riskyChoicesN; ambigChoicesN];
% %     choices_allN = [choices4N,choices_rnaN];
%     
%     all_data_subject = [valueP; choices_allP ;valueP; choices_allN];
%     
%     xlFile = [path 'choice_data.xls'];
%     dlmwrite(xlFile, subject , '-append', 'roffset', 1, 'delimiter', ' ');  
%     dlmwrite(xlFile, all_data_subject, 'coffset', 1, '-append', 'delimiter', '\t');
%     
end

fclose(fid1);
