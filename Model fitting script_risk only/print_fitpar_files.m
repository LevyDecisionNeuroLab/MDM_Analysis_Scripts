% This script is meant to take parameters files and print some data to Excel
clearvars
close all

%% Define conditions
fitparwave = 'Behavior data fitpar_1220042017';
fitbyrating = false;

%% Setup
root = 'D:\Ruonan\Projects in the lab\MDM Project\Medical Decision Making Imaging\MDM_imaging\Behavioral Analysis';
function_path = fullfile(root,'MDM_Analysis_Scripts','Model fitting script_risk only');
addpath(function_path)
data_path = fullfile(root, 'PTB Behavior Log/'); % Original log from PTB
subjects = getSubjectsInDir(data_path, 'subj'); %function
exclude = [2581]; % TEMPORARY: subjects incomplete data (that the script is not ready for)
subjects = subjects(~ismember(subjects, exclude));
% subjects = [2585];

path = fullfile(root, 'Behavior fitpar files', fitparwave,filesep);
cd(path);

% defining monetary values
valueP = [5 8 12 25];

output_file1 = ['param_unconstrained_' fitparwave '.txt'];
% might not need
% output_file2 = 'choice_data.txt';
% output_file3 = 'choice_prob.txt';

% results file
fid1 = fopen([output_file1],'w')

if ~fitbyrating
    fprintf(fid1,'subject\tmonetary\n')
    fprintf(fid1,'\talpha\talphase\tgamma\tgammase\tLL\tr2\tAIC\tBIC\tmodel\texitFlag\toptimizer\n')
else
    fprintf(fid1,'subject\tmonetary\t\t\t\t\t\t\t\t\tmedical\n')
    fprintf(fid1,'\talpha\talphase\tgamma\tgammase\tLL\tr2\tAIC\tBIC\tmodel\texitFlag\toptimizer\talpha\talphase\tgamma\tgammase\tLL\tr2\tAIC\tBIC\tmodel\texitFlag\toptimizer\n')
end

% Fill in subject numbers separated by commas
% subjects = {'87','88'};
for s = 1:length(subjects)
%% Print parameters
    subject = subjects(s); 
    % load mon file for subject and extract params & choice data
    load(['MDM_MON_' num2str(subject) '_fitpar.mat']);
    choiceMatrixP = Datamon.choiceMatrix;

    riskyChoicesP = choiceMatrixP.riskProb;
    riskyChoicesPC = choiceMatrixP.riskCount;

    alphaP = Datamon.alpha;
    alphaseP = Datamon.MLE.se(2);
    gammaP = Datamon.gamma;
    gammaseP = Datamon.MLE.se(1);
    LLP = Datamon.MLE.LL;
    r2P = Datamon.MLE.r2;
    AICP = Datamon.MLE.AIC;
    BICP = Datamon.MLE.BIC;
    modelP = Datamon.MLE.model;
    exitFlagP = Datamon.MLE.exitflag;
    optimizerP = Datamon.MLE.optimizer;
    

    % load med file for subject and extract params & choice data
    load(['MDM_MED_' num2str(subject) '_fitpar.mat']);
    choiceMatrixN = Datamed.choiceMatrix;

    riskyChoicesN = choiceMatrixN.riskProb;
    riskyChoicesNC = choiceMatrixN.riskCount;

    if fitbyrating
        svByLottN = Datamed.svByLott;
        svRefN = Datamed.svRef;
        alphaN = Datamed.alpha;
        alphaseN = Datamed.MLE.se(2);
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
        fprintf(fid1,'%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%f\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%f\t%s\n',...
                num2str(subject),alphaP,alphaseP,gammaP,gammaseP,LLP,r2P,AICP,BICP,modelP,exitFlagP,optimizerP,...
                                 alphaN,alphaseN,gammaN,gammaseN,LLN,r2N,AICN,BICN,modelN,exitFlagN,optimizerN)
    else
        %write into param text file
        fprintf(fid1,'%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%f\t%s\n',...
                num2str(subject),alphaP,alphaseP,gammaP,gammaseP,LLP,r2P,AICP,BICP,modelP,exitFlagP,optimizerP)
    end
  
    %% print true and model fitting choice prob by trial types
    riskIndicesP = Datamon.ambigs == 0;
    choiceMatrixModelP = create_choice_matrix(Datamon.vals(riskIndicesP),Datamon.probs(riskIndicesP),Datamon.choiceModeled);
    
    % Plot risky choice prob
    screensize = get(groot, 'Screensize');
    fig = figure('Position', [screensize(3)/4 screensize(4)/22 screensize(3)/2 screensize(4)*10/12]);
    % actual choice
    barplot = bar([choiceMatrixP.riskProb(1,:),choiceMatrixP.riskProb(2,:),choiceMatrixP.riskProb(3,:)],'FaceColor','y');
    hold on
    xticklabels({'r25-$5','r25-$8','r25-$12','r25-$25','r50-$5','r50-$8','r50-$12','r50-$25','r75-$5','r75-$8','r75-$12','r75-$25'})
    ylim([0 1.1])
    yticks([0:0.1:1.1])
    % predicted by model
    plot([choiceMatrixModelP.riskProb(1,:),choiceMatrixModelP.riskProb(2,:),choiceMatrixModelP.riskProb(3,:)],'LineStyle','none','Marker','o' )
    title(['Subject ' num2str(subject) ' monetary risky choice probability, model:' Datamon.MLE.model])
    
    % save figure
    saveas(fig,['Subject ' num2str(subject), ' mon choice prob-' Datamon.MLE.model])
    
    if fitbyrating
        riskIndicesN = Datamed.ambigs == 0;
        choiceMatrixModelN = create_choice_matrix(Datamed.vals(riskIndicesN),Datamed.probs(riskIndicesN),Datamed.choiceModeled);

        % Plot risky choice prob
        screensize = get(groot, 'Screensize');
        fig = figure('Position', [screensize(3)/4 screensize(4)/22 screensize(3)/2 screensize(4)*10/12]);
        % actual choice
        barplot = bar([choiceMatrixN.riskProb(1,:),choiceMatrixN.riskProb(2,:),choiceMatrixN.riskProb(3,:)],'FaceColor','y');
        hold on
        xticklabels({'r25-sl','r25-mod','r25-maj','r25-rec','r50-sl','r50-mod','r50-maj','r50-rec','r75-sl','r75-mod','r75-maj','r75-rec'})
        ylim([0 1.1])
        yticks([0:0.1:1.1])
        % predicted by model
        plot([choiceMatrixModelN.riskProb(1,:),choiceMatrixModelN.riskProb(2,:),choiceMatrixModelN.riskProb(3,:)],'LineStyle','none','Marker','o' )
        title(['Subject ' num2str(subject) ' medical risky choice probability, model:' Datamed.MLE.model])

        % save figure
        saveas(fig,['Subject ' num2str(subject), ' med choice prob-' Datamed.MLE.model])

    end
    
    
    %% for Excel file - choice prob by lottery type
    
%     % Firt, combine choice data with and without $4
%     choices_allP = [riskyChoicesP; ambigChoicesP];
%     choices_allN = [riskyChoicesN; ambigChoicesN];
%     
%     all_data_subject = [valueP; choices_allP ;valueP; choices_allN];
%     
%     xlFile = ['choice_data.xls'];
%     dlmwrite(xlFile, subject , '-append', 'roffset', 1, 'delimiter', ' ');  
%     dlmwrite(xlFile, all_data_subject, 'coffset', 1, '-append', 'delimiter', '\t');
    
    %% for Excel file - subjective values
%     xlFile = ['SV_unconstrained_by_lottery.xls'];
%     dlmwrite(xlFile, subject, '-append', 'roffset', 1, 'delimiter', ' '); 
%     dlmwrite(xlFile, svRefUncstr, '-append', 'coffset', 1, 'delimiter', '\t');
%     dlmwrite(xlFile, svUncstrByLott, 'coffset', 1, '-append', 'delimiter', '\t');
    
    %% for Excel file - choice prob by uncertainty level
    
%     %exclude choices with $5 or slight improvement
%     %P is monetary, N is medical
%     riskyChoicesP = riskyChoicesP(:,2:size(riskyChoicesP,2));
%     riskyChoicesPC = riskyChoicesPC(:,2:size(riskyChoicesPC,2));
%     ambigChoicesP = ambigChoicesP(:,2:size(ambigChoicesP,2));
%     ambigChoicesPC = ambigChoicesPC(:,2:size(ambigChoicesPC,2));
%     riskyChoicesN = riskyChoicesN(:,2:size(riskyChoicesN,2));
%     riskyChoicesNC = riskyChoicesNC(:,2:size(riskyChoicesNC,2));
%     ambigChoicesN = ambigChoicesN(:,2:size(ambigChoicesN,2));
%     ambigChoicesNC = ambigChoicesNC(:,2:size(ambigChoicesNC,2));
% 
%     % monetary
%     riskyChoicesPT = riskyChoicesP .* riskyChoicesPC; % choice total counts = choice prob * trial counts
%     cpByRiskP = zeros(size(riskyChoicesP,1),1); % choice prob by risk level
%     for i = 1:size(cpByRiskP,1);
%       cpByRiskP(i) = sum(riskyChoicesPT(i,:))/sum(riskyChoicesPC(i,:));
%     end
%     cpRiskAllP = sum(riskyChoicesPT(:))/sum(riskyChoicesPC(:));
%     
%     ambigChoicesPT = ambigChoicesP .* ambigChoicesPC; % choice total counts = choice prob * tial counts
%     cpByAmbigP = zeros(size(ambigChoicesP,1),1); % choice prob by ambig level
%     for i = 1:size(cpByAmbigP,1);
%       cpByAmbigP(i) = sum(ambigChoicesPT(i,:))/sum(ambigChoicesPC(i,:));
%     end
%     cpAmbigAllP = sum(ambigChoicesPT(:))/sum(ambigChoicesPC(:));
%     
%     ambigAttP = ambigChoicesP - riskyChoicesP(2,:); 
%     ambigAttByAmbigP = nanmean(ambigAttP.'); % model free ambig attitude by ambig level
%     ambigAttAllP = nanmean(ambigAttP(:));
%     
%     AllPT = sum(riskyChoicesPT(:)) + sum(ambigChoicesPT(:)); % choice total counts = choice prob * tial counts
%     AllPC = sum(riskyChoicesPC(:)) + sum(ambigChoicesPC(:));
%     cpAllP = sum(AllPT(:))/AllPC;
% 
%    
%     
%     %Medical
%     riskyChoicesNT = riskyChoicesN .* riskyChoicesNC; % choice total counts = choice prob * tial counts
%     cpByRiskN = zeros(size(riskyChoicesN,1),1); % choice prob by risk level
%     for i = 1:size(cpByRiskN,1);
%       cpByRiskN(i) = sum(riskyChoicesNT(i,:))/sum(riskyChoicesNC(i,:));
%     end
%     cpRiskAllN = sum(riskyChoicesNT(:))/sum(riskyChoicesNC(:));
%     
%     ambigChoicesNT = ambigChoicesN .* ambigChoicesNC; % choice total counts = choice prob * tial counts
%     cpByAmbigN = zeros(size(ambigChoicesN,1),1); % choice prob by ambig level
%     for i = 1:size(cpByAmbigN,1);
%       cpByAmbigN(i) = sum(ambigChoicesNT(i,:))/sum(ambigChoicesNC(i,:));
%     end
%     cpAmbigAllN = sum(ambigChoicesNT(:))/sum(ambigChoicesNC(:));
%     
%     ambigAttN = ambigChoicesN - riskyChoicesN(2,:); 
%     ambigAttByAmbigN = nanmean(ambigAttN.'); % model free ambig attitude by ambig lavel
%     ambigAttAllN = nanmean(ambigAttN(:));
%     
%     AllNT = sum(riskyChoicesNT(:)) + sum(ambigChoicesNT(:)); % choice total counts = choice prob * tial counts
%     AllNC = sum(riskyChoicesNC(:)) + sum(ambigChoicesNC(:));
%     cpAllN = sum(AllNT(:))/AllNC;
% 
%     
% 
%    cp_title = {'subject ID', 'r25', 'r50','r75','rAll','a24', 'a50','a74','aAll','a24-r50', 'a50-r50','a74-r50','a-r50 All','All',...
%                              'r25', 'r50','r75','rAll','a24', 'a50','a74','aAll','a24-r50', 'a50-r50','a74-r50','a-r50 All','All'};
%    cp_data_subject = [subject,cpByRiskP.',cpRiskAllP,cpByAmbigP.',cpAmbigAllP,ambigAttByAmbigP,ambigAttAllP,cpAllP...
%                               cpByRiskN.',cpRiskAllN,cpByAmbigN.',cpAmbigAllN,ambigAttByAmbigN,ambigAttAllN,cpAllN];
%                           
%     xlFile = ['choice_prob_without5.xls'];
%     dlmwrite(xlFile, cp_data_subject, 'coffset', 1, '-append', 'delimiter', '\t');

end

fclose(fid1);
