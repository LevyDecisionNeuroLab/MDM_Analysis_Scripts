% This script is meant to take parameters files and print some data to Excel
clearvars
close all

%% Define conditions
fitparwave = 'Behavior data fitpar_09300219';
includeAmbig = true;
model = 'ambigSVPar';

%% Setup
root = 'D:\Ruonan\Projects in the lab\MDM Project\Medical Decision Making Imaging\MDM_imaging\Behavioral Analysis';
function_path = fullfile(root,'MDM_Analysis_Scripts','Model fitting script');
addpath(function_path)
data_path = fullfile(root, 'PTB Behavior Log/'); % Original log from PTB
subjects = getSubjectsInDir(data_path, 'subj'); %function
exclude = [2581]; % TEMPORARY: subjects incomplete data (that the script is not ready for)
subjects = subjects(~ismember(subjects, exclude));
% subjects = [2073 2582 2587 2597 2651 2663 2665 2666 2550 2585 2596 2600 2655 2659 2660 2664];

path = fullfile(root, 'Behavior fitpar files', fitparwave,filesep);
cd(path);

% defining monetary values
valueP = [5 8 12 25];

output_file1 = ['param_' fitparwave '.txt'];
% might not need
% output_file2 = 'choice_data.txt';
% output_file3 = 'choice_prob.txt';

% results file
fid1 = fopen([output_file1],'w')

% fprintf(fid1,'subject\tmonetary\t\t\t\t\t\t\t\t\t\t\t\t\t\tmedical\n')
fprintf(fid1,'id\talpha_mon\tbeta_mon\tgamma_mon\tval1_mon\tval2_mon\tval3_mon\tval4_mon\tLL_mon\tr2_adj_mon\tAIC_mon\tBIC_mon\tmodel_mon\texitFlag_mon\toptimizer_mon\talpha_med\tbeta_med\tgamma_med\tval1_med\tval2_med\tval3_med\tval4_med\tLL_med\tr2_adj_med\tAIC_med\tBIC_med\tmodel_med\texitFlag_med\toptimizer_med\n')


% Fill in subject numbers separated by commas
% subjects = {'87','88'};
for s = 1:length(subjects)
%% Print parameters
    subject = subjects(s); 
    % load mon file for subject and extract params & choice data
    load(['MDM_MON_' num2str(subject) '_fitpar.mat']);
    choiceMatrixP = Datamon.choiceMatrix;

    riskyChoicesP = choiceMatrixP.riskProb;
    ambigChoicesP = choiceMatrixP.ambigProb;
    riskyChoicesPC = choiceMatrixP.riskCount;
    ambigChoicesPC = choiceMatrixP.ambigCount;
    
%     svByLottP = Datamon.svByLott;
%     svRefP = Datamon.svRef;

    betaP = Datamon.beta;
%     betaseP = Datamon.MLE.se(2);
    gammaP = Datamon.gamma;
%     gammaseP = Datamon.MLE.se(1);

    if strcmp(model, 'ambigNriskValPar')
        alphaP = Datamon.alpha;
        %     alphaseP = Datamon.MLE.se(3);
    elseif strcmp(model, 'ambigSVPar')
         alphaP = NaN;
    end
    
    val1P = Datamon.val_par(1);
    %     val1seP = Datamon.MLE.se(4);
    val2P = Datamon.val_par(2);
    %     val2seP = Datamon.MLE.se(5);
    val3P = Datamon.val_par(3);
    %     val3seP = Datamon.MLE.se(6);
    val4P = Datamon.val_par(4);
    %     val4seP = Datamon.MLE.se(7);      
    LLP = Datamon.MLE.LL;
    r2_adjP = Datamon.MLE.r2_adj;
    AICP = Datamon.MLE.AIC;
    BICP = Datamon.MLE.BIC;
    modelP = Datamon.MLE.model;
    exitFlagP = Datamon.MLE.exitflag;
%     optimizerP = Datamon.MLE.optimizer;
    optimizerP = 'BADS';
    

    % load med file for subject and extract params & choice data
    load(['MDM_MED_' num2str(subject) '_fitpar.mat']);
    choiceMatrixN = Datamed.choiceMatrix;

    riskyChoicesN = choiceMatrixN.riskProb;
    ambigChoicesN = choiceMatrixN.ambigProb;
    riskyChoicesNC = choiceMatrixN.riskCount;
    ambigChoicesNC = choiceMatrixN.ambigCount;

    if strcmp(model, 'ambigNriskValPar')
        alphaN = Datamed.alpha;
        %     alphaseP = Datamon.MLE.se(3);
    elseif strcmp(model, 'ambigSVPar')
         alphaN = NaN;
    end
    
    betaN = Datamed.beta;
%     betaseN = Datamed.MLE.se(2);
    gammaN = Datamed.gamma;
%     gammaseN = Datamed.MLE.se(1);
    val1N = Datamed.val_par(1);
%     val1seN = Datamed.MLE.se(4);
    val2N = Datamed.val_par(2);
%     val2seN = Datamed.MLE.se(5);
    val3N = Datamed.val_par(3);
%     val3seN = Datamed.MLE.se(6);
    val4N = Datamed.val_par(4);
%     val4seN = Datamed.MLE.se(7);
    LLN = Datamed.MLE.LL;
    r2_adjN = Datamed.MLE.r2_adj;
    AICN = Datamed.MLE.AIC;
    BICN = Datamed.MLE.BIC;
    modelN = Datamed.MLE.model;
    exitFlagN = Datamed.MLE.exitflag;
%     optimizerN = Datamed.MLE.optimizer; 
    optimizerN = 'BADS';

    %write into param text file
    fprintf(fid1,'%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%d\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%d\t%s\n',...
        num2str(subject),alphaP,betaP,gammaP,val1P,val2P,val3P,val4P,LLP,r2_adjP,AICP,BICP,modelP,exitFlagP,optimizerP,...
                         alphaN,betaN,gammaN,val1N,val2N,val3N,val4N,LLN,r2_adjN,AICN,BICN,modelN,exitFlagN,optimizerN);                     
  
    %% print true and model fitting choice prob by trial types
%     choiceMatrixModelP = create_choice_matrix(Datamon.vals,Datamon.ambigs,Datamon.probs,Datamon.choiceModeled);
%     
%     % Plot risky choice prob
%     screensize = get(groot, 'Screensize');
%     fig = figure('Position', [screensize(3)/4 screensize(4)/22 screensize(3)/2 screensize(4)*10/12]);
%     ax1 = subplot(2,1,1);
%     % actual choice
%     barplot = bar(ax1,[choiceMatrixP.riskProb(1,:),choiceMatrixP.riskProb(2,:),choiceMatrixP.riskProb(3,:)],'FaceColor','y');
%     hold on
%     xticklabels({'r25-$5','r25-$8','r25-$12','r25-$25','r50-$5','r50-$8','r50-$12','r50-$25','r75-$5','r75-$8','r75-$12','r75-$25'})
%     ylim([0 1.1])
%     yticks([0:0.1:1.1])
%     % predicted by model
%     plot(ax1,[choiceMatrixModelP.riskProb(1,:),choiceMatrixModelP.riskProb(2,:),choiceMatrixModelP.riskProb(3,:)],'LineStyle','none','Marker','o' )
%     title(ax1,['Subject ' num2str(subject) ' monetary risky choice probability, model:' Datamon.MLE.model])
%     
%     if includeAmbig
%         % Plot ambiguous choice prob
%         ax2 = subplot(2,1,2);
%         % actual choice
%         barplot = bar(ax2,[choiceMatrixP.ambigProb(1,:),choiceMatrixP.ambigProb(2,:),choiceMatrixP.ambigProb(3,:)],'FaceColor','y');
%         hold on
%         xticklabels({'a24-$5','a24-$8','a24-$12','a24-$25','a50-$5','a50-$8','a50-$12','a50-$25','a74-$5','a74-$8','a74-$12','a74-$25'})
%         ylim([0 1.1])
%         yticks([0:0.1:1.1])
%         % predicted by model
%         plot(ax2,[choiceMatrixModelP.ambigProb(1,:),choiceMatrixModelP.ambigProb(2,:),choiceMatrixModelP.ambigProb(3,:)],'LineStyle','none','Marker','o' )
%         title(ax2,['Subject ' num2str(subject) ' monetary ambiguous choice probability, model:' Datamon.MLE.model])
%     
%         % save figure
%         saveas(fig,['Subject ' num2str(subject), ' mon choice prob-' Datamon.MLE.model])
%     end
%     
%     choiceMatrixModelN = create_choice_matrix(Datamed.vals,Datamed.ambigs,Datamed.probs,Datamed.choiceModeled);
% 
%     % Plot risky choice prob
%     screensize = get(groot, 'Screensize');
%     fig = figure('Position', [screensize(3)/4 screensize(4)/22 screensize(3)/2 screensize(4)*10/12]);
%     ax1 = subplot(2,1,1);
%     % actual choice
%     barplot = bar(ax1,[choiceMatrixN.riskProb(1,:),choiceMatrixN.riskProb(2,:),choiceMatrixN.riskProb(3,:)],'FaceColor','y');
%     hold on
%     xticklabels({'r25-sl','r25-mod','r25-maj','r25-rec','r50-sl','r50-mod','r50-maj','r50-rec','r75-sl','r75-mod','r75-maj','r75-rec'})
%     ylim([0 1.1])
%     yticks([0:0.1:1.1])
%     % predicted by model
%     plot(ax1,[choiceMatrixModelN.riskProb(1,:),choiceMatrixModelN.riskProb(2,:),choiceMatrixModelN.riskProb(3,:)],'LineStyle','none','Marker','o' )
%     title(ax1,['Subject ' num2str(subject) ' medical risky choice probability, model:' Datamed.MLE.model])
% 
%     if includeAmbig
%         % Plot ambiguous choice prob
%         ax2 = subplot(2,1,2);
%         % actual choice
%         barplot = bar(ax2,[choiceMatrixN.ambigProb(1,:),choiceMatrixN.ambigProb(2,:),choiceMatrixN.ambigProb(3,:)],'FaceColor','y');
%         hold on
%         xticklabels({'a24-sl','a24-mod','a24-maj','a24-rec','a50-sl','a50-mod','a50-maj','a50-rec','a74-sl','a74-mod','a74-maj','a74-rec'})
%         ylim([0 1.1])
%         yticks([0:0.1:1.1])
%         % predicted by model
%         plot(ax2,[choiceMatrixModelN.ambigProb(1,:),choiceMatrixModelN.ambigProb(2,:),choiceMatrixModelN.ambigProb(3,:)],'LineStyle','none','Marker','o' )
%         title(ax2,['Subject ' num2str(subject) ' medical ambiguous choice probability, model:' Datamed.MLE.model])
%     end
% 
%     % save figure
%     saveas(fig,['Subject ' num2str(subject), ' med choice prob-' Datamed.MLE.model])

    
    %% for Excel file - subjective values
%     xlFile = ['SV_unconstrained_by_lottery.xls'];
%     dlmwrite(xlFile, subject, '-append', 'roffset', 1, 'delimiter', ' '); 
%     dlmwrite(xlFile, svRefUncstr, '-append', 'coffset', 1, 'delimiter', '\t');
%     dlmwrite(xlFile, svUncstrByLott, 'coffset', 1, '-append', 'delimiter', '\t');
   

end

fclose(fid1);
