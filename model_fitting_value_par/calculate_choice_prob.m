% calculate model-free choice probability

clearvars
close all

%% Set up loading & subject selection
root = 'D:\Ruonan\Projects in the lab\MDM Project\Medical Decision Making Imaging\MDM_imaging\Behavioral Analysis';
data_path = fullfile(root, 'PTB Behavior Log/'); % root of folders is sufficient
% rating_filename = fullfile(root, 'Behavior Analysis/MDM_Rating.csv');
out_path = fullfile(root, 'Behavior Analysis');

choiceFile_mon = fullfile(out_path, 'choice_data_mon_11122019.csv');
choiceFile_med = fullfile(out_path, 'choice_data_med_11122019.csv');
countFile_mon = fullfile(out_path, 'choice_count_data_mon_11122019.csv');
countFile_med = fullfile(out_path, 'chocie_count_data_med_11122019.csv');
cpbylevelFile = fullfile(out_path, 'nonpar_11122019.csv');

% cp_title = ['id','is_med', 'r25', 'r50','r75','rAll','a24', 'a50','a74','aAll','a24_r50', 'a50_r50','a74_r50','a_r50 All','All','error'];
% dlmwrite(cpbylevelFile, cp_title, 'coffset', 1, '-append', 'delimiter', ',');    

addpath(genpath(data_path)); % generate path for all the subject data folder

% get all subjects number in the folder
subjects = getSubjectsInDir(data_path, 'subj');
exclude = [2581]; % TEMPORARY: subjects incomplete data (that the script is not ready for)
subjects = subjects(~ismember(subjects, exclude));
% subjects = [2585];

% all values
vals_mon = [5,8,12,25];
vals_med = [1,2,3,4];

for subj_idx = 1:length(subjects)
  domains = {'MON', 'MED'};

  for domain_idx = 1:length(domains)
    
    subjectNum = subjects(subj_idx);
    domain = domains{domain_idx};
    
    % directly use load/save violate transparency, use a function instead 
    Data = load_mat(subjectNum, domain);

    %% Refine variables
    
    % exlucde missing response
    include_indices = Data.choice ~= 0;

    choice = Data.choice(include_indices);
    values = Data.vals(include_indices);
    ambigs = Data.ambigs(include_indices);
    probs  = Data.probs(include_indices);
    
    % Side with lottery is counterbalanced across subjects 
    % code 0 as reference choice, 1 as lottery choice
    % if sum(choice == 2) > 0 % Only if choice has not been recoded yet. RJ-Not necessary
    % RJ-If subject do not press 2 at all, the above if condition is problematic
      if Data.refSide == 2
          choice(choice == 2) = 0;
          choice(choice == 1) = 1;
      elseif Data.refSide == 1 % Careful: rerunning this part will make all choices 0
          choice(choice == 1) = 0;
          choice(choice == 2) = 1;
      end
    
    % choice data for $5 only, risky and ambiguous trials, for stochastic
    % dominance error
    idx_only5 = and(Data.choice ~= 0, Data.vals' == 5);
    choice5 = Data.choice(idx_only5);
    values5 = Data.vals(idx_only5);
    ambigs5 = Data.ambigs(idx_only5);
    probs5  = Data.probs(idx_only5);
    
    if Data.refSide == 2
        choice5(choice5 == 2) = 0;
        choice5(choice5 == 1) = 1;
    elseif Data.refSide == 1 % Careful: rerunning this part will make all choices 0
        choice5(choice5 == 1) = 0;
        choice5(choice5 == 2) = 1;
    end
    
    choice_prob_5= sum(choice5)/length(choice5);
    
    %% Create choice matrices
    % One matrix per condition. Matrix values are binary (0 for sure
    % choice, 1 for lottery). Matrix dimensions are prob/ambig-level
    % x payoff values. Used for graphing and some Excel exports.
 
    choiceMatrix = create_choice_matrix(values,ambigs,probs,choice);
    
    
    %% print choice prob and count by lottery type
    
    % choice probability
    riskyChoices = choiceMatrix.riskProb;
    ambigChoices = choiceMatrix.ambigProb;
    % count of trials (with responses)
    riskyChoicesC = choiceMatrix.riskCount;
    ambigChoicesC = choiceMatrix.ambigCount;  
    
    
    choices_all = [riskyChoices; ambigChoices];
    counts_all = [riskyChoicesC; ambigChoicesC];
    
    if strcmp(domain, 'MON')
        all_data_subject = [vals_mon; choices_all];
        dlmwrite(choiceFile_mon, subjectNum , '-append', 'roffset', 1, 'delimiter', ',');  
        dlmwrite(choiceFile_mon, all_data_subject, 'coffset', 1, '-append', 'delimiter', ',');
        all_data_subjectC = [vals_mon; counts_all];
        dlmwrite(countFile_mon, subjectNum , '-append', 'roffset', 1, 'delimiter', ',');  
        dlmwrite(countFile_mon, all_data_subjectC, 'coffset', 1, '-append', 'delimiter', ',');

    elseif strcmp(domain, 'MED')
        all_data_subject = [vals_med; choices_all];
        dlmwrite(choiceFile_med, subjectNum , '-append', 'roffset', 1, 'delimiter', ',');  
        dlmwrite(choiceFile_med, all_data_subject, 'coffset', 1, '-append', 'delimiter', ',');
        all_data_subjectC = [vals_med; counts_all];
        dlmwrite(countFile_med, subjectNum , '-append', 'roffset', 1, 'delimiter', ',');  
        dlmwrite(countFile_med, all_data_subjectC, 'coffset', 1, '-append', 'delimiter', ',');
        
    end    
     
        
    %% print choice prob by uncertainty level
    
    %exclude choices with $5 or slight improvement
    %P is monetary, N is medical
    riskyChoices_no5 = riskyChoices(:,2:size(riskyChoices,2));
    riskyChoicesC_no5 = riskyChoicesC(:,2:size(riskyChoicesC,2));
    ambigChoices_no5 = ambigChoices(:,2:size(ambigChoices,2));
    ambigChoicesC_no5 = ambigChoicesC(:,2:size(ambigChoicesC,2));

    % monetary
    riskyChoicesT = riskyChoices_no5 .* riskyChoicesC_no5; % choice total counts = choice prob * trial counts
    cpByRisk = zeros(size(riskyChoices_no5,1),1); % choice prob by risk level
    for i = 1:size(cpByRisk,1)
      cpByRisk(i) = sum(riskyChoicesT(i,:))/sum(riskyChoicesC_no5(i,:));
    end
    cpRiskAll = sum(riskyChoicesT(:))/sum(riskyChoicesC_no5(:));
    
    ambigChoicesT = ambigChoices_no5 .* ambigChoicesC_no5; % choice total counts = choice prob * tial counts
    cpByAmbig = zeros(size(ambigChoices_no5,1),1); % choice prob by ambig level
    for i = 1:size(cpByAmbig,1)
      cpByAmbig(i) = sum(ambigChoicesT(i,:))/sum(ambigChoicesC_no5(i,:));
    end
    cpAmbigAll = sum(ambigChoicesT(:))/sum(ambigChoicesC_no5(:));
    
    ambigAtt = ambigChoices_no5 - riskyChoices_no5(2,:); 
    ambigAttByAmbig = nanmean(ambigAtt.'); % model free ambig attitude by ambig level
    ambigAttAll = nanmean(ambigAtt(:));
    
    AllT = sum(riskyChoicesT(:)) + sum(ambigChoicesT(:)); % choice total counts = choice prob * tial counts
    AllC = sum(riskyChoicesC_no5(:)) + sum(ambigChoicesC_no5(:));
    cpAll = sum(AllT(:))/AllC;
    

    if strcmp(domain, 'MON')
        is_med = 0;
    elseif strcmp(domain, 'MED')
        is_med = 1;
    end
    
    cp_title = ['id','is_med', 'r25', 'r50','r75','rAll','a24', 'a50','a74','aAll','a24_r50', 'a50_r50','a74_r50','a_r50 All','All','error']; % do not print. does not work for dlmwrite
    cp_data_subject = [subjectNum,is_med,cpByRisk.',cpRiskAll,cpByAmbig.',cpAmbigAll,ambigAttByAmbig,ambigAttAll,cpAll,choice_prob_5];
                          
    dlmwrite(cpbylevelFile, cp_data_subject, 'coffset', 1, '-append', 'delimiter', ',');    
    
  end
  
end





