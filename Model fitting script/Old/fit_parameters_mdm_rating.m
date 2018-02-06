% NOTE: Requires MATLAB optim library
%% Changes need to make for MDM analysis
%   -Separate monetary and medical
%   -Monetary should be the same, 
%% Set up loading + subject selection
% TODO: Maybe grab & save condition somewhere?
clear all
close all

root = 'D:\Ruonan\Projects in the lab\MDM Project\Medical Decision Making Imaging\MDM_imaging\Analysis';
data_path = fullfile(root, 'PTB Behavior Log/'); % root of folders is sufficient
rating_filename = fullfile(root, 'Behavior Analysis/MDM_Rating.csv');
fitpar_out_path = fullfile(root, 'Behavior data fitpar_01292017/');

addpath(genpath(data_path)); % generate path for all the subject data folder

% get all subjects number in the folder
subjects = getSubjectsInDir(data_path, 'subj');
exclude = []; % TEMPORARY: subjects incomplete data (that the script is not ready for)
subjects = subjects(~ismember(subjects, exclude));
%  subjects = [9999];

%% load the subjective rating
% column1-subj ID, c2-$0, c3-$5,c4-$8,c5-$12,c6-$25,c7-no effect, c8-slight,c9-moderate,c10-major,c11-recovery.
rating = csvread(rating_filename,1,0); %reads data from the file starting at row offset R1 and column offset C1. For example, the offsets R1=0, C1=0 specify the first value in the file.

for subj_idx = 1:length(subjects)
  domains = {'MON', 'MED'};
  
  for domain_idx = 1:length(domains)
    subjectNum = subjects(subj_idx);
    domain = domains{domain_idx};
    
    fname = sprintf('MDM_%s_%d.mat', domain, subjectNum);
    load(fname) % produces variable `Datamon` or 'Datamed' for convenience, change its name into 'Data'
    sprintf('MDM_%s_%d.mat', domain, subjectNum);
    
    %% Refine variables
    if strcmp(domain, 'MON') ==1;
        Data = Datamon;
        clear Datamon
    else
        Data = Datamed;
        clear Datamed
    end
    
    % Exclude non-responses
    include_indices = Data.choice ~= 0;

    choice = Data.choice(include_indices);
    values = Data.vals(include_indices);
    ambigs = Data.ambigs(include_indices);
    probs  = Data.probs(include_indices);
    
    % Side with lottery is counterbalanced across subjects 
    % -> code 0 as reference choice, 1 as lottery choice
    % TODO: Double-check this is so? - This is true(RJ)
    % TODO: Save in a different variable?
    % if sum(choice == 2) > 0 % Only if choice has not been recoded yet. RJ-Not necessary
    % RJ-If subject do not press 2 at all, the above if condition is problematic
      if Data.refSide == 2
          choice(choice == 2) = 0;
          choice(choice == 1) = 1;
      elseif Data.refSide == 1 % Careful: rerunning this part will make all choices 0
          choice(choice == 1) = 0;
          choice(choice == 2) = 1;
      end
    % end
    
    % choice data for $5 only, for rationality check only
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
    
    %% Prepare variables for model fitting & fit the model. Use subjective rating for both monetary and medical.
    fixed_valueP = 5; % Objective value of fixed reward
    
    % subjective rating
    if strcmp(domain, 'MON') ==1 % Monetary block
        subjRefVal = rating(find(rating(:,1)==subjectNum),3) * ones(length(choice), 1);
        %values = Data.vals(include_indices);
        subjRatings = ones(length(values),1);
        for i=1:length(subjRatings)
            subjRatings(i) = rating(find(rating(:,1)==subjectNum),1+find(rating(1,2:6)==values(i)));
        end
    else % Medical block
        subjRefVal = rating(find(rating(:,1)==subjectNum),8) * ones(length(choice), 1);
        %values = Data.vals(include_indices);
        subjRatings = ones(length(values),1);
        for i=1:length(subjRatings)
            subjRatings(i) = rating(find(rating(:,1)==subjectNum),6+find(rating(1,7:11)==values(i)));
        end
    end
    
%     % test if different scale of subjective rating affects fitting
%     subjRefVal = subjRefVal.*10;
%     subjRatings = subjRatings.*10;
    
    fixed_prob = 1;   % prb of fixed reward 
    ambig = unique(ambigs(ambigs > 0)); % All non-zero ambiguity levels 
    prob = unique(probs); % All probability levels
    base = 0; % ? % TODO: Find out meaning -- undescribed in function. RJ-another parm in the model. Not used.

    model = 'ambigNrisk';
    b0 = [-1 .5 .5]'; % ? % TODO: Find out meaning -- undescribed in function
    refProb = fixed_prob  * ones(length(choice), 1);

    % TODO: Check that the calculation in the called function makes sense
    % Two versions of function:
    %       fit_ambgiNrisk_model: unconstrained
    %       fit_ambigNrisk_model_Constrained: constrained on alpha and beta
    [info, p] = fit_ambigNrisk_model_Constrained(choice, ...
        subjRefVal', ...
        subjRatings', ...
        refProb', ...
        probs', ...
        ambigs', ...
        model, ...
        b0, ...
        base);

    slopeP = info.b(1);
    aP = info.b(3);
    bP = info.b(2);
    r2P = info.r2;
 
    %% Create choice matrices

    % One matrix per condition. Matrix values are binary (0 for sure
    % choice, 1 for lottery). Matrix dimensions are prob/ambig-level
    % x payoff values. Used for graphing and some Excel exports.

    % Inputs: 
    %  Data
    %   .values, .ambigs, .probs, .choices (filtered by include_indices and transformed)
    %  ambig, prob (which are subsets of ambigs and probs, ran through `unique`)
    %
    % Outputs:
    %  ambigChoicesP
    %  riskyChoicesP
    %
    % Side-effects:
    %  one graph generated per-subject-domain
    %  .ambigChoicesP and .riskyChoicesP saved into `fitpar` file

    % Ambiguity levels by payoff values
    valueP = unique(values(ambigs > 0)); % each lottery payoff value under ambiguity
    ambigChoicesP = zeros(length(ambig), length(valueP)); % each row an ambiguity level
    ambigChoicesC = zeros(length(ambig), length(valueP)); % count of trials, each row an ambiguity level
    for i = 1:length(ambig)
        for j = 1:length(valueP)
            selection = find(ambigs == ambig(i) & values == valueP(j));
            if ~isempty(selection)
                ambigChoicesC(i,j) = length(selection);
                ambigChoicesP(i, j) = mean(choice(selection)); % Have already got rid of NaN choices, do not need to worry about nanmean.
            else
                ambigChoicesP(i, j) = NaN;
            end
        end
    end
    
    
    % Risk levels by payoff values
    valueP = unique(values(ambigs == 0));
    riskyChoicesP = zeros(length(prob), length(valueP));
    riskyChoicesC = zeros(length(prob), length(valueP));    
    for i = 1:length(prob)
        for j = 1:length(valueP)
            selection = find(probs == prob(i) & values == valueP(j) & ambigs == 0);
            if ~isempty(selection)
                riskyChoicesC(i, j) = length(selection);
                riskyChoicesP(i, j) = mean(choice(selection));
            else
                riskyChoicesP(i, j)=NaN;
            end
        end
    end
    
    
    %% Graph
%    colors =   [255 0 0;
%     180 0 0;
%     130 0 0;
%     52 181 233;
%     7 137 247;
%     3 85 155;
%     ]/255;
% 
%     figure
%     counter=5;
%     for i=1:3
%         subplot(3,2,counter)
%         plot(valueP,ambigChoicesP(i,:),'--*','Color',colors(3+i,:))
%         legend([num2str(ambig(i)) ' ambiguity'])
%         if counter==1
%             title(['Beta = ' num2str(bP)])
%         end
%         ylabel('Chose Lottery')
%         if counter==5
%         xlabel('Lottery Value ($)')
%         end
%         counter=counter-2;
%     end
% 
%     counter=2;
%     for i=1:3
%         subplot(3,2,counter)
%         plot(valueP,riskyChoicesP(i,:),'--*','Color',colors(i,:))
%         legend([num2str(prob(i)) ' probability'])
%         if counter==2
%             title(['Alpha = ' num2str(aP)])
%         end
%             if counter==6
%         xlabel('Lottery Value ($)')
%             end
%         counter=counter+2;
%     end
% 
%     set(gcf,'color','w');
%     figName=['RA_GAINS_' num2str(subjectNum) '_fitpar'];
%     exportfig(gcf,figName,'Format','eps','bounds','tight','color','rgb','LockAxes',1,'FontMode','scaled','FontSize',1,'Width',4,'Height',2,'Reference',gca);

    %% Save generated values
    Data.riskyChoices = riskyChoicesP; % Choice probabilaty
    Data.ambigChoices = ambigChoicesP;
    Data.ambigChoicesC = ambigChoicesC; % Trial counts
    Data.riskyChoicesC = riskyChoicesC;
    
    Data.choiceProb5 = choice_prob_5;
    
    Data.subjRatings = subjRatings;
    
        Data.MLE = info;
        Data.alpha = info.b(3);
        Data.beta = info.b(2);
        Data.gamma = info.b(1);
    
    % save data struct for the two domains
    if strcmp(domain, 'MON') ==1;
        Datamon = Data;
        clear Data
        save(fullfile(fitpar_out_path, ['MDM_' domain '_' num2str(subjectNum) '_fitpar_subjratg_constrained.mat']), 'Datamon')
    else
        Datamed = Data;
        clear Data
        save(fullfile(fitpar_out_path, ['MDM_' domain '_' num2str(subjectNum) '_fitpar_subjratg_constrained.mat']), 'Datamed')
    end

  end
end

