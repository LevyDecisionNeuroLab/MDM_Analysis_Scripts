% NOTE: Requires MATLAB optim library

% run this file to fit all possible models to each individual subject
% model fitting results saved in MLE structures
% subjective ratings are also saved in the *_fitpar.mat

clearvars
close all

%% Define conditions
fitparwave = 'Behavior data fitpar_1226012017'; % folder to save all the fitpar data structures
fitbyrating = true; % whether to use subjective ratings to fit
model = 'ambigNrisk'; % which utility function
includeAmbig = true;
search = 'single';

%% Set up loading & subject selection
root = 'D:\Ruonan\Projects in the lab\MDM Project\Medical Decision Making Imaging\MDM_imaging\Behavioral Analysis';
data_path = fullfile(root, 'PTB Behavior Log/'); % root of folders is sufficient
rating_filename = fullfile(root, 'Behavior Analysis/MDM_Rating.csv');
fitpar_out_path = fullfile(root, 'Behavior fitpar files',fitparwave);

% if folder does not exist, create folder
if exist(fitpar_out_path)==0
    mkdir(fullfile(root, 'Behavior fitpar files'),fitparwave)
end

addpath(genpath(data_path)); % generate path for all the subject data folder

% get all subjects number in the folder
subjects = getSubjectsInDir(data_path, 'subj');
exclude = [2581]; % TEMPORARY: subjects incomplete data (that the script is not ready for)
subjects = subjects(~ismember(subjects, exclude));
% subjects = [2585];

% load subjective ratings
% column1-subj ID, c2-$0, c3-$5,c4-$8,c5-$12,c6-$25,c7-no effect, c8-slight,c9-moderate,c10-major,c11-recovery.
rating = csvread(rating_filename,1,0); %reads data from the file starting at row offset R1 and column offset C1. For example, the offsets R1=0, C1=0 specify the first value in the file.

%% Individual subject fitting

for subj_idx = 1:length(subjects)
  domains = {'MON', 'MED'};

  for domain_idx = 1:length(domains)
    subjectNum = subjects(subj_idx);
    domain = domains{domain_idx};
    
    fname = sprintf('MDM_%s_%d.mat', domain, subjectNum);
    load(fname) % produces variable `Datamon` or 'Datamed' for convenience, change its name into 'Data'
%     sprintf('MDM_%s_%d.mat', domain, subjectNum);
    
    if strcmp(domain, 'MON') ==1
        Data = Datamon;
        clear Datamon
    else
        Data = Datamed;
        clear Datamed
    end
    
    %% Load subjective ratings
    % prepare subjective rating for each trial
    if strcmp(domain, 'MON') ==1 % Monetary block
        subjRefRatings = rating(find(rating(:,1)==subjectNum),3) * ones(length(Data.choice), 1);
        %values = Data.vals(include_indices);
        subjRatings = ones(length(Data.vals),1);
        for i=1:length(subjRatings)
            subjRatings(i) = rating(find(rating(:,1)==subjectNum),1+find(rating(1,2:6)==Data.vals(i)));
        end
    else % Medical block
        subjRefRatings = rating(find(rating(:,1)==subjectNum),8) * ones(length(Data.choice), 1);
        %values = Data.vals(include_indices);
        subjRatings = ones(length(Data.vals),1);
        for i=1:length(subjRatings)
            subjRatings(i) = rating(find(rating(:,1)==subjectNum),6+find(rating(1,7:11)==Data.vals(i)));
        end
    end
    
    %% Refine variables
    
    if includeAmbig
        % Exclude non-responses
        include_indices = Data.choice ~= 0;
    else
        % Exclude ambiguious trials (fit only risky trials)
        include_indices = Data.ambigs' ~= 0 & Data.choice ~= 0;
    end

    choice = Data.choice(include_indices);
    values = Data.vals(include_indices);
    ambigs = Data.ambigs(include_indices);
    probs  = Data.probs(include_indices);
    ratings = subjRatings(include_indices);
    refRatings = subjRefRatings(include_indices);
    
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
    
    %% Fitting, whether to use outcome magnitude or rating is conditioned on the variable 'fitbyrating' 
    
    % define fitting values if fitting by rating
    if fitbyrating     
        fixed_valueP = refRatings(1);
        fitrefVal = refRatings;
        fitVal = ratings;
    end
    
    % define fitting values if fit by objective value in the monetary domain
    if ~fitbyrating && strcmp(domain, 'MON') ==1 
        fixed_valueP = 5; % Value of fixed reward
        fitrefVal = fixed_valueP * ones(length(choice), 1);
        fitVal = values;
    end
    
    % fit the model
    if (~fitbyrating && strcmp(domain, 'MON') ==1) || fitbyrating
        fixed_prob = 1;   % prb of fixed reward 
        refProb = fixed_prob  * ones(length(choice), 1);
        fixed_ambig = 0;
        ambig = unique(ambigs(ambigs > 0)); % All non-zero ambiguity levels 
        prob = unique(probs); % All probability levels
        base = 0; % ? % TODO: Find out meaning -- undescribed in function. RJ-another parm in the model. Not used.

        if strcmp(search, 'grid')
            % grid search
            % range of each parameter
            if strcmp(model,'ambigNrisk')
                slopeRange = -4:0.2:1;
                bRange = -2:0.2:2;
            end
            % three dimenstions
            [b1, b2] = ndgrid(slopeRange, bRange);
            % all posibile combinatinos of three parameters
            b0 = [b1(:) b2(:)];
        elseif strcmp(search,'single')
            % single search
            b0 = [-1 0.5]; % starting point of the search process, [gamma, beta, alpha]
        elseif strcmp(search, 'random')
            % independently randomized multiple search starting points
            bstart = [-1 0 1]; % starting point of the search process, [gamma, beta, alpha]
            itr = 100; % 100 iteration of starting point
            b0 = zeros(itr,length(bstart));
            for i = 1:itr
                % gamma: negative, around -1, so (-2,0)
                % beta: [-1,1] possible to be larger than 1?
                % alpha: restricted to 1
                b0(i,:) = bstart + [-1+2*rand(1) -1+2*rand(1)]; % randomize search starting point, slope, beta, alpha
            end
        end
        

        % Two versions of function:
        %       fit_ambgiNrisk_model: unconstrained
        %       fit_ambigNrisk_model_Constrained: constrained on alpha and beta
        
        % Unconstrained fitting
        [info, p] = fit_ambigNrisk_model(choice, ...
            fitrefVal', ...
            fitVal', ...
            refProb', ...
            probs', ...
            ambigs', ...
            model, ...
            b0, ...
            base);

        slope = info.b(1);
        b = info.b(2);
        r2 = info.r2;
        
        % choice probability for each trial based on fitted model parameters
        % should not using the model fitting inputs, but rather also
        % include missing response trials. So IMPORTANTLY, use all trials!
        if (~fitbyrating && strcmp(domain, 'MON') ==1)
            choiceModeled = choice_prob_ambigNrisk(base,fixed_valueP * ones(length(Data.vals), 1)',Data.vals',...
                fixed_prob  * ones(length(Data.vals), 1)',Data.probs',Data.ambigs',info.b,model);
        elseif fitbyrating
            choiceModeled = choice_prob_ambigNrisk(base,fixed_valueP * ones(length(subjRatings), 1)',subjRatings',...
                fixed_prob  * ones(length(Data.vals), 1)',Data.probs',Data.ambigs',info.b,model);
        end
            
        
        
        % calculate subject values by unconstrained fit
        if ~fitbyrating && strcmp(domain, 'MON') ==1 
            for reps = 1:length(Data.choice)
              sv(reps, 1) = ambig_utility(0, ...
                  Data.vals(reps), ...
                  Data.probs(reps), ...
                  Data.ambigs(reps), ...
                  b, ...
                  model);
            end
        elseif fitbyrating
           for reps = 1:length(Data.choice)
              sv(reps, 1) = ambig_utility(0, ...
                  subjRatings(reps), ...
                  Data.probs(reps), ...
                  Data.ambigs(reps), ...
                  b, ...
                  model);
            end            
        end
        
        svRef = ambig_utility(0, ...
              fixed_valueP, ...
              fixed_prob, ...
              fixed_ambig, ...
              b, ...
              model);
    end   
       
    %% Chosen SV (cv), chosen reward magnitude (cr), chosen subective rating (CRating)
    
    % All choices
    choiceAll = Data.choice;
    valuesAll = Data.vals;
    refValue = 5;
    ambigsAll = Data.ambigs;
    probsAll  = Data.probs;
    ratingsAll = subjRatings;
    refRatingsAll = subjRefRatings;
    % mark miss-response
    choiceAll(choiceAll==0) = NaN;

    % Side with lottery is counterbalanced across subjects 
    % code 0 as reference choice, 1 as lottery choice
    % if sum(choice == 2) > 0 % Only if choice has not been recoded yet. RJ-Not necessary
    % RJ-If subject do not press 2 at all, the above if condition is problematic
      if Data.refSide == 2
          choiceAll(choiceAll == 2) = 0;
          choiceAll(choiceAll == 1) = 1;
      elseif Data.refSide == 1 % Careful: rerunning this part will make all choices 0
          choiceAll(choiceAll == 1) = 0;
          choiceAll(choiceAll == 2) = 1;
      end
    
    
      if (~fitbyrating && strcmp(domain, 'MON') ==1) || fitbyrating
          %chosen subjective value of all trials
          cv = zeros(length(choiceAll),1);
          cv(cv == 0) = NaN; % this makes sure the missed-trial is marked by NaN
          cv(choiceAll == 0) = svRef;
          cv(choiceAll ==1 ) = sv(choiceAll == 1);
      end
      
      % chosen reward magnitude for both monetary and medical
      cr = zeros(length(choiceAll),1);
      cr(cr == 0) = NaN;
      cr(choiceAll ==0) = refValue;
      cr(choiceAll ==1) = valuesAll(choiceAll ==1 );

      % chosen subjective rating for both monetary and medical
      cRating = zeros(length(choiceAll),1);
      cRating(cRating == 0) = NaN;
      cRating(choiceAll ==0) = refRatingsAll(1);
      cRating(choiceAll ==1) = ratingsAll(choiceAll ==1 );
    
    %% Create choice matrices
    % One matrix per condition. Matrix values are binary (0 for sure
    % choice, 1 for lottery). Matrix dimensions are prob/ambig-level
    % x payoff values. Used for graphing and some Excel exports.
 
    choiceMatrix = create_choice_matrix(values,ambigs,probs,choice);
    
    %% Create matrix for subjective value
    valueP = unique(values(ambigs == 0));

    if (~fitbyrating && strcmp(domain, 'MON') ==1) || fitbyrating
        svByLott = zeros(length(prob)+length(ambig), length(valueP));
        for i = 1:length(prob)+length(ambig)
            for j = 1:length(valueP)
                if i < length(prob)+1
                   svByLott(i,j) = ambig_utility(0, ...
                                               valueP(j), ...
                                               prob(i), ...
                                               0, ...
                                               b, ...
                                               model); 
                else
                   svByLott(i,j) = ambig_utility(0, ...
                                               valueP(j), ...
                                               0.5, ...
                                               ambig(i-length(prob)), ...
                                               b, ...
                                               model);  
                end
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
%             title(['Beta = ' num2str(b_uncstr)])
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
%             title(['Alpha = ' num2str(a_uncstr)])
%         end
%             if counter==6
%         xlabel('Lottery Value ($)')
%             end
%         counter=counter+2;
%     end
% 
%     set(gcf,'color','w');
%     figName=['RA_GAINS_' num2str(subjectNum) '_fitpar'];
% %     exportfig(gcf,figName,'Format','eps','bounds','tight','color','rgb','LockAxes',1,'FontMode','scaled','FontSize',1,'Width',4,'Height',2,'Reference',gca);


%% graph with fitted lines
% 
%     xP = 0:0.1:max(valueP);
%     uFP = fixed_prob * (fixed_valueP).^a_uncstr;
%      
%    figure
%      
%     % risk pos
%     for i = 1 :length(prob)
%         plot(valueP,riskyChoicesP(i,:),'o','MarkerSize',8,'MarkerEdgeColor',colors([1 1 1])...
%             ,'MarkerFaceColor',colors(i,:),'Color',colors(i,:));
%           hold on
%         % logistic function
%         uA = prob(i) * xP.^a_uncstr;
%         p = 1 ./ (1 + exp(slope_uncstr*(uA-uFP)));
% 
%         plot(xP,p,'-','LineWidth',4,'Color',colors(i,:));
%         axis([0 25 0 1])
%         set(gca, 'ytick', [0 0.5 1])
%         set(gca,'xtick', [0 5 10 15 20 25])
%         set(gca,'FontSize',25)
%         set(gca,'LineWidth',3)
%         set(gca, 'Box','off')
% 
% 
%     end
% %     title(['  alpha gain = ' num2str(a_uncstr)]);
%     
%     figure
%     % ambig pos
%     for i = 1:length(ambig)
%         plot(valueP,ambigChoicesP(i,:),'o','MarkerSize',8,'MarkerEdgeColor',colors([1 1 1]),'MarkerFaceColor',colors(length(prob)+i,:));
%          hold on
% % 
%         % logistic function
%         uA = (0.5 - b_uncstr.*ambig(i)./2) * xP.^a_uncstr;
%         p = 1 ./ (1 + exp(slope_uncstr*(uA-uFP)));
% 
% 
%         plot(xP,p,'-','LineWidth',2,'Color',colors(length(prob)+i,:));
%         axis([0 25 0 1])
%         set(gca, 'ytick', [0 0.5 1])
%         set(gca,'xtick', [0 5 10 15 20 25])
%         set(gca,'FontSize',25)
%         set(gca,'LineWidth',3)
%         set(gca, 'Box','off')
% 
%     end
% %     title([ '  beta gain = ' num2str(b_uncstr)]);

    %% Save generated values
    Data.choiceMatrix = choiceMatrix;
    Data.choiceProb5 = choice_prob_5;
    Data.subjRatings = subjRatings;
    
    % choices per each trial, 0-ref,1-lottery
    Data.choiceLott = choiceAll;
    Data.choiceModeled = choiceModeled;

    
    % chosen reward magnitude and chosen subjective ratins
    Data.chosenVal = cr;
    Data.chosenRating = cRating;
    
    
    if (~fitbyrating && strcmp(domain, 'MON') ==1) || fitbyrating
        Data.MLE = info;
        Data.beta = info.b(2);
        Data.gamma = info.b(1);
        Data.r2 = info.r2;
        Data.sv = sv;
        Data.svRef = svRef;
        Data.svByLott = svByLott;
        Data.svChosen = cv;
    end
    

    % save data struct for the two domains
    if strcmp(domain, 'MON') ==1
        Datamon = Data;
        clear Data
        save(fullfile(fitpar_out_path, ['MDM_' domain '_' num2str(subjectNum) '_fitpar.mat']), 'Datamon')
    else
        Datamed = Data;
        clear Data
        save(fullfile(fitpar_out_path, ['MDM_' domain '_' num2str(subjectNum) '_fitpar.mat']), 'Datamed')
    end

  end
end

