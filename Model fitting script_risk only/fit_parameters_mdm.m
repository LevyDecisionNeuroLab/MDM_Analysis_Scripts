% NOTE: Requires MATLAB optim library

% run this file to fit all possible models to each individual subject
% model fitting results saved in MLE structures
% subjective ratings are also saved in the *_fitpar.mat

clearvars
close all

%% Define conditions
fitparwave = 'Behavior data fitpar_1220042017'; % folder to save all the fitpar data structures
fitbyrating = false; % whether to use subjective ratings to fit
model = 'riskAmbigPremium'; % which utility function
includeAmbig = false;
search = 'grid';

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
        include_indices = Data.ambigs' == 0 & Data.choice ~= 0;
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
        prob = unique(probs); % All probability levels
        base = 0; % ? % TODO: Find out meaning -- undescribed in function. RJ-another parm in the model. Not used.

        if strcmp(search, 'grid')
            % grid search
            % range of each parameter
            if strcmp(model,'ambigNrisk')
                slopeRange = -4:0.2:1;
                aRange = 0:0.2:4;
            else
                slopeRange = -4:0.2:1;
                aRange = -2:0.2:2;
            end
            % three dimenstions
            [b1, b2] = ndgrid(slopeRange, aRange);
            % all posibile combinatinos of three parameters
            b0 = [b1(:) b2(:)];
        elseif strcmp(search,'single')
            % single search
            b0 = [-1 0.5]; % starting point of the search process, [gamma, alpha]
        elseif strcmp(search, 'random')
            % independently randomized multiple search starting points
            bstart = [-1 1]; % starting point of the search process, [gamma, alpha]
            itr = 100; % 100 iteration of starting point
            b0 = zeros(itr,length(bstart));
            for i = 1:itr
                % gamma: negative, around -1, so (-2,0)
                % alpha: (0,4)
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
            model, ...
            b0, ...
            base);

        slope = info.b(1);
        a = info.b(2);
        r2 = info.r2;
        
        % choice probability for each trial based on fitted model parameters
        % should not using the model fitting inputs, but rather also
        % include missing response trials. So IMPORTANTLY, use all trials!
        riskIndices = Data.ambigs' == 0;
        if (~fitbyrating && strcmp(domain, 'MON') ==1)
            choiceModeled = choice_prob_ambigNrisk(base,fixed_valueP * ones(sum(riskIndices),1)',Data.vals(riskIndices)',...
                fixed_prob  * ones(sum(riskIndices), 1)',Data.probs(riskIndices)',info.b,model);
        elseif fitbyrating
            choiceModeled = choice_prob_ambigNrisk(base,fixed_valueP * ones(sum(riskIndices), 1)',subjRatings(riskIndices)',...
                fixed_prob  * ones(sum(riskIndices),1)',Data.probs(riskIndices)',info.b,model);
        end
                
        % calculate subject values by unconstrained fit
        if ~fitbyrating && strcmp(domain, 'MON') ==1 
            sv = ambig_utility(0, ...
              Data.vals(riskIndices)', ...
              Data.probs(riskIndices)', ...
              a, ...
              model);
        elseif fitbyrating
            sv = ambig_utility(0, ...
              subjRatings(riskIndices)', ...
              Data.probs(riskIndices)', ...
              a, ...
              model);
        end
        
        svRef = ambig_utility(0, ...
              fixed_valueP, ...
              fixed_prob, ...
              a, ...
              model);
    end   
       
    %% Create choice matrices
    % One matrix per condition. Matrix values are binary (0 for sure
    % choice, 1 for lottery). Matrix dimensions are prob/ambig-level
    % x payoff values. Used for graphing and some Excel exports.
 
    choiceMatrix = create_choice_matrix(values,probs,choice);
    
    %% Create matrix for subjective value
    valueP = unique(values);

    if (~fitbyrating && strcmp(domain, 'MON') ==1) || fitbyrating
        svByLott = zeros(length(prob), length(valueP));
        for i = 1:length(prob)
            for j = 1:length(valueP)
               svByLott(i,j) = ambig_utility(0,valueP(j),prob(i),a,model); 
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
    Data.subjRatings = subjRatings;
    Data.choiceModeled = choiceModeled;

    
    
    if (~fitbyrating && strcmp(domain, 'MON') ==1) || fitbyrating
        Data.MLE = info;
        Data.alpha = info.b(2);
        Data.gamma = info.b(1);
        Data.r2 = info.r2;
        Data.sv = sv;
        Data.svRef = svRef;
        Data.svByLott = svByLott;
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

