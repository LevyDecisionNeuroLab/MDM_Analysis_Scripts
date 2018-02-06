function choiceMatrix = create_choice_matrix(vals,probs,choices)

    % Chreate choice probability matrix for risky trials only.Matrix dimensions are prob-level
    % x payoff values. Used for graphing and some Excel exports.

    % Input: 
    %       vals:       values by trial
    %       probs:      probability levels by trial
    %       choices:    choices prob by trial
    %
    % Output:
    %  choiceMatrix
    %       riskProb:       riskty trials, dimension = risk level * val level
    %       riskCount:      risky trials, count of registered trial number
    %                       in case subject missed trials
    %
    % Histroy
    % ruonan 12.08.2017 written
    
    probUniq = unique(probs); % All probability levels
    value = unique(vals);
    
    % Risk levels by payoff values
    riskyChoicesProb = zeros(length(probUniq), length(value));
    riskyChoicesCount = zeros(length(probUniq), length(value));    
    for i = 1:length(probUniq)
        for j = 1:length(value)
            selection = find(probs == probUniq(i) & vals == value(j));
            if ~isempty(selection)
                riskyChoicesCount(i, j) = length(selection);
                riskyChoicesProb(i, j) = mean(choices(selection));
            else
                riskyChoicesProb(i, j)=NaN;
            end
        end
    end
    
    % output
    choiceMatrix = struct;
    choiceMatrix.riskProb = riskyChoicesProb;
    choiceMatrix.riskCount = riskyChoicesCount;
end