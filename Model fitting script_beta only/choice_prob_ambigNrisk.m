% CHOICE_PROB_AMBIGNRISK                Binary logit choice probability
% 
%     p = choice_prob_ambigNrisk(vF,vA,pF,pA,AL,beta,model);
%
%     INPUTS
%     vF    - values of fixed option
%     vA    - values of ambiguous option
%     pF    - probability of fixed option
%     pA    - probability of ambiguous option
%     AL    - ambiguity level
%     beta  - Parameters corresponding to MODEL
%     model - String indicating which model to fit; currently valid are:
%               'ambigNrisk' - (p-beta(2)*AL/2)*v^beta(3*)
%
%     OUTPUTS
%     p     - choice probabilities for the *SHORTER* option
%

%
%     REVISION HISTORY:
%     brian 03.14.06 written
%     ifat  12.01.06 adapted for ambiguity and risk

function p = choice_prob_ambigNrisk(base,vF,vA,pF,pA,AL,beta,model)

    % single logistic
if (strcmp(model,'ambigNrisk') || strcmp(model,'ambiguity') || strcmp(model,'ambigPower') || strcmp(model,'ambigSubjRate')) || strcmp(model,'riskAmbigPremiumVar') || strcmp(model, 'riskAmbigPremiumStd') || strcmp(model,'riskPremium') || strcmp(model,'riskAmbigPremium')
    % this is the one we are currently using!
    uF = ambig_utility(base,vF,pF,zeros(size(vF)),beta(2),model); %fixed non-ambiguous
    uA = ambig_utility(base,vA,pA,AL,beta(2),model); % ambiguous
    slope = beta(1);
end
    %s = ones(size(uA)); %sign(uA);
    p = 1 ./ (1 + exp(slope*(uA-uF)));

return

