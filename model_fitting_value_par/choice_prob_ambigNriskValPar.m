function p = choice_prob_ambigNriskValPar(base,vF,vA,pF,pA,AL,beta,model,vals)
% CHOICE_PROB_AMBIGNRISK                Binary logit choice probability
%
%     INPUTS
%     vF    - values of fixed option
%     vA    - value of the ambiguous option
%     pF    - probability of fixed option
%     pA    - probability of ambiguous option
%     AL    - ambiguity level
%     beta  - Parameters corresponding to MODEL: slope, beta, alpha, val1,
%     val2, val3, val4, with 0 as the reference ($0)
%     model - String indicating which model to fit; currently valid are:
%               'ambigNrisk' - (p-beta(2)*AL/2)*v^beta(3*)
%     vals   - all values in the experiment paradigm
%
%     OUTPUTS
%     p     - choice probabilities for the *SHORTER* option
%
%     REVISION HISTORY:
%     brian 03.14.06 written
%     ifat  12.01.06 adapted for ambiguity and risk
%     ruonan 8.19.19 include values as parameters
    
if (strcmp(model,'ambigNriskValPar')) 
    % get the values as parameters
    vals_par = beta(4:7); % beta(4:7) and val dimension must agree

    % select value
    % get the index of vF in val
    [~, vF_par_idx] = ismember(vF, vals);
    % get the value parameter
    vF_par = vals_par(vF_par_idx)'; % change dimeinsion into 1 by n
    % calculate sv
    uF = ambig_utility(base,vF_par,pF,zeros(size(vF)),beta(3),beta(2),model); %fixed non-ambiguous

    % get the index of vA in val
    [~, vA_par_idx] = ismember(vA, vals);
    % get the value parameter
    vA_par = vals_par(vA_par_idx)';
    uA = ambig_utility(base,vA_par,pA,AL,beta(3),beta(2),model);
    
    slope = beta(1);
elseif strcmp(model,'ambigSVPar')
    % this model does not have alpha as parameter
    alpha = 0; % no alpha, but has to pass to the function
    
    % [slope, beta, val1, val2, val3, val4]
    % get the values as parameters
    vals_par = beta(3:6);

    % select value
    % get the index of vF in val
    [~, vF_par_idx] = ismember(vF, vals);
    % get the value parameter
    vF_par = vals_par(vF_par_idx)'; % change dimeinsion into 1 by n
    
    % calculate sv
    uF = ambig_utility(base,vF_par,pF,zeros(size(vF)),alpha,beta(2),model); %fixed non-ambiguous

    % get the index of vA in val
    [~, vA_par_idx] = ismember(vA, vals);
    % get the value parameter
    vA_par = vals_par(vA_par_idx)';
    uA = ambig_utility(base,vA_par,pA,AL,alpha,beta(2),model); 
    
    slope = beta(1);  
end
%s = ones(size(uA)); %sign(uA);
p = 1 ./ (1 + exp(slope*(uA-uF)));

return

