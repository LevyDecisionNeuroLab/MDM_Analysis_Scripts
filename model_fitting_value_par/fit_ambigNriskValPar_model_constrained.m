function [info,p] = fit_ambigNriskValPar_model_constrained(choice,vF,vA,pF,pA,AL,model,b0,base,vals)
% FIT_AMBIGNRISK_MODEL      Fit a variety of probabilistic ambiguity models
% 
%
%     Fits a binary logit model by maximum likelihood.
%
%     INPUTS
%     choice      - Dependent variable. The data should be *ungrouped*,
%                   such that CHOICE is a column of 0s and 1s, where 1 indicates 
%                   a choice of the FIXED option.
%     vF          - value of fixed lottery
%     vA          - value of ambiguous lottery
%     pF          - probability of fixed lottery
%     pA          - probability of ambiguous lottery
%     AL          - ambiguity level
%     b0          - search starting point, size: search iteration x
%                   parameter numbers (slope, beta, alpha, val1, val2, val3, val4)
%     model       - String indicating which model to fit; currently valid are:
%                   'ambigNriskValPar'   - power with subjective probability, estimates both risk and ambiguity coefficients, vals as parameters   
%     vals        - all vals by experiment design
%                    
%                   Multiple models can be fit by passing in a cell array
%                   of strings. 
%     
%
%     OUTPUTS
%     info       - data structure with following fields:
%                     .nobs      - number of observations
%                     .nb        - number of parameters
%                     .optimizer - function minimizer used
%                     .exitflag  - see FMINSEARCH
%                     .b         - fitted parameters; note that for all the
%                                  available models, the first element of B
%                                  is a noise term for the logistic
%                                  function, the remaining elements are
%                                  parameters for the selected discount
%                                  functions. eg., for model='exp', B(2) is
%                                  the time constant of the exponential
%                                  decay.
%                     .LL        - log-likelihood evaluated at maximum
%                     .LL0       - restricted (minimal model) log-likelihood
%                     .AIC       - Akaike's Information Criterion 
%                     .BIC       - Schwartz's Bayesian Information Criterion 
%                     .r2        - pseudo r-squared
%                   This is a struct array if multiple models are fit.
%     p           - Estimated choice probabilities evaluated at the values
%                   delays specified by the inputs vS, vR, dS, dL. This is
%                   a cell array if multiple models are fit.
%
%     EXAMPLES
%     see TEST_FAKE_DATA_AMBIGUITTY, TEST_FAKE_DATA, TEST_JOE_DATA, and TEST_KENWAY_DATA

%
%     REVISION HISTORY:
%     brian 03.10.06 written
%     brian 03.14.06 added fallback to FMINSEARCH, multiple fit capability
%     ifat  12.01.06 adapted for ambiguity and risk + CI
%     ruonan 12.04.17 change single point search to grid search
%     ruonan 8.19.19 add values as parameters to fit and set constraints



% If multiple model fits requested, loop and pack 
% if iscell(model)
%    for i = 1:length(model)
%       [info(i),p{i}] = fit_ambigNrisk_model(choice,vF,vA,pF,pA,AL,model{i},b0,base,vals);
%    end
%    return;
% end

thresh = 0.05;
nobs = length(choice);

for i = 1 : size(b0,1)
        
    b00 = b0(i,:)'; % search starting point
   
    optimizer = 'fmincon';
        
    OPTIONS = optimset('Display','off','TolCon',1e-6,'TolFun',1e-5,'TolX',1e-5,...
        'DiffMinChange',1e-4,'Maxiter',100000,'MaxFunEvals',20000);
    
    % constraint: 0<val1<val2<val3<val4
    % ineq_mat*b00 <= ineq
    % linear inequality constraints matrix, matrix different for different
    % models
    if strcmp(model, 'ambigNriskValPar')
        ineq_mat = [0 0 0 -1 0 0 0; 0 0 0 1 -1 0 0; 0 0 0 0 1 -1 0; 0 0 0 0 0 1 -1];
    elseif strcmp(model, 'ambigSVPar')
        ineq_mat = [0 0 -1 0 0 0; 0 0 1 -1 0 0; 0 0 0 1 -1 0; 0 0 0 0 1 -1];
    end
    
    % linear inequality constraints
    ineq = [0 0 0 0];
    % lower bounds
    lb = [];
    % upper bounds
    ub = [];
    
    [b,negLL,exitflag,convg] = fmincon(@local_negLL,b00,ineq_mat,ineq,[],[],lb,ub,[],OPTIONS,choice,vF,vA,pF,pA,AL,model,base,vals);
    %X = fmincon(fun,x0,A,B,Aeq,Beq,lb,ub,nonlcon,Opts,a,b)
    
    % using bads
    nonbcon = 
    [b,negLL,exitflag,convg] = bads(@local_negLL,b00,OPTIONS,choice,vF,vA,pF,pA,AL,model,base,vals, nonbcon);
    
    % Unrestricted log-likelihood
    LL = -negLL;
    if i == 1
        info.LL = LL;
    end

    if i == 1 || (i ~=1 && LL > info.LL)% first iteration; and if a later iteration renders larger likelihood, replace info
        % Choice probabilities (for VARIED)
        p = choice_prob_ambigNriskValPar(base,vF,vA,pF,pA,AL,b,model,vals);

        % Restricted log-likelihood
        LL0 = sum((choice==1).*log(0.5) + (1 - (choice==1)).*log(0.5)); % assuming no predictors, the chance of choosing and not choosing the lottery are both 50%

        % Confidence interval, requires Hessian from FMINUNC
        try
            invH = inv(-H);
            se = sqrt(diag(-invH));
        catch
        end

        info.nobs = nobs;
        info.nb = length(b);
        info.model = model;
        info.optimizer = optimizer;
        info.exitflag = exitflag;
        info.b = b;

        try
            info.se = se;
            info.ci = [b-se*norminv(1-thresh/2) b+se*norminv(1-thresh/2)]; % Wald confidence
            info.tstat = b./se;
        catch
        end

        info.LL = LL;
        info.LL0 = LL0;
        info.AIC = -2*LL + 2*length(b);
        info.BIC = -2*LL + length(b)*log(nobs);
        % https://stats.idre.ucla.edu/other/mult-pkg/faq/general/faq-what-are-pseudo-r-squareds/
        info.r2 = 1 - LL/LL0; % McFadden's Pseudo r squared = 1-LLmodel/LLwithoutModel (LL is negative)
        info.r2_adj = 1-((LL - length(b))/LL0); % adjusted r2 
    end    

end




%----- LOCAL FUNCTIONS
% This is the function to minimize, sum of -log-likelihood.
function sumerr = local_negLL(beta,choice,vF,vA,pF,pA,AL,model,base,vals)
% estimated likelihood of chooseing the lottery
p = choice_prob_ambigNriskValPar(base,vF,vA,pF,pA,AL,beta,model,vals); 

% Trap log(0)
ind = p == 1;
p(ind) = 0.9999;
ind = p == 0;
p(ind) = 0.0001;

% Log-likelihood
% If lottery is chosen, err=log(p), log likelihood of choosing the lottery
% If reference is chosen, err=log(1-p), log likelihood of choosing the reference
% Because 0<p<1, log(p) and log(1-p) < 0. sum of these likelihoods is negative 
% The sum of these likelihood should be maximized
err = (choice==1).*log(p) + (1 - (choice==1)).*log(1-p);

% Sum of -log-likelihood. sumerr is a positive value, should me minimized
sumerr = -sum(err);



