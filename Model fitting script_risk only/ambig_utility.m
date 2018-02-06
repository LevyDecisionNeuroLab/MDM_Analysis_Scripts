%% History
% Ruonan 10.23.2017: Add financial utility model
% Ruonan 12.19.2017: reduced models to fit only risky trials


function y = ambig_utility(base,v,p,alpha,model)

if (strcmp(model,'ambiguity') || strcmp(model,'ambigNrisk')) || strcmp(model,'ambigNriskFixSlope')
    % the model we are using
    y = p .* v .^alpha + (1-p) .* base .^alpha;
elseif strcmp(model,'ambigPower')
    y = p .* v .^alpha; % change that
elseif strcmp(model,'discounting')
    %y = v ./ (1 + alpha.*log(1+(1-p+beta.*AL./2)./(p-beta.*AL./2)));
    y = v ./ (1 + alpha.*(1-p)./p);
    %y = v ./ (1 + alpha.*(1-p)./p);
elseif strcmp(model,'riskAmbigPremiumVar')
    y = v .* p - alpha .* v .^ 2 .*p .* (1 - p) ; %sv = ev - risk premium on risk - risk premium on ambig
elseif strcmp(model,'riskAmbigPremiumStd')
    y = v .* p - alpha .* v .* sqrt(p .* (1 - p)); %sv = ev - risk premium on risk - risk premium on ambig
elseif strcmp(model,'riskAmbigPremium')
    y = v .* p - alpha .* p .* (1 - p); %sv = ev - risk premium on risk - risk premium on ambig, but no quadratic term of V
end


