%% History
% Written: Ruonan 10.23.2017

function y = ambig_utility(base, v,p,AL,alpha,beta,model)
if strcmp(model,'riskAmbigPremium')
    y = v .* p - alpha .* v .^ 2 .*p .* (1 - p) - beta .* AL .^ 2; %sv = ev - risk premium on risk - risk premium on ambig
elseif strcmp(model,'riskPremium')
    y = v .* p - alpha .* v .^ 2 .*(p-beta.*AL./2) .* (1 - p + beta.*AL./2); %sv = ev - risk premium, assume a single guess of probability
elseif strcmp(model,'riskAmbigPremium2')
    y = v .* p - alpha .* p .* (1 - p) - beta .* AL .^ 2; %sv = ev - risk premium on risk - risk premium on ambig, but no quadratic term of V
elseif strcmp(model, '')
    % another model, distribution of risk? 
end

end



