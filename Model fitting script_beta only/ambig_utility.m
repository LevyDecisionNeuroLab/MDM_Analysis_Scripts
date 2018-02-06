%% History
% Ruonan 10.23.2017: Add financial utility model


function y = ambig_utility(base,v,p,AL,beta,model)

if (strcmp(model,'ambiguity') || strcmp(model,'ambigNrisk')) || strcmp(model,'ambigNriskFixSlope')
    % the model we are using
    y = (p - beta .* (AL./2)) .* v  + (1-p - beta .* (AL./2)) .* base;
end
end


