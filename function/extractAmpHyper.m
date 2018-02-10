function V = extractAmpHyper(model,options)
%% EXTRACT HYPERPARAMETERS IN THE MODEL
switch options.inferType
    case 2
        V = [];
    case 6
        V = [model.prior.a;model.prior.b];
end
end