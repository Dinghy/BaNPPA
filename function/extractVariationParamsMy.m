function [V] = extractVariationParamsMy(model,options)
%EXTRACTVARIATIONPARAMSMY Summary of this function goes here
% 1, single model; 2, mixture model; 3, mixture+supervised part
V=[];
for c=1:model.K
    V=[V;model.var.Mu{c};model.var.L{c}(tril(true(size(model.var.L{c}))))];
end
end

