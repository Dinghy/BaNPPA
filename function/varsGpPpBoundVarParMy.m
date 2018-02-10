function [f,df] = varsGpPpBoundVarParMy(V,model,modelD,options)
% VARSGPPPCONBOUNDE Summary of this function goes here
% Detailed explanation goes here
model = returnVariationParamsMy(model,V,options);

% store result
modelD.f=0;

% initialization to store gradient
modelD.U=model.U;
modelD.K=model.K;
modelD.var.Mu=cell(1,model.K);
modelD.var.L=cell(1,model.K);
modelD.diagKnn = cell(1,model.K);
% parametrization for amplitude

for c=1:model.K
    modelD.var.Mu{c}=zeros(model.m,1);
    modelD.var.L{c}=zeros(model.m,model.m);
    modelD.diagKnn{c} = sum(exp(2*model.GP.logtheta{c}(end-1:end)));
    modelD.InvKmmMu{c} = modelD.Kmm{c}\model.var.Mu{c};
    modelD.InvKmmSigma{c} = modelD.Kmm{c}\model.var.Sigma{c};
    modelD.InvKmmL{c} = modelD.Kmm{c}\model.var.L{c};
end

for u=1:model.U
    % 1 point value and derivatives
    modelD = calPointValue(model,modelD,u,1,options);
end
for c=1:model.K
    modelD.var.L{c}=modelD.var.L{c}*model.var.L{c};
end
% 2 Integration over Tlim--SEard case(1 dimension test)
modelD = calIntegral(model,options,modelD,1);

% 3 KL divergence and gradient
modelD = calDivergence(model,modelD,1);

% Get derivatives in vector
f = -modelD.f/model.U;
df = -extractVariationParamsMy(modelD,options)/model.U;

end
