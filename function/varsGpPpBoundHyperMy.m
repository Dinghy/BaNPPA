function [f,df] = varsGpPpBoundHyperMy(V,model,options)
% VARSGPPPBOUNDHYPER Summary of this function goes here
%   Detailed explanation goes here
debug=0;
model = returnHyperParamsMy(model,options,V);


% parameters f,df
modelD.K=model.K;
modelD.GP.logtheta=cell(1,model.K);
modelD.f=0;
for c=1:model.K
    modelD.GP.logtheta{c}=zeros(size(model.GP.logtheta{c}));
end
% parametrization for amplitude
modelD = varsTransform(model,modelD,options);
% other useful parameters
modelD.numPar = length(model.GP.logtheta{1});
modelD.kmmDx=cell(1,model.K);
modelD.diagKnn=cell(1,model.K);
modelD.InvKmmMu = cell(1,model.K);
modelD.InvKmmSigma = cell(1,model.K);
modelD.InvKmmSigmaInvKmm = cell(1,model.K);
modelD.InvKmmL = cell(1,model.K);
modelD.Kmm = cell(1,model.K);
for c= 1:model.K
    modelD.Kmm{c} = feval(model.cov{:}, model.GP.logtheta{c},model.Xm);
    modelD.InvKmmMu{c} = modelD.Kmm{c}\model.var.Mu{c};
    modelD.InvKmmSigma{c} = modelD.Kmm{c}\model.var.Sigma{c};
    modelD.InvKmmSigmaInvKmm{c} = modelD.InvKmmSigma{c}/modelD.Kmm{c};
    modelD.InvKmmL{c} = modelD.Kmm{c}\model.var.L{c};
end
for c=1:model.K
    modelD.diagKnn{c} = sum(exp(2*model.GP.logtheta{c}(end-1:end)));
    modelD.kmmDx{c} = cell(modelD.numPar,1);
    for k=1:modelD.numPar
        modelD.kmmDx{c}{k} = feval(model.cov{:}, model.GP.logtheta{c}, model.Xm, [], k);
    end
end

for u=1:model.U
    %% pre-computation
     modelD = calPointValue(model,modelD,u,2,options);
end
%% 2 Integration over Tlim--SEard case(1 dimension test) Diff:50
modelD = calIntegral(model,options,modelD,2);
%% 3 KL divergence and gradient Diff:3
modelD = calDivergence(model,modelD,2);
f=-modelD.f/model.U;
modelD.derivative = 1;
df = -extractHyperParamsMy(modelD,options)/model.U;

end


