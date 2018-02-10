function [f,df] = varsGpPpAmpHyper(V,model,modelD,options)
%% GIVEN PRE COMPUTATION OPTIMIZE THE HYPER-PARAMETERS IN GAMMA PROCESS
model = returnAmpHyper(V,model,options);
modelD = varsTransform(model,modelD,options);
modelD.f = 0;
switch options.inferType
    case 6
        modelD.prior.a = 0;
        modelD.prior.b = 0;
        for u = 1:model.U
            modelD.f = modelD.f+model.prior.a*log(model.prior.b)-gammaln(model.prior.a)...
                -model.prior.b*model.amp.s(u)+(model.prior.a-1)*log(model.amp.s(u));
            modelD.prior.a = modelD.prior.a + log(model.prior.b)-psi(model.prior.a)+log(model.amp.s(u));
            modelD.prior.b = modelD.prior.b + model.prior.a/model.prior.b-model.amp.s(u);
        end  
end
f = - modelD.f/model.U/model.K;
df = - extractAmpHyper(modelD,options)/model.U/model.K;
end