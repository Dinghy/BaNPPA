function [f,df] = varsGpPpBoundAmp(V,model,modelD,options,user)
%VARSGPPPCONBOUNDE Summary of this function goes here
%   Detailed explanation goes here
model = returnAmp(V,model,options,user);
% store result
modelD.f=0;
modelD.K = model.K;
switch options.inferType
    case 6
        modelD.amp.a = zeros(size(model.amp.a));
        modelD.amp.b = zeros(size(model.amp.b));
    case 2
        modelD.amp = zeros(size(model.amp));
end
modelD = varsTransform(model,modelD,options);

% calculation
u = user;
modelD = calPointValue(model,modelD,u,3,options);

% 2 Integration over Tlim--SEard case(1 dimension test)
switch options.inferType
    case 2
        modelD.f=modelD.f-sum(modelD.mAmp_2(:,u).*modelD.vfInt);
        modelD.amp(:,u)=modelD.amp(:,u)-modelD.vfInt;
end

%% return variables
f = - modelD.f;
df = -extractAmp(modelD,options,u);
end
