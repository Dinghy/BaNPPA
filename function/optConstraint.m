function [LB,UB,LBHyper,UBHyper] = optConstraint(model,options)
% OPTCONSTRAINT Summary of this function goes here
%   Detailed explanation goes here
nEps = 1e-2;
D=length(model.GP.logtheta{1})-2;
M=model.m;

LB0 = -inf(M+M*(M+1)/2,1); % Lower Bound
UB0 = inf(M+M*(M+1)/2,1);  % Upper Bound
cnt=1;
for j = 1:M % L
    LB0(M+cnt)=nEps;
    cnt=cnt+(M-j+1);
end

% store as many as the cluster
LB=[];% Lower Bound
UB=[];% Upper Bound
for c=1:options.cluster
    LB=[LB;LB0];
    UB=[UB;UB0];
end
%% Bound on Hyperparameters
LBHyper = -inf(model.K*D+1,1);
UBHyper = inf(model.K*D+1,1);
LBHyper(end) = log(1e-8*model.prior.gbase); % 1% of real data
UBHyper(end) = log(0.01*model.prior.gbase); % 1% of real data

end

