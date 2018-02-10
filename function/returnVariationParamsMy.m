function model = returnVariationParamsMy(model,V,options)
% UNTITLED Summary of this function goes here
% 1, single model; 2, mixture model; 3, mixture+supervised part
U = model.U;
M = model.m;
M2= M*(M+1)/2;

nStart = 1;
for c=1:model.K
    model.var.Mu{c} = V(nStart:nStart+M-1);
    mtmp = tril(ones(M));
    nStart=nStart+M;
    mtmp(mtmp==1)=V(nStart:nStart+M2-1);
    model.var.L{c} = mtmp;
    model.var.Sigma{c} = mtmp*mtmp';
    nStart=nStart+M2;
end

end

