function modelD = calPre(model)
%  PRE CALCULATION FOR AMPLITUDE OPTIMIZATION AND AMPLITUDE HYPERPARAMETERS 
%  OPTIMIZATION
%  Kmm, Knm, inv(Kmm)*Mu...
modelD.Kmm = cell(1,model.K);
modelD.userKnm = cell(model.U,model.K);
modelD.userKnmInvKmm = cell(model.U,model.K);
modelD.diagKnn = cell(1,model.K);
for c= 1:model.K
    modelD.diagKnn{c} = sum(exp(2*model.GP.logtheta{c}(end-1:end)));
    modelD.Kmm{c} = feval(model.cov{:}, model.GP.logtheta{c},model.Xm);
    modelD.InvKmmMu{c} = modelD.Kmm{c}\model.var.Mu{c};
    modelD.InvKmmSigma{c} = modelD.Kmm{c}\model.var.Sigma{c};
    modelD.InvKmmL{c} = modelD.Kmm{c}\model.var.L{c};
end
modelD.mTmp_basef = cell(1,model.U);
modelD.mTmp_baseLogf = cell(1,model.U);
modelD.vVar = cell(1,model.U);
for u=1:model.U
    nData = size(model.data{u}.feature,1);
    vVar = zeros(model.K,nData);
    vMean2DVar = zeros(model.K,nData);
    for c=1:model.K
        if iscell(model.cov{2}{1});
            modelD.userKnm{u,c} = feval(model.cov{2}{1}{:},model.GP.logtheta{c}(1:end-1),model.data{u}.feature,model.Xm);
        else
            modelD.userKnm{u,c} = feval(model.cov{2}{1},model.GP.logtheta{c}(1:end-1),model.data{u}.feature,model.Xm);
        end
        modelD.userKnmInvKmm{u,c} = modelD.userKnm{u,c}/modelD.Kmm{c};
        vMean = (modelD.userKnm{u,c}*modelD.InvKmmMu{c})';
        mtmp = modelD.userKnm{u,c}*modelD.InvKmmL{c};
        vVar(c,:) = modelD.diagKnn{c}-sum(modelD.userKnmInvKmm{u,c}.*modelD.userKnm{u,c},2)+sum(mtmp.^2,2);
        vMean2DVar(c,:) = vMean.^2/2./vVar(c,:);
    end
    mGG = queryGz(vMean2DVar,model.gz,model.dgz);
    modelD.vVar{u} = vVar;
    modelD.mTmp_basef{u}= 0.5*exp(-model.const.C)*exp(-mGG).*vVar;
    modelD.mTmp_baseLogf{u} = -model.const.C-mGG+log(vVar)+log(0.5);
end
end