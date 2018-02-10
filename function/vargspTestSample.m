function [f,arrFUser] = vargspTestSample(data,model,options)
%VARGSPTEST BY SAMPLING
debug = 0;
% store result
modelD.f = 0;

D = length(model.Tend);
C = 0.5772156649;
% basic calculation
modelD = varsTransform(model,modelD,options);
switch options.inferType
    case 2
        modelD.sampleTime = 1;
    case 6
        modelD.sampleTime = 500;
end

modelD.Kmm = cell(1,model.K);
modelD.mPhi = cell(1,model.K);
modelD.diagKnn = cell(1,model.K);
modelD.vfInt = zeros(model.K,1);
for c=1:model.K
    nGamma = exp(2*model.GP.logtheta{c}(D+1));
    mTmp=model.var.Sigma{c}+model.var.Mu{c}*model.var.Mu{c}';
    modelD.Kmm{c} = feval(model.cov{:}, model.GP.logtheta{c},model.Xm);
    modelD.diagKnn{c} = sum(exp(2*model.GP.logtheta{c}(end-1:end)));
    modelD.mPhi{c} = MatrixPhi(model.GP.logtheta{c}, model);
    invKmmPhi = modelD.Kmm{c}\modelD.mPhi{c};
    invKmmPhiInvKmm = invKmmPhi/modelD.Kmm{c};
    modelD.vfInt(c) = nGamma*prod(model.Tend-model.Tstart)-trace(invKmmPhi)+...
        trace(invKmmPhiInvKmm*mTmp);
end
modelD.Knm = cell(model.U,model.K);
modelD.KnmInvKmm = cell(model.U,model.K);
modelD.L = cell(model.U,model.K);
modelD.Mu = cell(model.U,model.K);
modelD.vfpoint = cell(1,model.U);
for u=1:model.U
    nData = length(data{u}.X);
    vMean2DVar = zeros(model.K,nData);
    vVar = zeros(model.K,nData);
    for c=1:model.K
        modelD.Knm{u,c} = feval(model.cov{2}{1},model.GP.logtheta{c}(1:end-1),data{u}.X,model.Xm);
        modelD.KnmInvKmm{u,c} = modelD.Knm{u,c}/modelD.Kmm{c};
        vMean = (modelD.KnmInvKmm{u,c}*model.var.Mu{c})';
        mtmp = modelD.KnmInvKmm{u,c}*model.var.L{c};
        vVar(c,:) = modelD.diagKnn{c}-sum(modelD.KnmInvKmm{u,c}.*modelD.Knm{u,c},2)+sum(mtmp.^2,2);
        vMean2DVar(c,:) = vMean.^2/2./vVar(c,:);
    end
    mGG = queryGz(vMean2DVar,model.gz,model.dgz);
    modelD.vfpoint{u} = 0.5*exp(-C-mGG).*vVar;
end

modelD.sampleMu = cell(1,model.K);
arrFUser = zeros(1,model.U);
for u =1:model.U
    %sDetail.integral{u} = zeros(model.K,1);
    %sDetail.point{u} = zeros(length(data{u}.X),1);
    %sDetail.each{u} = 0;
end
f = 0;
for s = 1:modelD.sampleTime
    for u=1:model.U
        switch options.inferType
            case 2
                vAmp = modelD.mAmp_2(:,u); 
            case 6
                mtau = betarnd(model.amp.a(:,u),model.amp.b(:,u));
                m_tau = 1-mtau;
                m_tau = [1;cumprod(m_tau,1)];
                mtau = [mtau;1];
                vAmp = model.amp.s(u)*mtau.*m_tau;
        end
        f = f+sum(log(modelD.vfpoint{u}'*vAmp))-sum(vAmp.*modelD.vfInt);
        %sDetail.integral{u} = sDetail.integral{u}+vAmp.*modelD.vfInt;
        %sDetail.point{u} = sDetail.point{u}+modelD.vfpoint{u}'*vAmp;
        arrFUser(u) = arrFUser(u)+sum(log(modelD.vfpoint{u}'*vAmp))-sum(vAmp.*modelD.vfInt);
    end
end
f = f/modelD.sampleTime;
for u =1:model.U
    arrFUser(u) = arrFUser(u)/modelD.sampleTime;
    %sDetail.point{u} = sDetail.point{u}/modelD.sampleTime;
    %sDetail.each{u} = sDetail.each{u}/modelD.sampleTime;
end
end



