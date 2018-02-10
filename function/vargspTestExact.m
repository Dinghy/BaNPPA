function [f,arrFUser] = vargspTestExact(data,model,options,fig,arruser)
% VARGSPTEST BY SAMPLING
debug = 0;
% store result
modelD.f = 0;

D = length(model.Tend);
C = 0.5772156649;
% basic calculation
modelD = varsTransform(model,modelD,options);
modelD.sampleTime = 200;


modelD.Kmm = cell(1,model.K);
modelD.diagKnn = cell(1,model.K);
modelD.vfInt = zeros(model.K,1);
for c=1:model.K
    nGamma = exp(2*model.GP.logtheta{c}(D+1));
    mTmp=model.var.Sigma{c}+model.var.Mu{c}*model.var.Mu{c}';
    modelD.Kmm{c} = feval(model.cov{:}, model.GP.logtheta{c},model.Xm);
    model.var.L{c} = chol(model.var.Sigma{c},'lower');
    mPhi = MatrixPhi(model.GP.logtheta{c}, model);
    invKmmPhi = modelD.Kmm{c}\mPhi;
    invKmmPhiInvKmm = invKmmPhi/modelD.Kmm{c};
    modelD.vfInt(c) = nGamma*prod(model.Tend-model.Tstart)-trace(invKmmPhi)+...
        trace(invKmmPhiInvKmm*mTmp);
end
arrFUser = zeros(1,model.U);

modelD.Knm = cell(model.U,model.K);
modelD.KnmInvKmm = cell(model.U,model.K);
modelD.L = cell(model.U,model.K);
modelD.vfpoint = cell(1,model.U);
for u=1:model.U
    for c=1:model.K
        modelD.Knm{u,c} = feval(model.cov{:}, model.GP.logtheta{c},data{u}.X,model.Xm);
        modelD.KnmInvKmm{u,c} = modelD.Knm{u,c}/modelD.Kmm{c};
        modelD.sigma{u,c} = sum(exp(2*model.GP.logtheta{c}(end-1:end)))-sum(modelD.KnmInvKmm{u,c}.*modelD.Knm{u,c},2);
        modelD.sigma{u,c} = sqrt(modelD.sigma{u,c});
    end
end

if nargin == 5
    fnRate=zeros(model.K,size(model.plot.Xm,1));
    for c=1:model.K
        Knm=feval(model.cov{2}{1},model.GP.logtheta{c}(1:end-1),model.plot.Xm,model.Xm);
        fnRate(c,:)=Knm*(modelD.Kmm{c}\model.var.Mu{c});
        fnRate(c,:) = fnRate(c,:).^2;
    end
    figure(fig);hold on;
    f = 0;
    vM = zeros(model.m,1);
    mM = eye(model.m);
    sampleM = cell(1,model.K);
    for c=1:model.K
        sampleM{c} = zeros(size(model.var.Mu{c}));
    end
    for s = 1:modelD.sampleTime
        for c=1:model.K
            sampleM{c} = model.var.Mu{c}+model.var.L{c}*mvnrnd(vM,mM)';
        end
        for iU = 1:length(arruser)
            u = arruser(iU);
            nData = length(data{u}.X);
            vSampleN = zeros(nData,1);
            vAmp = modelD.mAmp(:,u);
            for c=1:model.K
                vSampleN = vSampleN + vAmp(c)*(modelD.KnmInvKmm{u,c}*sampleM{c}+modelD.sigma{u,c}.*randn(nData,1)).^2;
            end
            f = f+sum(log(vSampleN));
            % plot for test
            fnRate_u = fnRate'*vAmp;
            plot(model.plot.Xm,fnRate_u,'r');
            plot(data{u}.X,vSampleN,'b.');
        end
    end
else
    f = 0;
    vM = zeros(model.m,1);
    mM = eye(model.m);
    sampleM = cell(1,model.K);
    for c=1:model.K
        sampleM{c} = zeros(size(model.var.Mu{c}));
    end
    for s = 1:modelD.sampleTime
        for c=1:model.K
            sampleM{c} = model.var.Mu{c}+model.var.L{c}*mvnrnd(vM,mM)';
        end
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
            
            nData = length(data{u}.X);
            vSampleN = zeros(nData,1);
%             vAmp = modelD.mAmp(:,u);
            for c=1:model.K
                vSampleN = vSampleN + vAmp(c)*(modelD.KnmInvKmm{u,c}*sampleM{c}+modelD.sigma{u,c}.*randn(nData,1)).^2;
            end
            f = f+sum(log(vSampleN))-sum(vAmp.*modelD.vfInt);
            arrFUser(u) = arrFUser(u)+sum(log(vSampleN))-sum(vAmp.*modelD.vfInt);
        end
    end
end
f = f/modelD.sampleTime;
for u=1:model.U
    arrFUser(u) = arrFUser(u)/modelD.sampleTime;
end
end



