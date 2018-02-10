function [vF,vDiv] = calLikelihood(model,options)
vDiv = zeros(model.K,1);
D = length(model.Tend);
vF = zeros(1,5);
 % fprintf('Point %.4f\tIntegral %.4f\tPenalty %.4f\tKLA %.4f\tKLMu %.4f\n',vF(1),vF(2),vF(3),vF(4),vF(5));
modelD.f = 0;
modelD = varsTransform(model,modelD,options);
modelD.InvKmmMu = cell(1,model.K);
modelD.InvKmmSigma = cell(1,model.K);
modelD.InvKmmL = cell(1,model.K);
modelD.diagKnn=cell(1,model.K);
modelD.Kmm = cell(1,model.K);
for c= 1:model.K
    modelD.diagKnn{c} = sum(exp(2*model.GP.logtheta{c}(end-1:end)));
    modelD.Kmm{c} = feval(model.cov{:}, model.GP.logtheta{c},model.Xm);
    modelD.InvKmmMu{c} = modelD.Kmm{c}\model.var.Mu{c};
    modelD.InvKmmSigma{c} = modelD.Kmm{c}\model.var.Sigma{c};
    modelD.InvKmmL{c} = modelD.Kmm{c}\model.var.L{c};
end
switch options.inferType
    case 2
        for u=1:model.U
            nData = size(model.data{u}.feature,1);
            vMean = zeros(model.K,nData);
            vMean2DVar = zeros(model.K,nData);
            vVar = zeros(model.K,nData);
            for c=1:model.K
                modelD.Knm{c} = feval(model.cov{2}{1},model.GP.logtheta{c}(1:end-1),model.data{u}.feature,model.Xm);
                modelD.KnmInvKmm{c} = modelD.Knm{c}/modelD.Kmm{c};
                vMean(c,:) = modelD.Knm{c}*modelD.InvKmmMu{c};
                vVar(c,:) = modelD.diagKnn{c}+sum(modelD.KnmInvKmm{c}.*(modelD.Knm{c}*(modelD.InvKmmSigma{c}-eye(model.m))),2);
                vMean2DVar(c,:) = vMean(c,:).^2/2./vVar(c,:);
            end
            mGG = queryGz(vMean2DVar,model.gz,model.dgz);
            mTmp_basef = 0.5*exp(-model.const.C)*exp(-mGG).*vVar;
            mTmp = repmat(modelD.mAmp_2(:,u),1,model.NumData(u)).*mTmp_basef;
            vTmpSum1 = sum(mTmp,1);
            vF(1) = vF(1) +sum(log(vTmpSum1));
        end
        
        vfInt = zeros(model.K,1);
        for c=1:model.K
            nGamma = exp(2*model.GP.logtheta{c}(D+1));
            InvKmmSigmaMu=modelD.InvKmmSigma{c}+modelD.InvKmmMu{c}*model.var.Mu{c}';
            % calculate basics
            mPhi = MatrixPhi(model.GP.logtheta{c}, model);
            invKmmPhi = modelD.Kmm{c}\mPhi;
            % record integral
            vfInt(c)=nGamma*prod(model.Tend-model.Tstart)-trace(invKmmPhi)+...
                trace(invKmmPhi*InvKmmSigmaMu);
        end
        for c=1:model.K
            nAmp=sum(modelD.mAmp_2(c,:));
            vF(2)=vF(2)-nAmp*vfInt(c);
        end
        % KL
        for c=1:model.K
            InvKmmMuG = modelD.Kmm{c}\(model.var.Mu{c}-model.prior.g);
            InvKmmSigmaMuG=modelD.InvKmmSigma{c}+InvKmmMuG*(model.var.Mu{c}-model.prior.g)';
            vF(5) = vF(5)+0.5*(model.m-trace(InvKmmSigmaMuG))+sum(log(diag(model.var.L{c})))-0.5*logdet(modelD.Kmm{c});
        end
    case 6
        for u=1:model.U
            nData = size(model.data{u}.feature,1);
            vMean = zeros(model.K,nData);
            vMean2DVar = zeros(model.K,nData);
            vVar = zeros(model.K,nData);
            for c=1:model.K
                modelD.Knm{c} = feval(model.cov{2}{1},model.GP.logtheta{c}(1:end-1),model.data{u}.feature,model.Xm);
                modelD.KnmInvKmm{c} = modelD.Knm{c}/modelD.Kmm{c};
                vMean(c,:) = modelD.Knm{c}*modelD.InvKmmMu{c};
                vVar(c,:) = modelD.diagKnn{c}+sum(modelD.KnmInvKmm{c}.*(modelD.Knm{c}*(modelD.InvKmmSigma{c}-eye(model.m))),2);
                vMean2DVar(c,:) = vMean(c,:).^2/2./vVar(c,:);
            end
            mGG = queryGz(vMean2DVar,model.gz,model.dgz);
            mTmp_basef = 0.5*exp(-model.const.C)*exp(-mGG).*vVar;
            mTmp = repmat(modelD.mAmp(:,u),1,nData).*mTmp_basef;
            vTmpSum1 = sum(mTmp,1);
            vF(1) = vF(1) +sum(log(vTmpSum1))+nData*log(model.amp.s(u));
        end
        
        vfInt = zeros(model.K,1);
        for c=1:model.K
            nGamma = exp(2*model.GP.logtheta{c}(D+1));
            InvKmmSigmaMu=modelD.InvKmmSigma{c}+modelD.InvKmmMu{c}*model.var.Mu{c}';
            % calculate basics
            mPhi = MatrixPhi(model.GP.logtheta{c}, model);
            invKmmPhi = modelD.Kmm{c}\mPhi;
            % record integral
            vfInt(c)=nGamma*prod(model.Tend-model.Tstart)-trace(invKmmPhi)+...
                trace(invKmmPhi*InvKmmSigmaMu);
        end
        for u=1:model.U
            vF(2) = vF(2)+model.prior.a*log(model.prior.b)-gammaln(model.prior.a)-model.amp.s(u)*(sum(vfInt.*modelD.mMap(:,u))+model.prior.b)...
                +(model.prior.a-1)*log(model.amp.s(u));
        end
        
        % adding the lagrangian term
        if options.constraint == 1
            
            vF(3) = vF(3)-sum(model.la.lambda.*(vfInt-model.A))-0.5*model.la.c*sum((vfInt-model.A).^2);
            vDiv = vfInt-model.A;
        end
        % KL Theta
        vF(4) = vF(4) + model.U*(model.K-1)*log(model.prior.alpha)+...
            sum(sum(gammaln(model.amp.a)+gammaln(model.amp.b)-gammaln(model.amp.a+model.amp.b)))...
            +sum(sum((model.prior.alpha-model.amp.b).*(modelD.psiTauB-modelD.psiTauAB)...
            -(model.amp.a-1).*(modelD.psiTauA-modelD.psiTauAB)));
        % KL f
        for c=1:model.K
            InvKmmMuG = modelD.Kmm{c}\(model.var.Mu{c}-model.prior.g);
            InvKmmSigmaMuG=modelD.InvKmmSigma{c}+InvKmmMuG*(model.var.Mu{c}-model.prior.g)';
            vF(5) = vF(5)+0.5*(model.m-trace(InvKmmSigmaMuG))+sum(log(diag(model.var.L{c})))-0.5*logdet(modelD.Kmm{c});
        end
end

end