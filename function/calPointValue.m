function modelD = calPointValue(model,modelD,u,bflag,options)
% CALCULATE POINT VALUE AT EAmodel.const.CH POISTION
% u--user
% bflag -- 0 update f
%       -- 1 dmu,dL
%       -- 2 dkernel
%       -- 3 damp
% Time multiplied By number of elements
debug = 0;
nData = size(model.data{u}.feature,1);
vMean = zeros(model.K,nData);
vMean2DVar = zeros(model.K,nData);
vVar = zeros(model.K,nData);
modelD.Knm = cell(1,model.K);
modelD.KnmInvKmm=cell(1,model.K);

if bflag == 3
    vVar = modelD.vVar{u};
    mTmp_basef = modelD.mTmp_basef{u};
else
    for c=1:model.K
        if isfield(modelD,'userKnm')
            modelD.Knm{c} = modelD.userKnm{u,c};
            modelD.KnmInvKmm{c} = modelD.userKnmInvKmm{u,c};
        else
            switch bflag
                case 0
                    modelD.Knm{c} = feval(model.cov{2}{1},model.GP.logtheta{c}(1:end-1),model.data{u}.feature,model.Xm);
                case 2 % calculate the derivatives together
                    sRes = calCovSEard(model.GP.logtheta{c}(1:end-1),model.data{u}.feature,model.Xm);
                    modelD.Knm{c} = sRes.Knm;
                    modelD.DknmDx{c} = sRes.DknmDx;
            end
            modelD.KnmInvKmm{c} = modelD.Knm{c}/modelD.Kmm{c};
        end
        vMean(c,:) = modelD.Knm{c}*modelD.InvKmmMu{c};
        vVar(c,:) = modelD.diagKnn{c}+sum(modelD.KnmInvKmm{c}.*(modelD.Knm{c}*(modelD.InvKmmSigma{c}-eye(model.m))),2);
        vMean2DVar(c,:) = vMean(c,:).^2/2./vVar(c,:);
    end
    [mGG,mDGG] = queryGz(vMean2DVar,model.gz,model.dgz);
    mTmp_basef = 0.5*exp(-model.const.C)*exp(-mGG).*vVar;
end
switch options.inferType
    case 2
        mTmp = repmat(modelD.mAmp_2(:,u),1,model.NumData(u)).*mTmp_basef;
        vTmpSum1 = sum(mTmp,1);
        modelD.f = modelD.f+sum(log(vTmpSum1));
        
        
        mTmp_sum = repmat(vTmpSum1,model.K,1);
        mTmp = mTmp./mTmp_sum;
        switch bflag
            case 1
                mMutmp = mTmp.*mDGG;
                mLtmp = 2*mTmp.*(-mDGG.*vMean2DVar+1)./vVar;
                for c=1:model.K
                    vLtmp=mLtmp(c,:)';
                    vMutmp=(mMutmp(c,:)./vVar(c,:))';
                    mMulti=repmat(vLtmp,1,model.m)';
                    modelD.var.Mu{c} = modelD.var.Mu{c}+modelD.KnmInvKmm{c}'*(vMutmp.*vMean(c,:)');
                    modelD.var.L{c}=modelD.var.L{c}+(mMulti.*modelD.KnmInvKmm{c}')*modelD.KnmInvKmm{c};
                end
            case 2
                mDMutmp = mTmp.*mDGG;
                mDLtmp = 2*mTmp.*(-mDGG.*vMean2DVar+1)./vVar;
                for c=1:model.K
                    DbDxtmp = (2*modelD.Knm{c}*modelD.InvKmmSigma{c}-modelD.Knm{c})/modelD.Kmm{c};
                    vDMutmp=mDMutmp(c,:)';
                    vDLtmp=mDLtmp(c,:)';
                    vDMutmp=vDMutmp./vVar(c,:)';
                    
                    for k=1:modelD.numPar
                        
                        if k>=modelD.numPar-1
                            DknnDx=2*exp(2*model.GP.logtheta{c}(k))*ones(model.NumData(u),1);
                        else % except last parameter noise
                            DknnDx=zeros(model.NumData(u),1);
                        end
                        %% HAHA
                        DknmDx = modelD.DknmDx{c}{k};
                        DknmkmmInvDx = (-modelD.KnmInvKmm{c}*modelD.kmmDx{c}{k} + DknmDx);
                        DbDx = DknnDx -sum(modelD.KnmInvKmm{c}.*DknmDx,2)+...
                            sum(DknmkmmInvDx.*DbDxtmp,2);
                        DaDx = DknmkmmInvDx*modelD.InvKmmMu{c};
                        modelD.GP.logtheta{c}(k)=modelD.GP.logtheta{c}(k)+...
                            (vMean(c,:)'.*DaDx)'*vDMutmp+DbDx'*vDLtmp/2;
                    end
                end
            case 3
                modelD.amp(:,u) = modelD.amp(:,u)+sum(mTmp_basef./mTmp_sum,2);
        end

    case 6 % Gamma Process Prior
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mTmp = repmat(modelD.mAmp(:,u),1,model.NumData(u)).*mTmp_basef;
        vTmpSum1 = sum(mTmp,1);
        modelD.f = modelD.f+sum(log(vTmpSum1));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mTmp_sum = repmat(vTmpSum1,model.K,1);
        mTmp = mTmp./mTmp_sum;
        switch bflag
            case 1
                mMutmp = mTmp.*mDGG;
                mLtmp = 2*mTmp.*(-mDGG.*vMean2DVar+1)./vVar;
                for c=1:model.K
                    vLtmp=mLtmp(c,:)';
                    vMutmp=(mMutmp(c,:)./vVar(c,:))';
                    mMulti=repmat(vLtmp,1,model.m)';
                    modelD.var.Mu{c} = modelD.var.Mu{c}+modelD.KnmInvKmm{c}'*(vMutmp.*vMean(c,:)');
                    modelD.var.L{c}=modelD.var.L{c}+(mMulti.*modelD.KnmInvKmm{c}')*modelD.KnmInvKmm{c};
                end
            case 2
                mDMutmp = mTmp.*mDGG;
                mDLtmp = 2*mTmp.*(-mDGG.*vMean2DVar+1)./vVar;
                for c=1:model.K
                    DbDxtmp = 2*modelD.Knm{c}*modelD.InvKmmSigmaInvKmm{c}-modelD.KnmInvKmm{c};
                    vDMutmp=mDMutmp(c,:)';
                    vDLtmp=mDLtmp(c,:)';
                    vDMutmp=vDMutmp./vVar(c,:)';
                    for k=1:modelD.numPar
                        if k>=modelD.numPar-1
                            DknnDx=2*exp(2*model.GP.logtheta{c}(k))*ones(model.NumData(u),1);
                        else % except last parameter noise
                            DknnDx=zeros(model.NumData(u),1);
                        end
                        DknmDx = modelD.DknmDx{c}{k};
                        DknmkmmInvDx = -modelD.KnmInvKmm{c}*modelD.kmmDx{c}{k} + DknmDx;
                        DbDx = DknnDx -sum(modelD.KnmInvKmm{c}.*DknmDx,2)+...
                            sum(DknmkmmInvDx.*DbDxtmp,2);
                        DaDx = DknmkmmInvDx*modelD.InvKmmMu{c};
                        modelD.GP.logtheta{c}(k)=modelD.GP.logtheta{c}(k)+...
                            (vMean(c,:)'.*DaDx)'*vDMutmp+DbDx'*vDLtmp/2;
                    end
                end
            case 3
                modelD.f = modelD.f+(model.K-1)*log(model.prior.alpha)+...
                    sum(gammaln(model.amp.a(:,u))+gammaln(model.amp.b(:,u))-gammaln(model.amp.a(:,u)+model.amp.b(:,u)))...
                    +sum((model.prior.alpha-model.amp.b(:,u)).*(modelD.psiTauB(:,u)-modelD.psiTauAB(:,u))...
                    -(model.amp.a(:,u)-1).*(modelD.psiTauA(:,u)-modelD.psiTauAB(:,u)));
                
                modelD.amp.a(:,u) = modelD.amp.a(:,u)...
                    +(1+sum(mTmp(1:model.K-1,:),2)-model.amp.a(:,u)).*psi(1,model.amp.a(:,u))...
                    +(model.amp.a(:,u)+model.amp.b(:,u)-1-model.prior.alpha-model.NumData(u)+...
                    [0;sum(cumsum(mTmp(1:model.K-2,:),1),2)]).*psi(1,model.amp.a(:,u)+model.amp.b(:,u));
                modelD.amp.b(:,u) = modelD.amp.b(:,u)...
                    +(model.prior.alpha-model.amp.b(:,u)+model.NumData(u)-sum(cumsum(mTmp(1:model.K-1,:),1),2)).*psi(1,model.amp.b(:,u))...
                    +(model.amp.a(:,u)+model.amp.b(:,u)-1-model.prior.alpha-model.NumData(u)+...
                    [0;sum(cumsum(mTmp(1:model.K-2,:),1),2)]).*psi(1,model.amp.a(:,u)+model.amp.b(:,u));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                vfIntMap = modelD.vfInt.*modelD.mMap(:,u);
                modelD.f = modelD.f-model.amp.s(u)*(model.prior.b+sum(vfIntMap));
                
                vAinvB = model.amp.b(:,u)./model.amp.a(:,u);
                vfIntMapCumSum = sum(vfIntMap)-cumsum(vfIntMap(1:model.K-1));
                modelD.amp.a(:,u) = modelD.amp.a(:,u)-model.amp.s(u)*...
                    (vfIntMap(1:model.K-1).*vAinvB.*modelD.mInvAmpAB(:,u)-...
                    vfIntMapCumSum.*modelD.mInvAmpAB(:,u));
                modelD.amp.b(:,u) = modelD.amp.b(:,u)-model.amp.s(u)*...
                    (-vfIntMap(1:model.K-1).*modelD.mInvAmpAB(:,u)...
                    +vfIntMapCumSum.*modelD.mInvAmpAB(:,u)./vAinvB);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
end
end