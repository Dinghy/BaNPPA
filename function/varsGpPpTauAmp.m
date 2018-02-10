function model = varsGpPpTauAmp(model,options,modelD)
%%  CLOSED FORM UPDATE OF TAU
if nargin == 2
    for u=1:model.U
        mTmp = 1e-4+rand(model.K,model.NumData(u));
        mTmp = mTmp./repmat(sum(mTmp,1),model.K,1);
        model.amp.a(:,u) = 1+sum(mTmp(1:model.K-1,:),2);
        model.amp.b(:,u) = model.prior.alpha+model.NumData(u)-sum(cumsum(mTmp(1:model.K-1,:),1),2);
    end
else
    modelD = varsTransform(model,modelD,options);
    for u=1:model.U
        mTmp_basef = modelD.mTmp_basef{u};
        mTmp = repmat(modelD.mAmp(:,u),1,model.NumData(u)).*mTmp_basef;
        vTmpSum1 = sum(mTmp,1);
        mTmp = mTmp./repmat(vTmpSum1,model.K,1);
        model.amp.a(:,u) = 1+sum(mTmp(1:model.K-1,:),2);
        model.amp.b(:,u) = model.prior.alpha+model.NumData(u)-sum(cumsum(mTmp(1:model.K-1,:),1),2);
    end
end
end