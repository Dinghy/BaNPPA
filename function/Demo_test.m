function Demo_test(model,options)
modelD = calPre(model);
modelD = varsTransform(model,modelD,options);
% one method
mtau = model.amp.a./(model.amp.a+model.amp.b);
m_tau = 1-mtau;
m_tau = [ones(1,model.U);cumprod(m_tau,1)];
mtau = [mtau;ones(1,model.U)];
vampA = mtau.*m_tau;
mTmpA = cell(1,model.U);
for u = 1:model.U
    mTmp_basef = modelD.mTmp_basef{u};
    mTmpA{u} = repmat(vampA(:,u),1,model.NumData(u)).*mTmp_basef;
    vTmpSum1 = sum(mTmpA{u},1);
    mTmpA{u} = mTmpA{u}./repmat(vTmpSum1,model.K,1);
end
% the other method
mtau = log(model.amp.a)-log(model.amp.a+model.amp.b);
m_tau = log(model.amp.b)-log(model.amp.a+model.amp.b);
m_tau = [zeros(1,model.U);cumsum(m_tau,1)];
mtau = [mtau;zeros(1,model.U)];
vampB = mtau+m_tau;
mTmpB = cell(1,model.U);
for u = 1:model.U
    mTmp_baseLogf = modelD.mTmp_baseLogf{u};
    mTmpB{u} = repmat(vampB(:,u),1,model.NumData(u))+mTmp_baseLogf;
    mTmpB{u} = mTmpB{u}-repmat(max(mTmpB{u},[],1),model.K,1);
    mTmpB{u} = exp(mTmpB{u});
    vTmpSum1 = sum(mTmpB{u},1);
    mTmpB{u} = mTmpB{u}./repmat(vTmpSum1,model.K,1);
end
%% result
disp(sum(sum(abs(vampA-exp(vampB)))));
disp(sum(cellfun(@(x,y)sum(sum(abs(x-y))),mTmpA,mTmpB)));
end