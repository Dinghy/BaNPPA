function modelD = varsTransform(model,modelD,options)
% TRANSORM PARAMETERS INTO PROPER FORM
% initialization
switch options.inferType
    case 2
        modelD.mAmp_2 = model.amp;
        modelD.mAmp = model.amp;
        modelD.mMap = model.amp;
    case 6
        % expectation log
        modelD.psiTauA = psi(model.amp.a);
        modelD.psiTauB = psi(model.amp.b);
        modelD.psiTauAB = psi(model.amp.a+model.amp.b);
        modelD.mAmp = zeros(model.K,model.U);
        modelD.mAmp(1:model.K-1,:) = modelD.mAmp(1:model.K-1,:)+...
            modelD.psiTauA-modelD.psiTauAB;
        modelD.mAmp(2:model.K,:) = modelD.mAmp(2:model.K,:)+...
            cumsum(modelD.psiTauB,1)-cumsum(modelD.psiTauAB,1);
        modelD.mAmp = exp(modelD.mAmp);
        % expectation
        modelD.mInvAmpAB = 1./(model.amp.a+model.amp.b);
        mtau = model.amp.a.*modelD.mInvAmpAB;
        m_tau = 1-mtau;
        m_tau = [ones(1,model.U);cumprod(m_tau,1)];
        mtau = [mtau;ones(1,model.U)];
        modelD.mMap = mtau.*m_tau;
        % mAmp size: K*U
end