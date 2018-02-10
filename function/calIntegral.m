function modelD = calIntegral(model,options,modelD,bflag)
% CALCULATE INTEGRATION PART AND UPDATE DERIVATIVES
% bflag = 0 NOTHING (Test
%       = 1 dMU,dL  (VarPar
%       = 2 dKernel (Hypers
%       = 3 Outer
switch options.inferType
    case 2
        D = length(model.Tend);
        if bflag == 3 % outer
            % initialization
            modelD.f=0;
            modelD.vfInt = zeros(model.K,1);
            modelD.invKmmPhi = cell(1,model.K);
        end
        
        for c=1:model.K
            nGamma = exp(2*model.GP.logtheta{c}(D+1));
            switch bflag
                case {0,3}
                    mPhi = MatrixPhi(model.GP.logtheta{c}, model);
                    invKmmPhi = modelD.Kmm{c}\mPhi;
                    modelD.invKmmPhi{c} = invKmmPhi;
                case 1
                    invKmmPhi = modelD.invKmmPhi{c};
                case 2
                    [mPhi,DmPhidx] = MatrixPhi(model.GP.logtheta{c}, model);
                    invKmmPhi = modelD.Kmm{c}\mPhi;
            end
            
            InvKmmSigmaMu=modelD.InvKmmSigma{c}+modelD.InvKmmMu{c}*model.var.Mu{c}';
            
            nfInt=nGamma*prod(model.Tend-model.Tstart)-trace(invKmmPhi)+...
                trace(invKmmPhi*InvKmmSigmaMu);
            
            switch bflag
                case {0,3}
                    modelD.vfInt(c) = nfInt;
                case 1
                    nAmp=sum(modelD.mAmp_2(c,:));
                    modelD.f=modelD.f-nAmp*nfInt;
                    modelD.var.Mu{c}= modelD.var.Mu{c}-2*nAmp*invKmmPhi*modelD.InvKmmMu{c};
                    modelD.var.L{c} = modelD.var.L{c}-2*nAmp*invKmmPhi*modelD.InvKmmL{c};
                case 2
                    nAmp=sum(modelD.mAmp_2(c,:));
                    modelD.f=modelD.f-nAmp*nfInt;
                    mtmp_a = (eye(model.m)-InvKmmSigmaMu)/modelD.Kmm{c};
                    mtmp_b = (2*invKmmPhi*InvKmmSigmaMu-invKmmPhi)/modelD.Kmm{c};
                    for k=[1:modelD.numPar-2,modelD.numPar]
                        modelD.GP.logtheta{c}(k)=modelD.GP.logtheta{c}(k)+nAmp*trace(mtmp_a*DmPhidx{k})...
                            +nAmp*trace(mtmp_b*modelD.kmmDx{c}{k});
                    end
            end
        end
    case 6 % Gamma Process Prior with or without regularization
        D = length(model.Tend);
        
        if bflag == 3 % outer
            % initialization
            modelD.f=0;
            modelD.vfInt = zeros(model.K,1);
            modelD.invKmmPhi = cell(1,model.K);
        else
            switch bflag
                case 1
                    dMutmp = cell(1,model.K);
                    dLtmp = cell(1,model.K);
                case 2
                    dftmp = zeros(model.K,modelD.numPar);
            end
        end
        vfInt = zeros(model.K,1);
        for c=1:model.K
            nGamma = exp(2*model.GP.logtheta{c}(D+1));
            InvKmmSigmaMu=modelD.InvKmmSigma{c}+modelD.InvKmmMu{c}*model.var.Mu{c}';
            % calculate basics
            switch bflag
                case {0,3}
                    mPhi = MatrixPhi(model.GP.logtheta{c}, model);
                    invKmmPhi = modelD.Kmm{c}\mPhi;
                    modelD.invKmmPhi{c} = invKmmPhi;
                case 1
                    invKmmPhi = modelD.invKmmPhi{c};
                    dMutmp{c} = 2*invKmmPhi*modelD.InvKmmMu{c};
                    dLtmp{c} = 2*invKmmPhi*modelD.InvKmmL{c};
                case 2
                    [mPhi,DmPhidx] = MatrixPhi(model.GP.logtheta{c}, model);
                    invKmmPhi = modelD.Kmm{c}\mPhi;
                    mtmp_a = (eye(model.m)-InvKmmSigmaMu)/modelD.Kmm{c};
                    mtmp_b = (2*invKmmPhi*InvKmmSigmaMu-invKmmPhi)/modelD.Kmm{c};
                    for k=[1:modelD.numPar-2,modelD.numPar]
                        dftmp(c,k) = trace(mtmp_a*DmPhidx{k})...
                            +trace(mtmp_b*modelD.kmmDx{c}{k});
                    end
            end
            % record integral
            vfInt(c)=nGamma*prod(model.Tend-model.Tstart)-trace(invKmmPhi)+...
                trace(invKmmPhi*InvKmmSigmaMu);        
        end
        switch bflag
            case {0,3}
                modelD.vfInt = vfInt;
            case 1
                vRatio = zeros(model.K,1);
                for u=1:model.U
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    LD = model.prior.b+sum(vfInt.*modelD.mMap(:,u));
                    modelD.f = modelD.f-model.amp.s(u)*LD;
                    vRatio = vRatio+model.amp.s(u)*modelD.mMap(:,u);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
                for c=1:model.K
                    modelD.var.Mu{c} = modelD.var.Mu{c}-vRatio(c)*dMutmp{c};
                    modelD.var.L{c} = modelD.var.L{c}-vRatio(c)*dLtmp{c};
                end
                if options.constraint == 1
                    for c=1:model.K
                        % adding the lagrangian term
                        modelD.f = modelD.f-model.la.lambda(c)*(vfInt(c)-model.A)-0.5*model.la.c*(vfInt(c)-model.A)^2;
                        nRatio = model.la.lambda(c)+model.la.c*(vfInt(c)-model.A);
                        modelD.var.Mu{c} = modelD.var.Mu{c}-nRatio*dMutmp{c};
                        modelD.var.L{c} = modelD.var.L{c}-nRatio*dLtmp{c};
                    end
                end
            case 2
                vRatio = zeros(model.K,1);
                for u=1:model.U
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    LD = model.prior.b+sum(vfInt.*modelD.mMap(:,u));
                    modelD.f = modelD.f-model.amp.s(u)*LD;
                    vRatio = vRatio+model.amp.s(u)*modelD.mMap(:,u);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
                for c=1:model.K
                    modelD.GP.logtheta{c} = modelD.GP.logtheta{c} + vRatio(c)*dftmp(c,:)';
                end
                if options.constraint==1
                    for c=1:model.K
                        % adding the lagrangian term
                        modelD.f = modelD.f-model.la.lambda(c)*(vfInt(c)-model.A)-0.5*model.la.c*(vfInt(c)-model.A)^2;
                        nRatio = model.la.lambda(c)+model.la.c*(vfInt(c)-model.A);
                        modelD.GP.logtheta{c} = modelD.GP.logtheta{c}+nRatio*dftmp(c,:)';
                    end
                end
        end
end
end