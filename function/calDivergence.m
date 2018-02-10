function modelD = calDivergence(model,modelD,bflag)
%  CALCULATE DIVIGENCE AND UPDATE DERIVATIVES
%  bflag = 0 NOTHING
%       = 1 dMU,dL 
%       = 2 dKernel
for c=1:model.K
    InvKmmMuG = modelD.Kmm{c}\(model.var.Mu{c}-model.prior.g);
    InvKmmSigmaMuG=modelD.InvKmmSigma{c}+InvKmmMuG*(model.var.Mu{c}-model.prior.g)';
    
    modelD.f=modelD.f+0.5*(model.m-trace(InvKmmSigmaMuG));
    switch bflag
        case 1
            modelD.f = modelD.f+sum(log(diag(model.var.L{c})));
            modelD.var.Mu{c} = modelD.var.Mu{c}-InvKmmMuG;
            modelD.var.L{c} = modelD.var.L{c}+(diag(1./diag(model.var.L{c}))-modelD.InvKmmL{c});
        case 2
            modelD.f = modelD.f-0.5*logdet(modelD.Kmm{c});
            mtmp = InvKmmSigmaMuG/modelD.Kmm{c};
            modelD.GP.logtheta{c} = modelD.GP.logtheta{c}-0.5*cellfun(@(x)trace(modelD.Kmm{c}\x-mtmp*x),modelD.kmmDx{c});
    end
end
end