function [LB,UB] = varsGpPpQueryRange(model,options,fig)
modelD.Kmm = cell(1,model.K);
modelD.diagKnn = cell(1,model.K);
for c=1:model.K
    modelD.Kmm{c} = feval(model.cov{:}, model.GP.logtheta{c},model.Xm);
    modelD.diagKnn{c} = sum(exp(2*model.GP.logtheta{c}(end-1:end)));
end

modelD.KnmInvKmm=cell(1,model.K);
LB_save = zeros(model.U,1);
UB_save = zeros(model.U,1);
for u=1:model.U
    nData = size(model.data{u}.feature,1);
    vMean2DVar = zeros(model.K,nData);
    for c=1:model.K
        
        modelD.Knm{c} = feval(model.cov{2}{1},model.GP.logtheta{c}(1:end-1),model.data{u}.feature,model.Xm);
        modelD.KnmInvKmm{c} = modelD.Knm{c}/modelD.Kmm{c};
        vMean = modelD.KnmInvKmm{c}*model.var.Mu{c};
        mtmp = modelD.KnmInvKmm{c}*model.var.L{c};
        vVar = modelD.diagKnn{c}-sum(modelD.KnmInvKmm{c}.*modelD.Knm{c},2)+sum(mtmp.^2,2);
        vMean2DVar(c,:) = vMean.^2/2./vVar;
    end
    
    vlog = log10(vMean2DVar);
    LB_save(u) = min(min(vlog));
    UB_save(u) = max(max(vlog));
end
LB = min(LB_save);
UB = min(UB_save);
if nargin == 3
    figure(fig);
    subplot(1,2,1);plot(LB_save,'b-+');
    subplot(1,2,2);plot(UB_save,'b-+');
end
end




