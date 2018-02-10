clear;close all;clc;
addpath(genpath('./'));
addpath(genpath('.././'));

idataset = 9;
switch idataset
    case 8
        iclu = 3;
        isave = 11;
        strSource = 'H';
    case 9
        iclu = 2;
        isave = 12;
        strSource = 'I';
end
    
load(['Source',num2str(idataset),'Size200Infer2Cluster',num2str(iclu),'model.mat']);
model = model_save{ibest};
% basic calculation
cPlot.fnRate=zeros(model.K,size(model.plot.Xm,1));
cPlot.fnLow = zeros(model.K,size(model.plot.Xm,1));
cPlot.fnHigh = zeros(model.K,size(model.plot.Xm,1));
for c= 1:model.K
    c
    modelD.diagKnn{c} = sum(exp(2*model.GP.logtheta{c}(end-1:end)));
    modelD.Kmm{c} = feval(model.cov{:}, model.GP.logtheta{c},model.Xm);
    modelD.InvKmmMu{c} = modelD.Kmm{c}\model.var.Mu{c};
    modelD.InvKmmSigma{c} = modelD.Kmm{c}\model.var.Sigma{c};
    modelD.InvKmmL{c} = modelD.Kmm{c}\model.var.L{c};

    Knm=feval(model.cov{2}{1},model.GP.logtheta{c}(1:end-1),model.plot.Xm,model.Xm);
    vMean=Knm*(modelD.Kmm{c}\model.var.Mu{c});
    KnmInvKmm = Knm/modelD.Kmm{c};
    vVar = modelD.diagKnn{c}...
        +sum(KnmInvKmm.*(Knm*(modelD.InvKmmSigma{c}-eye(model.m))),2);
    cPlot.fnRate(c,:) = vMean.^2+vVar;
    cPlot.fnLow(c,:) = ncx2inv(0.05,1,vMean.^2./vVar).*vVar;
    cPlot.fnHigh(c,:) = ncx2inv(0.95,1,vMean.^2./vVar).*vVar;
end
cPlot.K = model.K;
cPlot.vT = model.plot.Xm;
cPlot.logtheta = model.GP.logtheta;
% True
load(['dataSource1',strSource,'.mat']);
TrueK = length(data.logtheta);
cPlot.TrueK = TrueK;
cPlot.TrueMean = cell(1,TrueK);
cPlot.TrueLow = cell(1,TrueK);
cPlot.TrueHigh = cell(1,TrueK);
for c = 1:TrueK
    c
    [~,fMean,fHigh,fLow] = fGpInt(data.GP{c},cPlot.vT);
    cPlot.TrueMean{c} = fMean;
    cPlot.TrueLow{c} = fLow;
    cPlot.TrueHigh{c} = fHigh;
end
save(['.\MiddleMat\Figure',num2str(isave),'.mat'],'cPlot');
