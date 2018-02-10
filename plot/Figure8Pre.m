%% plot for one experiment
clear;close all;clc;
addpath(genpath('./'));
addpath(genpath('.././'));
dataSource = 5;
iplot = 14;
iRow = 3;
strTitle = 'Synthetic B';
strMethod = {'BaNPPA','BaNPPA-'};
dataU = 250;
nylim = 0.5;
arr_infer = [6,7];

cPlot = cell(1,length(arr_infer));
for i=1:length(arr_infer)
    cPlot{i}.NER = zeros(iplot,5);
    cPlot{i}.UNER = zeros(iplot,5);
    cPlot{i}.fInt = zeros(iplot,5);
end

for i = 1:length(arr_infer)
    iInfer = arr_infer(i);
    load(['Alpha',num2str(dataSource),'Size',num2str(dataU),'Infer',num2str(iInfer),'.mat']);
    
    load(['Alpha',num2str(dataSource),'Source',num2str(dataSource),'Size',num2str(dataU),...
        'Infer',num2str(iInfer),'Cluster',num2str(iplot),'Model.mat']);
    for isave = 1:5
        [i,isave]
        
        model = model_save{isave};
        modelD = calPre(model);
        modelD = calIntegral(model,options,modelD,3);
        modelD = varsTransform(model,modelD,options);
        
        mMap = zeros(model.K,model.U);
        switch options.inferType
            case 2
                mMap = modelD.mAmp_2;
            case {5,6}
                mMap = modelD.mMap;
        end
        mMapO = mMap;
        mMapO = mMapO./repmat(sum(mMapO,1),model.K,1);
        vOccO = sum(mMapO,2)/model.U;
        
        
        mMap = mMap.*repmat(modelD.vfInt,1,model.U);
        mMap = mMap./repmat(sum(mMap,1),model.K,1);
        vOcc = sum(mMap,2)/model.U;
        
        cPlot{i}.UNER(:,isave) = vOccO;
        cPlot{i}.NER(:,isave) = vOcc;
        cPlot{i}.fInt(:,isave) = modelD.vfInt;
        cPlot{i}.best = ibest;
        % record the first four basis functions for the best result
        if isave == ibest
            [~,I] = sort(vOcc,'descend');
            vT = model.plot.Xm;
            fnRate=zeros(model.K,size(model.plot.Xm,1));
            for ic=1:model.K
                c=I(ic);
                Knm=feval(model.cov{2}{1},model.GP.logtheta{c}(1:end-1),model.plot.Xm,model.Xm);
                fnRate(ic,:)=Knm*(modelD.Kmm{c}\model.var.Mu{c});
                fnRate(ic,:) = fnRate(ic,:).^2;
            end
            cPlot{i}.fnRate=fnRate;
            cPlot{i}.vT = vT;
        end
    end
end
save('.\MiddleMat\Figure7.mat','strTitle','strMethod','cPlot');     