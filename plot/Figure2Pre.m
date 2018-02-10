clear;close all;clc;
addpath(genpath('.././'));
addpath(genpath('./'));
iplot = 14;
ifontsize = 16;
ilinewidth = 2;
iInfer = 6;

cPlot = cell(2,2);
for idataSource = 1:2
    dataSource = idataSource+1;
    switch dataSource
        case 2
            iSize = 500;
        case 3
            iSize = 600;
    end
    ifound = 1;
    str_folder = '.\ServerResult\';
    load(['dataSource',num2str(dataSource),'.mat']);
    load(['Source',num2str(dataSource),'Size',num2str(iSize),'Infer',num2str(iInfer),'.mat']);
    load(['Source',num2str(dataSource),'Size',num2str(iSize),'Infer',num2str(iInfer),'Cluster',num2str(iplot),'Model.mat']);
    model = model_save{ibest};
    modelD = calPre(model);
    modelD = calIntegral(model,options,modelD,3);
    modelD = varsTransform(model,modelD,options);
    
    vT = model.plot.Xm;
    fnRate=zeros(model.K,size(model.plot.Xm,1));
    for c=1:model.K
        Knm=feval(model.cov{2}{1},model.GP.logtheta{c}(1:end-1),model.plot.Xm,model.Xm);
        fnRate(c,:)=Knm*(modelD.Kmm{c}\model.var.Mu{c});
        fnRate(c,:) = fnRate(c,:).^2;
    end

    mMap = modelD.mMap;
    mMap = mMap.*repmat(modelD.vfInt,1,model.U);
    mMap = mMap./repmat(sum(mMap,1),model.K,1);
    
    % find two cases
    switch dataSource
        case 2
            arrIn = [1,6];
        case 3
            arrIn = [1,2];
    end
    ifoundstart = ifound;
    for i=1:length(model.NumData)
        bplot = 0;
        if mMap(arrIn(1),i) > 0.3 && mMap(arrIn(2),i) < 0.3 && ifound == ifoundstart && length(model.data{i}.feature) > 70
            bplot = 1;
        elseif mMap(arrIn(1),i) > 0.1 &&mMap(arrIn(2),i) > 0.2 && ifound == ifoundstart+1 && length(model.data{i}.feature) > 80
            bplot = 1;
        end
        if bplot == 1
            vAmp = mMap(:,i);
            vAmp = model.amp.s(i)*vAmp;
            fnRate_u = fnRate'*vAmp;
            
            cPlot{idataSource,ifound}.vT = vT;
            cPlot{idataSource,ifound}.fnRate = fnRate_u;
            cPlot{idataSource,ifound}.X = data.train{i}.X;
            if dataSource==2
                cPlot{idataSource,ifound}.hour = data.hour(i);
            end
            ifound = ifound+1;
        end
        if ifound == ifoundstart+2
            break
        end
    end
end
save('.\MiddleMat\Figure2.mat','cPlot');
