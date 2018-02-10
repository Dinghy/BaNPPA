%% plot for one experiment
clear;close all;clc;
addpath(genpath('./'));
addpath(genpath('.././'));

strTitle = {'Synthetic A','Synthetic B','Microblog','Citation'};
strMethod = {'LPPA','BaNPPA-','BaNPPA'};
arrDataset = [1,5,2,3];
iplotBasis = 4;
iplot = 14;

cPlot = cell(4,3);
for j=1:4
    for i=1:3
        cPlot{j,i}.NER = zeros(iplot,5);
        cPlot{j,i}.UNER = zeros(iplot,5);
        cPlot{j,i}.Map = cell(1,5);
        cPlot{j,i}.FInt = zeros(iplot,5);
        cPlot{j,i}.fBasis = cell(1,5);
    end
end

for idataSource = 1:4
    dataSource = arrDataset(idataSource);
    
    switch dataSource
        case 1
            strSource = 'Synthetic A';
            dataU = 200;
        case 2
            strSource = 'Microblog';
            dataU = 500;
        case 3
            strSource = 'Citation';
            dataU = 600;
        case 5
            strSource = 'Synthetic B';
            dataU = 250;
    end
    nylim = 0.5;
    arr_infer = [2,7,6];
    ifontsize = 13;
    ilinewidth = 2;
    
    strtitle = {'LPPA','BNPPA'};
    if iplotBasis == 4
        strstyle = {':','-.','-','--'};
    end
    
    for i = 1:length(arr_infer)
        iInfer = arr_infer(i);
        load(['Source',num2str(dataSource),'Size',num2str(dataU),'Infer',num2str(iInfer),'.mat']);
        load(['Source',num2str(dataSource),'Size',num2str(dataU),...
                'Infer',num2str(iInfer),'Cluster',num2str(iplot),'Model.mat']);
        for isave = 1:5
            [dataSource,i,isave]
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
            
            cPlot{idataSource,i}.Map{isave} = mMap;
            cPlot{idataSource,i}.UNER(:,isave) = vOccO;
            cPlot{idataSource,i}.NER(:,isave) = vOcc;
            cPlot{idataSource,i}.FInt(:,isave) = modelD.vfInt;
            cPlot{idataSource,i}.Best = ibest;
            
            %% record the Basis
            vT = model.plot.Xm;
            fnRate=zeros(model.K,size(model.plot.Xm,1));
            for c=1:model.K
                Knm=feval(model.cov{2}{1},model.GP.logtheta{c}(1:end-1),model.plot.Xm,model.Xm);
                fnRate(c,:)=Knm*(modelD.Kmm{c}\model.var.Mu{c});
                fnRate(c,:) = fnRate(c,:).^2;
            end
            cPlot{idataSource,i}.fBasisT = vT;
            cPlot{idataSource,i}.fBasis{isave}=fnRate;
        end
    end    
end
%%
strTitle = {'Synthetic A','Synthetic B','Microblog','Citation'};
strMethod = {'LPPA','BaNPPA-','BaNPPA'};
save('.\MiddleMat\Figure456.mat','cPlot','strTitle','strMethod');