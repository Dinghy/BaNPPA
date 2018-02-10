clear;close all;clc;
addpath(genpath('./'));
addpath(genpath('.././'));

cPlot.Hyper = cell(4,3);
cPlot.Map = cell(4,3);
cPlot.nStartHyper = zeros(1,4);
arrInfer = [2,7,6];
strTitle = {'Synthetic A','Synthetic B','Microblog','Citation'};
strMethod = {'LPPA','BaNPPA-NC','BaNPPA'};
for idataset = 1:4
    switch idataset
        case 1
            index = 1;
            iU = 200;
            dd = log((60)')/2;
        case 2
            index = 5;
            iU = 250;
            dd = log(4.3081);
        case 3
            index = 2;
            iU = 500;
            dd = log((12)')/2;
        case 4
            index = 3;
            iU = 600;
            dd = log((20)')/2;
    end
    cPlot.nStartHyper(idataset) = exp(dd);
    for i=1:3
        [idataset,i]
        infer = arrInfer(i);
        load(['Source',num2str(index),'Size',num2str(iU),'Infer',num2str(infer),'Cluster14Model.mat']);
        model = model_save{ibest};
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
        mMap = mMap.*repmat(modelD.vfInt,1,model.U);
        mMap = mMap./repmat(sum(mMap,1),model.K,1);
        vOcc = sum(mMap,2)/model.U;
        
        cPlot.NER{idataset,i} = vOcc;
        
        arrHyper = cellfun(@(x)exp(x(1)),model.GP.logtheta);
        cPlot.Hyper{idataset,i} = arrHyper;
    end
end
cPlot.strTitle = strTitle;
cPlot.strMethod = strMethod;
save('.\MiddleMat\Figure9.mat','cPlot');