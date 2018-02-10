%% plot execution time v.s. Train Likelihood/alpha
clear;close all;clc;
addpath(genpath('./'));
addpath(genpath('.././'));
arrSet = [1,5,2,3];
strTitle = {'Synthetic A','Synthetic B','Microblog','Citation'};
strMethod = {'LPPA','BaNPPA-','BaNPPA'};

cPlot = cell(4,3);

for idataSet = 1:4
    dataSet = arrSet(idataSet);
    switch dataSet
        case 1
            U = 200;
        case 2
            U = 500;
        case 3
            U = 600;
        case 5
            U = 250;
    end
    arrInfer = [2,7,6];
    %% plot the train likelihood
    
    for i = 1:length(arrInfer)
        infer = arrInfer(i);
        load(['Rec',num2str(dataSet),'Size',num2str(U),...
            'Infer',num2str(infer),'Cluster14Model.mat']);
        arrplot = mSaveTime > 0;
        vX = mSaveTime(arrplot);
        vY = mSaveTrain(arrplot);
        vAlpha = mSaveAlpha(arrplot);
        cPlot{idataSet,i}.vX = vX;
        cPlot{idataSet,i}.vY = vY;
        cPlot{idataSet,i}.vAlpha = vAlpha;
    end
end
save('.\MiddleMat\Supp_Figure3a.mat','cPlot','strTitle','strMethod');