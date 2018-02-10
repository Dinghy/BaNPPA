%% plot for one experiment
clear;close all;clc;
addpath(genpath('./'));
addpath(genpath('.././'));
%%
arrDataset = [1,5,2,3];
strTitle = {'Synthetic A','Synthetic B','Microblog','Citation'};
strMethod = {'BaNPPA','BaNPPA-','LPPA'};

cPlot = cell(4,3);
for idata = 1:length(arrDataset)
    dataSource = arrDataset(idata);
    switch dataSource
        case 1
            vX = [2,4,5,6,8,10,12,14];
            nSize = 200;
            arr_infer = [6,7,2];
        case 2
            vX = [4,5,6,7,8,10,12,14];
            nSize = 500;
            arr_infer = [6,7,2];
        case 3
            vX = [2,4,6,8,10,12,14];
            nSize = 600;
            arr_infer = [6,7,2];
        case 5
            vX = [4,5,6,8,10,12,14];
            nSize = 250;
            arr_infer = [6,7,2];
    end
    h = cell(1,6);
    for i = 1:length(arr_infer)
        iInfer = arr_infer(i);
        load(['Source',num2str(dataSource),'Size',num2str(nSize),'Infer',num2str(iInfer),'.mat']);
        options.repeat = 5;
        cPlot{idata,i}.Train = mLogTrain;
        cPlot{idata,i}.Test = mLogTest;
        cPlot{idata,i}.vX = vX;
    end
end
save('.\MiddleMat\Figure3a.mat','cPlot','strTitle','strMethod');
