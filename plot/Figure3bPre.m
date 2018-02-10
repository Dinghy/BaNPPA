%% Varying Alpha
clear;close all;clc;
addpath(genpath('./'));addpath(genpath('.././'));

arrAlpha = [1.1,2,4,6,8];
arrDataset = [1,5,2,3];
cPlot = cell(4,3);
strTitle = {'Synthetic A','Synthetic B','Microblog','Citation'};
strMethod = {'LPPA','BaNPPA-','BaNPPA'};

for idata = 1:length(arrDataset)
    dataset = arrDataset(idata);
    switch dataset
        case 1
            U = 200;
            arrInfer = [7,6];
        case 2
            U = 500;
            arrInfer = [7,6];
        case 3
            U = 600;
            arrInfer = [7,6];
        case 5
            U = 250;
            arrInfer = [7,6];
    end
    load(['Source',num2str(dataset),'Size',num2str(U),'Infer2.mat']);
    cPlot{idata,1}.Test = mLogTest;
    
    for i=1:length(arrInfer)
        iInfer = arrInfer(i);
        load(['Alpha',num2str(dataset),'Size',num2str(U),'Infer',num2str(iInfer),'.mat']);
        cPlot{idata,i+1}.Test = mLogTest;
    end
    
end
save('.\MiddleMat\Figure3b.mat','cPlot','strTitle','strMethod','arrAlpha');
