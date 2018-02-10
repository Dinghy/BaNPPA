%% Varying Alpha
clear;close all;clc;
addpath(genpath('./'));
arrAlpha = [1.1,2,3,4,5];

cPlot = cell(3,3);
strTitle = {'Synthetic C','Synthetic D','Synthetic E'};
strMethod = {'LPPA','BaNPPA-','BaNPPA'};

for dataset = 1:3
    iSource = dataset+3;
    switch dataset
        case 1
            U = 200;
            arrInfer = [7,6];
        case 2
            U = 200;
            arrInfer = [7,6];
        case 3
            U = 200;
            arrInfer = [7,6];
    end
    load(['.././ServerAdditional/Source',num2str(iSource),'Size',num2str(U),'Infer2.mat']);
    cPlot{dataset,1} = mLogTest;
    
    for i=1:length(arrInfer)
        iInfer = arrInfer(i);
        load(['.././ServerAdditional/Alpha',num2str(iSource),'Size',num2str(U),'Infer',num2str(iInfer),'.mat']);
        cPlot{dataset,i+1} = mLogTest;
    end
end
save('./MiddleMat/Supp_Figure4.mat','arrAlpha','cPlot');
