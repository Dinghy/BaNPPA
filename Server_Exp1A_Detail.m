%% synthetic
clear; close all;clc;
addpath(genpath('./'));
%%
load('gzdgz.mat');
load('dataSource1A.mat');
data.dataSource = 1;
options.m = 18;
arr_cluster = [2,4,5,6,8,10,12,14];

options.repeat = 1;
options.lenCluster = length(arr_cluster);
options.detail = 1;

for i=length(arr_cluster)
    ncluster = arr_cluster(i);
    options.AmpIter = 5;
    options.inferType = 2;
    options.constraint = 1;
    runServer(data,options,gz,dgz,ncluster,i);
    
    if i>=length(arr_cluster)
        options.AmpIter = 10;
        options.inferType = 6;
        options.constraint = 1;     % add constraint
        runServer(data,options,gz,dgz,ncluster,i);

        options.inferType = 6;
        options.constraint = 0;     % no constraint
        runServer(data,options,gz,dgz,ncluster,i);
    end
end