%% synthetic
clear; close all;clc;
addpath(genpath('./'));
%%
load('gzdgz.mat');
load('dataSource1A.mat');
data.dataSource = 1;
options.m = 18;
arr_cluster = [2,4,5,6,8,10,12,14];

options.repeat = 5;
options.lenCluster = length(arr_cluster);
options.detail = 0;

i=length(arr_cluster);
ncluster = arr_cluster(i);
options.AmpIter = 10;
options.inferType = 6;
options.constraint = 1;     % add constraint
runServerAlpha(data,options,gz,dgz,ncluster,i);

options.inferType = 6;
options.constraint = 0;     % no constraint
runServerAlpha(data,options,gz,dgz,ncluster,i);
