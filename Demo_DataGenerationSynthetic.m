% clear;close all;clc;addpath(genpath('./'));
% %% Test 1A
% dataSource = 1;
% 
% data.Individual = 4;
% data.U = 400;
% data.XBeg = 0;
% data.XEnd = 60;
% data.gamma = [1.2,1.0,0.8,0.6];
% data.flambda = cell(1,data.U);
% data.train = cell(1,data.U);
% data.test = cell(1,data.U);
% data.amp = zeros(data.Individual,data.U);
% u = 1;
% fbasic = @(x)exp(-(x-5).^2/30)+exp(-(x-45).^2/30);
% while u<= data.U
%     disp(u);
%     data.Amp(:,u)= gamrnd(data.gamma,1);
%     data.Amp(:,u)= gamrnd(2,3)*data.Amp(:,u)/sum(data.Amp(:,u));
%     data.flambda{u} = @(x)fbasic(repmat(-10*(0:data.Individual-1),length(x),1)+repmat(x,1,data.Individual))*data.Amp(:,u); % x is column
%     X = sampleIPP(data.flambda{u},data.XBeg,data.XEnd,max(sum(data.Amp(:,u)))+2,1)';
%     if length(X)>=10
%         index = rand(size(X))>=0.5;
%         data.train{u}.X = X(index);
%         data.test{u}.X = X(~index);
%         u=u+1;
%     end
% end
% figure;imshow(kron(data.Amp./repmat(sum(data.Amp,1),data.Individual,1),ones(10)));
% save('dataSource1A.mat','data');
% options.m = 20;

%% Settings
clear;close all;clc;addpath(genpath('./'));
iDataset = 2;
data.XBeg = 0;
strSave = {'A','B','C','D','E'};
switch iDataset
    case 1
        data.U = 200;
        data.Individual = 4;
        data.gamma = [1.2,1.0,0.8,0.6];
        data.XEnd = 60;
        fbasic = @(x)exp(-(x+5).^2/10)+exp(-(x-35).^2/10);
        arrPeak = [1,1,2,2];
    case 2
        data.U = 250;
        data.Individual = 6;
        data.gamma = [1.2,1.0,0.8,0.6,0.5,0.5];
        data.XEnd = 80;
        fbasic = @(x)exp(-(x+5).^2/10)+exp(-(x-55).^2/10);
        arrPeak = [1,1,2,2,2,2];
    case 3
        data.U = 250;
        data.Individual = 6;
        data.gamma = [0.8,0.4,0.2,0.2,0.2,0.2];
        data.XEnd = 60;
        fbasic = @(x)exp(-(x+5).^2/10)+exp(-(x-55).^2/10);
        arrPeak = ones(1,6);
    case 4
        data.U = 250;
        data.Individual = 8;
        data.gamma = [0.8,0.4,0.4,0.2,0.2,0.2,0.2,0.2];
        data.XEnd = 80;
        fbasic = @(x)exp(-(x+5).^2/10);
        arrPeak = ones(1,8);
    case 5
        data.U = 250;
        data.Individual = 10;
        data.gamma = [0.8,0.6,0.4,0.4,0.4,0.2,0.2,0.2,0.1,0.1];
        data.XEnd = 100;
        fbasic = @(x)exp(-(x+5).^2/10);
        arrPeak = ones(1,10);
end
data.flambda = cell(1,data.U);
data.train = cell(1,data.U);
data.test = cell(1,data.U);
data.amp = zeros(data.Individual,data.U);
u = 1;
cmap = jet(data.Individual);
figure;hold on;grid on;
vX = linspace(data.XBeg,data.XEnd,30*data.Individual);

for c=1:data.Individual
    subplot(data.Individual,1,c);
    y = arrPeak(c)*fbasic(vX-10*c);
    plot(vX,y,'Color',cmap(c,:));
end
%% Generate points
while u<= data.U
    disp(u);
    data.Amp(:,u)= gamrnd(data.gamma,1);
    data.Amp(:,u)= gamrnd(2,3)*data.Amp(:,u)/sum(data.Amp(:,u));
    data.flambda{u} = @(x)fbasic(repmat(-10*(1:data.Individual),length(x),1)+repmat(x,1,data.Individual))*(data.Amp(:,u).*arrPeak'); 
    nfMax = sum(data.Amp(:,u).*arrPeak')+2;
    X = sampleIPP(data.flambda{u},data.XBeg,data.XEnd,nfMax,1)';
    if length(X)>=10
        index = rand(size(X))>=0.5;
        data.train{u}.X = X(index);
        data.test{u}.X = X(~index);
        u=u+1;
    end
    clc;
end
figure;imshow(kron(data.Amp./repmat(sum(data.Amp,1),data.Individual,1),ones(10)));
save(['dataSource1',strSave{iDataset},'.mat'],'data');



