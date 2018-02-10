%% Figure 7 
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
clear;close all;clc;
addpath(genpath('./'));
load('Figure7.mat');
cNER = cell(1,2);
cUNER = cell(1,2);
cBasis = cell(1,2);
cFInt = cell(1,2);
for i=1:2
    cNER{i} = cPlot{i}.NER;
    cUNER{i} = cPlot{i}.UNER;
    cFInt{i} = cPlot{i}.fInt;
    cBasis{i} = cPlot{i}.fnRate;
end
vT = cPlot{1}.vT;
%% plot
iplot = 14;
ifontsize = 13;
ilinewidth = 2;
iplotBasis = 6;
fig = figure;
set(fig,'Position',[122,35,576,386]);
%% UNER
subplot('Position',[0.12,0.76,0.85,0.21]);
hold on;box on;grid on;set(gca,'fontsize',ifontsize);
hb = bar(1:iplot,[mean(cUNER{1},2),mean(cUNER{2},2)],1,'grouped');
hb(1).FaceColor = 'w';
hb(2).FaceColor = [0.5,0.5,0.5];
legend('BaNPPA','BaNPPA-NC','Location','best');
xlim([0,iplot+1]);ylim([0 0.3]);
ylabel('UNER');
set(gca,'YTick',[0,0.1,0.3]);
set(gca,'XTick',1:iplot);
% Finding the number of groups and the number of bars in each group
ngroups = 14;
nbars = 2;
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    mCov = cov(cUNER{i}');
    vCov = 1.96*sqrt(diag(mCov)/5);
    x= (1:ngroups)-groupwidth/2+(2*i-1)*groupwidth/(2*nbars);
    errorbar(x, mean(cUNER{i},2), vCov, 'k.');
end

%% NER
subplot('Position',[0.12,0.48,0.85,0.21]);
hold on;box on;grid on;set(gca,'fontsize',ifontsize);

hb = bar(1:iplot,[mean(cFInt{1},2),mean(cFInt{2},2)],1,'grouped');
hb(1).FaceColor = 'w';
hb(2).FaceColor = [0.5,0.5,0.5];
legend('BaNPPA','BaNPPA-NC','Location','best');
xlim([0,iplot+1]);ylim([0 200]);
ylabel('Volume');xlabel('Latent function index k');
set(gca,'YTick',[0,100,200]);
set(gca,'XTick',1:iplot);
% Finding the number of groups and the number of bars in each group
ngroups = 14;
nbars = 2;
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar

for i = 1:nbars
    % Calculate center of each bar
    mCov = cov(cFInt{i}');
    vCov = 1.96*sqrt(diag(mCov)/5);
    x= (1:ngroups)-groupwidth/2+(2*i-1)*groupwidth/(2*nbars);
    errorbar(x, mean(cFInt{i},2), vCov, 'k.');
end

cmap=jet(iplotBasis);cmap(5,:)=[1,0.8,0];
subplot('Position',[0.12,0.13,0.24,0.15]);hold on;grid on;box on;set(gca,'fontsize',ifontsize);
ylabel('E_qf^2_k(t)');title('BaNPPA-NC');xlabel('Time t');
for ic=1:iplotBasis
    if iplotBasis == 1
        plot(vT,cBasis{2}(ic,:),'k','Linewidth',ilinewidth);
    else
        plot(vT,cBasis{2}(ic,:),'Color',cmap(ic,:),'Linewidth',ilinewidth);
    end
end
xlim([0,80]);set(gca,'XTick',[0,40,80]);

subplot('Position',[0.41,0.13,0.24,0.15]);hold on;grid on;box on;set(gca,'fontsize',ifontsize);
title('BaNPPA');xlabel('Time t');
for ic=1:iplotBasis
    if iplotBasis == 1
        plot(vT,cBasis{1}(ic,:),'k','Linewidth',ilinewidth);
    else
        plot(vT,cBasis{1}(ic,:),'Color',cmap(ic,:),'Linewidth',ilinewidth);
    end
end
xlim([0,80]);set(gca,'XTick',[0,40,80]);


subplot('Position',[0.73,0.13,0.24,0.15]);hold on;grid on;box on;set(gca,'fontsize',ifontsize);
title('True');ylabel('Intensity');xlabel('Time t');
cmap = jet(6);cmap(5,:)=[1,0.8,0];
damp = [1,1,2,2,2,2];vX = linspace(0,80,100);
fbasic = @(x)exp(-(x-5).^2/10)+exp(-(x-65).^2/10);
for iu = 0:6-1
    plot(vX,fbasic(vX-iu*10)*damp(iu+1),'Color',cmap(iu+1,:),'Linewidth',ilinewidth);
end
xlim([0,80]);ylim([0 4]);set(gca,'XTick',[0,40,80]);


set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(fig,'.\Result\UNERNER_SynB','-dpdf','-r0');