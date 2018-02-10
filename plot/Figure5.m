%% Detailed plot for one experiment Setting A.
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
clear;close all;clc;addpath(genpath('./'));
load('Figure456.mat');
nylim = 0.5;
ifontsize = 13;
ilinewidth = 2;
iplotBasis = 4; % Best four
iplot= 14;
%% Figure 4 LPPA and BNPPA on Synthetic A data set, K = 14
cNER = cell(1,2);
cBasis = cell(1,2);
%LPPA
ibest = cPlot{1,1}.Best;
cNER{1} = cPlot{1,1}.NER; 
cBasis{1} = cPlot{1,1}.fBasis{ibest};
[~,I] = sort(cPlot{1,1}.NER(:,ibest),'descend');
Iorder_LPPA = I;
%BaNPPA
ibest = cPlot{1,3}.Best;
cNER{2} = cPlot{1,3}.NER; 
cBasis{2} = cPlot{1,3}.fBasis{cPlot{1,3}.Best};
[~,I] = sort(cPlot{1,3}.NER(:,ibest),'descend');
Iorder_BaNPPA = I;
vT = cPlot{1,1}.fBasisT;


fig = figure;
set(fig,'Position',[83,109,577,318]);
subplot('Position',[0.12,0.65,0.85,0.30]);ylim([0,nylim]);
hold on;box on;grid on;set(gca,'fontsize',ifontsize);
cMean = cell(1,2);
for i = 1:2
    cMean{i} = mean(cNER{i},2);
end
hb = bar(1:iplot,[cMean{1},cMean{2}],1,'grouped');
hb(1).FaceColor = 'w';
hb(2).FaceColor = [0.5,0.5,0.5];
legend('LPPA','BaNPPA');
xlim([0,iplot+1]);
ylabel('NER');xlabel('Latent function index k');
set(gca,'YTick',[0,0.25,0.5]);
set(gca,'XTick',1:iplot);

% Finding the number of groups and the number of bars in each group
ngroups = iplot;
nbars = 2;
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    % Calculate center of each bar
    mCov = cov(cNER{i}');
    vCov = 1.96*sqrt(diag(mCov)/5);
    x = (1:ngroups)-groupwidth/2+(2*i-1)*groupwidth/(2*nbars);
    errorbar(x, cMean{i}, vCov, 'k.');
end

cmap=jet(iplotBasis);cmap(3,:)=[1,0.8,0];
subplot('Position',[0.12,0.17,0.24,0.22]);
hold on;grid on;box on;set(gca,'fontsize',ifontsize);
ylabel('E_q[f^2_k(t)]');title('LPPA');xlabel('Time t');
xlim([0 60]);set(gca,'XTick',[0,20,40]);
for ic=1:iplotBasis
    c = Iorder_LPPA(ic);
    plot(vT,cBasis{1}(c,:),'Color',cmap(ic,:),'Linewidth',ilinewidth);
end
set(gca,'YTick',[0,0.2,0.5]);

subplot('Position',[0.41,0.17,0.24,0.22]);
hold on;grid on;box on;set(gca,'fontsize',ifontsize);
title('BaNPPA');xlabel('Time t');
xlim([0 60]);set(gca,'XTick',[0,20,40]);
for ic=1:iplotBasis
    c = Iorder_BaNPPA(ic);
    plot(vT,cBasis{2}(c,:),'Color',cmap(ic,:),'Linewidth',ilinewidth);
end

subplot('Position',[0.73,0.17,0.24,0.22]);
hold on;grid on;box on;set(gca,'fontsize',ifontsize);
title('True');xlabel('Time t');ylabel('Intensity');
xlim([0 60]);ylim([0,4]);set(gca,'XTick',[0,20,40]);
fbasic = @(x)exp(-(x-5).^2/10)+exp(-(x-45).^2/10);
vX = linspace(0,60,100);
cmap = jet(4);cmap(3,:)=[1,0.8,0];
damp = [1,1,2,2];
for iu = 0:4-1
    plot(vX,fbasic(vX-iu*10)*damp(iu+1),'Color',cmap(iu+1,:),'Linewidth',ilinewidth);
end

set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(fig,'.\Result\NERBasis_SynA','-dpdf','-r0');


