%%
clear;close all;clc;
load('Figure11.mat');

iwidth = 1.2;
cmap = [0,0,1;1,0,0;1,0.8,0];
fig=figure;
set(fig,'Position',[4,281,1184,296]);
vT = cPlot.vT;

subplot('Position',[0.05,0.20,0.43,0.7]);hold on;box on;grid on;
arrCenter = [5,10,20];
arrVar = [5,18,100];
hl = cell(1,cPlot.K);
for c=1:cPlot.TrueK
    hl{c} = plot(vT,cPlot.TrueMean{c},'Color',cmap(c,:),'LineWidth',iwidth);
    plot(vT,cPlot.TrueHigh{c},'Color',cmap(c,:),'LineStyle',':','LineWidth',iwidth);
    plot(vT,cPlot.TrueLow{c},'Color',cmap(c,:),'LineStyle',':','LineWidth',iwidth);
end
legend([hl{1},hl{2},hl{3}],'a = 3','a = 6','a = 24','Location','best');
set(gca,'FontSize',12);
title('True Intensity');xlabel('Time t');ylabel('Intensity');

x = [5,10,17.5];
y = [1.2,0.5,0.34];
subplot('Position',[0.53,0.20,0.43,0.7]);hold on;box on;grid on;
hl = cell(1,cPlot.K);
strKernel = cell(1,3);
cmapcol = [1,2,3];
for k=1:cPlot.K
    ck = cmapcol(k);
    hl{k} = plot(vT,cPlot.fnRate(k,:),'Color',cmap(ck,:),'LineWidth',iwidth);
    plot(vT,cPlot.fnLow(k,:),'Color',cmap(ck,:),'LineStyle',':','LineWidth',iwidth);
    plot(vT,cPlot.fnHigh(k,:),'Color',cmap(ck,:),'LineStyle',':','LineWidth',iwidth);
    strKernel{k} = ['a = ',num2str(exp(cPlot.logtheta{k}(1)))];
end
legend([hl{1},hl{2},hl{3}],strKernel{1},strKernel{2},strKernel{3},'Location','best');
title('LPPA (K=3)');xlabel('Time t');
set(gca,'FontSize',12);

set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(fig,'.\Result\HyperGPfBasis','-dpdf','-r0');
