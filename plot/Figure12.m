%%
clear;close all;clc;
load('Figure12.mat');

iwidth = 1.2;
cmap = [0,0,1;1,0,0;1,0.8,0];
fig=figure;
set(fig,'Position',[4,281,1184,296]);
vT = cPlot.vT;

subplot('Position',[0.05,0.20,0.43,0.7]);hold on;box on;grid on;
hl = cell(1,cPlot.K);
for c=1:cPlot.TrueK
    hl{c} = plot(vT,cPlot.TrueMean{c},'Color',cmap(c,:),'LineWidth',iwidth);
    plot(vT,cPlot.TrueHigh{c},'Color',cmap(c,:),'LineStyle',':','LineWidth',iwidth);
    plot(vT,cPlot.TrueLow{c},'Color',cmap(c,:),'LineStyle',':','LineWidth',iwidth);
end
legend([hl{1},hl{2}],'a = 3','a = 15','Location','best');
set(gca,'FontSize',12);
title('True Intensity');xlabel('Time t');ylabel('Intensity');

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
legend([hl{1},hl{2}],strKernel{1},strKernel{2},'Location','best');
title('LPPA (K=2)');xlabel('Time t');
set(gca,'FontSize',12);

set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(fig,'.\Result\HyperGP2fBasis','-dpdf','-r0');
