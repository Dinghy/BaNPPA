clear;close all;clc;
addpath(genpath('./'));

load('Figure9.mat');

fig = figure;set(fig,'Position',[9,213,1237,304]);
%%
arrInfer = [2,7,6];
strTitle = cPlot.strTitle;
strMethod = cPlot.strMethod;
caxis([0 0.5]);
cmap = colormap(jet);
colmin=1;colmax=0;
for idataset = 1:4
    for i=1:3
        colmin = min(colmin,min(cPlot.NER{idataset,i}));
        colmax = max(colmax,max(cPlot.NER{idataset,i}));
    end
end
colrange = colmax-colmin;
for idataset = 1:4
    subplot('Position',[0.05+(idataset-1)*0.23,0.10,0.2,0.8]);xlim([0 6]);
    imax = -inf;
    for i=1:3
        imax = max([imax,max(cPlot.Hyper{idataset,i})]);
    end
    ylim([0 imax+0.5]);hold on;grid on;box on;
    % line([0,6],cPlot.nStartHyper(idataset)*ones(1,2),'Color','r');
    for i=1:3
        vX = [(i-1)*2,2*i];
        for ik=1:14
            nIndex = floor((cPlot.NER{idataset,i}(ik)-colmin)/colrange*63)+1;
            line(vX,cPlot.Hyper{idataset,i}(ik)*ones(1,2),'Color',cmap(nIndex,:),'LineWidth',1);
        end
    end
    title(strTitle{idataset});
    if idataset == 1
        ylabel('Kernel length scale');
    end
    set(gca,'XTick',[1,3,5]);
    set(gca,'XTickLabel',{'LPPA','BaNPPA-NC','BaNPPA'});
    set(gca,'FontSize',12);
end

h=colorbar('Position',[0.95,0.10,0.02,0.8]);
caxis([colmin colmax]);
title(h,'NER');
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(fig,'.\Result\Hyper','-dpdf','-r0');
