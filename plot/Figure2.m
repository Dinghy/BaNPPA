clear;close all;clc;
addpath(genpath('./'));
iplot = 14;
ifontsize = 16;
ilinewidth = 2;
iInfer = 6;

load('Figure2.mat');
% Microblog data set
fig=figure;
for i=1:2
    switch i
        case 1
            arrPos1 = [0.12,12/18,0.83,5/18];
            arrPos2 = [0.12,11/18,0.83,1/18-0.01];
        case 2
            arrPos1 = [0.12,4/18,0.83,5/18];
            arrPos2 = [0.12,3/18,0.83,1/18-0.01];
    end
    % intensity plot
    subplot('Position',arrPos1);hold on;box on;grid on;
    plot(cPlot{1,i}.vT,cPlot{1,i}.fnRate,'k','Linewidth',ilinewidth);
    legend('BaNPPA','location','best');legend boxoff;
    xlim([3,15]);
    ylabel('Intensity');
    set(gca,'XTick',3:2:15);
    set(gca,'XTickLabel',[]);
    set(gca,'fontsize',ifontsize);
    switch i
        case 1
            text(3.5,35,'Active Hour','FontSize',ifontsize);
            set(gca,'YTick',0:20:40);
        case 2
            text(3.5,45,'Inactive Hour','FontSize',ifontsize);
            set(gca,'YTick',0:25:50);
    end
    % bar plot
    subplot('Position',arrPos2);box off;hold on;
    for ibar=1:length(cPlot{1,i}.X)
        iX=cPlot{1,i}.X(ibar);
        line([iX,iX],[-1,1],'Color','k','LineStyle','-');
    end
    xlim([3,15]);
    set(gca,'XTick',3:2:15);
    
    set(gca,'XTickLabel',cellstr(arrayfun(@(x)[num2str(mod(cPlot{1,i}.hour+x,24)),':00'],3:2:15,'UniformOutput',false)));
    set(gca,'YTickLabel',[]);
    set(gca,'YTick',[]);
    set(gca,'fontsize',ifontsize);
    if i== 2
        xlabel('Time t (Hour)');
    end
end
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(fig,'.\Result\EgMicroblog','-dpdf','-r0');

%% Citation data set
fig = figure;
for i=1:2
    switch i
        case 1
            arrPos1 = [0.12,12/18,0.83,5/18];
            arrPos2 = [0.12,11/18,0.83,1/18-0.01];
        case 2
            arrPos1 = [0.12,4/18,0.83,5/18];
            arrPos2 = [0.12,3/18,0.83,1/18-0.01];
    end
    % intensity plot
    subplot('Position',arrPos1);hold on;box on;grid on;
    plot(cPlot{2,i}.vT,cPlot{2,i}.fnRate,'k','Linewidth',ilinewidth);
    legend('BaNPPA','location','best');legend boxoff;
    xlim([0 21]);
    ylabel('Intensity');
    set(gca,'XTick',0:2:21);
    set(gca,'XTickLabel',[]);
    set(gca,'fontsize',ifontsize);
    switch i
        case 1
            set(gca,'YTick',0:25:50);
        case 2
            set(gca,'YTick',0:20:40);
    end
    % bar plot
    subplot('Position',arrPos2);box off;hold on;
    for ibar=1:length(cPlot{2,i}.X)
        iX=cPlot{2,i}.X(ibar);
        line([iX,iX],[-1,1],'Color','k','LineStyle','-');
    end
    xlim([0,21]);
    set(gca,'XTick',0:2:21);
    
    set(gca,'YTickLabel',[]);
    set(gca,'YTick',[]);
    set(gca,'fontsize',ifontsize);
    if i==2
        xlabel('Time t (Year)');
    end
end
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(fig,'.\Result\EgCitation','-dpdf','-r0');
