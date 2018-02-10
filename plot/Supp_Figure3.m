%% plot execution time v.s. Train Likelihood/alpha
clear;close all;clc;addpath(genpath('./'));
arrcolor = [1,0,0;0,0,1;1,0.8,0];
ilinewidth = 2;

load('Supp_Figure3a.mat');
fig = figure;set(fig,'Position',[0,200,1327,353]);
for idataSet = 1:4
    subplot('Position',[0.01+0.06*idataSet+(1-0.1)/(1+length(strTitle))*(idataSet-1),0.22,0.20,0.68]);
    hold on;box on;grid on;
    for i = 1:3
        vX = cPlot{idataSet,i}.vX;
        vY = cPlot{idataSet,i}.vY;
        plot(vX,vY,'Color',arrcolor(i,:),'LineWidth',ilinewidth);
    end
    switch idataSet
        case 1
            title('Synthetic A');
            ylabel('Train Likelihood');
        case 2
            title('Synthetic B');
            ylim([-1.7e4 -0.6e4]);
        case 3
            title('Microblog');
            ylim([0 1e5]);
        case 4
            title('Citation');
            ylim([0 2e5]);
            set(gca,'XTick',[0,1500,2500]);
    end
    
    xlabel('Time (Seconds)');
    
    set(gca,'fontsize',13);
end
leg1 = legend('LPPA(K=14)','BaNPPA-NC(K=14)','BaNPPA(K=14)','Orientation','horizontal');
set(leg1,'Position',[0,0,1,0.06]);
legend boxoff;
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(fig,'.\Result\Comparison_Time','-dpdf','-r0');