%%
clear;close all;clc;
load('Figure10.mat');
close all;
iwidth = 1.2;

fig=figure;set(fig,'Position',[4,381,1184,236]);
vT = cPlot.vT;
subplot('Position',[0.05,0.23,0.43,0.67]);hold on;box on;grid on;
arrCenter = [5,10,20];
arrVar = [5,18,100];
cmap = [0,0,1;1,0,0;1,0.8,0];
for i=1:2
    fBasic = @(x)3/sqrt(arrVar(i))*exp(-(x-arrCenter(i)).^2/arrVar(i))+...
        3/sqrt(arrVar(i))*exp(-(x-arrCenter(i)-20).^2/arrVar(i))+...
        3/sqrt(arrVar(i))*exp(-(x-arrCenter(i)-40).^2/arrVar(i));
    plot(vT,fBasic(vT),'Color',cmap(i,:),'LineWidth',iwidth);
end
fBasic = @(x)4.5/sqrt(arrVar(3))*exp(-(x-arrCenter(3)).^2/arrVar(3))+...
        4.5/sqrt(arrVar(3))*exp(-(x-arrCenter(3)-20).^2/arrVar(3));
plot(vT,fBasic(vT),'Color',cmap(3,:));
set(gca,'FontSize',12);
title('True Intensity');xlabel('Time t');ylabel('Intensity');ylim([0 2]);
cmap = [1,0.8,0;0,0,1;1,0,0];
x = [5,10,17.5];
y = [1.2,0.5,0.34];
subplot('Position',[0.53,0.23,0.43,0.67]);hold on;box on;grid on;
strTheta = cell(1,3);
hl = cell(1,3);
for k=1:cPlot.K
    hl{k} = plot(vT,cPlot.fnRate(k,:),'Color',cmap(k,:),'LineWidth',iwidth);
    plot(vT,cPlot.fnLow(k,:),'Color',cmap(k,:),'LineStyle',':','LineWidth',iwidth);
    plot(vT,cPlot.fnHigh(k,:),'Color',cmap(k,:),'LineStyle',':','LineWidth',iwidth);
    strTheta{k} = ['a = ',num2str(exp(cPlot.logtheta{k}(1)))];
end
leg = legend([hl{1},hl{2},hl{3}],strTheta{1},strTheta{2},strTheta{3},'Location','best');
set(leg,'Position',[0.73,0.60,0.1054,0.2621]);
title('LPPA (K=3)');xlabel('Time t');
set(gca,'FontSize',12);

set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(fig,'.\Result\HyperfBasis','-dpdf','-r0');