clear;close all;clc;
addpath(genpath('./'));

% plot the basis intensity functions for Synthetic A data set.
fbasic = @(x)exp(-(x-5).^2/10)+exp(-(x-45).^2/10);
iSize = 13;
fig = figure;
set(fig,'Position',[326,142,446,286]);
subplot('Position',[0.1,0.62,0.85,0.26]);hold on;box on;grid on;
vX = linspace(0,60,100);
cmap = jet(4);cmap(3,:)=[1,0.8,0];
damp = [1,1,2,2];
for iu = 0:4-1
    plot(vX,fbasic(vX-iu*10)*damp(iu+1),'Color',cmap(iu+1,:),'Linewidth',3);
end
set(gca,'fontsize',13);ylim([0,3]);
ylabel('Intensity');
title('Synthetic A: 4 latent functions');

% plot the basis intensity functions for Synthetic B data set.
vX = linspace(0,80,100);
cmap = jet(6);
damp = [1,1,2,2,2,2];
fbasic = @(x)exp(-(x-5).^2/10)+exp(-(x-65).^2/10);
subplot('Position',[0.1,0.17,0.85,0.26]);hold on;box on;grid on;
for iu = 0:6-1
    plot(vX,fbasic(vX-iu*10)*damp(iu+1),'Color',cmap(iu+1,:),'Linewidth',3);
end
set(gca,'fontsize',13);ylim([0,3]);
ylabel('Intensity');xlabel('Time $t\in\mathcal{T}$','Interpreter','LaTeX');
title('Synthetic B: 6 latent functions');
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(fig,'.\Result\SynBasis','-dpdf','-r0');

