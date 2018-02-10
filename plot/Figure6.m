%% Figure 6 LPPA and BaNPPA on Microblog data set, K = 14
clear;close all;clc;addpath(genpath('./'));
load('Figure456.mat');
nylim = 0.4;
ifontsize = 13;
ilinewidth = 2;
iplotBasis = 4; % Best four
iplot= 14;
ibarwidth = 0.8;
%% get the contents
cNER = cell(1,2);
cMap = cell(1,2);
cSort = cell(2,iplot);
%LPPA
ibest = cPlot{3,1}.Best;
cNER{1} = cPlot{3,1}.NER(:,ibest); 
cMap{1} = cPlot{3,1}.Map{ibest}(:,1:5:500);
[~,I]=max(cMap{1},[],1);
arrSort = [];
for c=1:iplot
    cSort{1,c} = find(I==c);
    vWeight = cMap{1}(c,cSort{1,c});
    [~,Imid] = sort(vWeight,'descend');
    cSort{1,c} = cSort{1,c}(Imid);
    arrSort=[arrSort,cSort{1,c}];
end
mMapRe = cMap{1}(:,arrSort);
cMap{1} = mMapRe;
%BaNPPA
ibest = cPlot{3,3}.Best;
cNER{2} = cPlot{3,3}.NER(:,ibest); 
cMap{2} = cPlot{3,3}.Map{ibest}(:,1:5:500);
[~,I]=max(cMap{2},[],1);
arrSort = [];
for c=1:iplot
    cSort{2,c} = find(I==c);
    vWeight = cMap{2}(c,cSort{2,c});
    [~,Imid] = sort(vWeight,'descend');
    cSort{2,c} = cSort{2,c}(Imid);
    arrSort=[arrSort,cSort{2,c}];
end
mMapRe = cMap{2}(:,arrSort);
cMap{2} = mMapRe;
%% plot 
close all;
vX = 1:iplot;
fig = figure;
set(fig,'Position',[83,109,620,480]);
subplot('Position',[0.1,0.70,0.37,0.23]);ylim([0,nylim]);xlim([0.5,14.5]);
hold on;box on;grid on;set(gca,'fontsize',ifontsize);
ylabel('NER');title('LPPA');set(gca,'XTick',1:14);set(gca,'XTickLabel',[]);
bar(vX,cNER{1},ibarwidth,'FaceColor','w');
%
subplot('Position',[0.53,0.70,0.37,0.23]);ylim([0,nylim]);xlim([0.5,14.5]);
hold on;box on;grid on;set(gca,'fontsize',ifontsize);
title('BaNPPA');set(gca,'XTick',1:14);set(gca,'XTickLabel',[]);
bar(vX,cNER{2},ibarwidth,'FaceColor','w');
%
subplot('Position',[0.1,0.15,0.37,0.53]);box on;
imagesc(cMap{1}');xlabel('Latent function index k');ylabel('Sequence ID d');set(gca,'fontsize',ifontsize);
%
subplot('Position',[0.53,0.15,0.37,0.53]);box on;
imagesc(cMap{2}');xlabel('Latent function index k');set(gca,'fontsize',ifontsize);
colorbar('Position',[0.92,0.15,0.03,0.53]);
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(fig,'.\Result\NERMatrix_Micro','-dpdf','-r0');