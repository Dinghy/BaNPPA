%% Figure 6 LPPA and BaNPPA on Microblog data set, K = 14
clear;close all;clc;addpath(genpath('./'));
load('Figure456.mat');
nylim = 0.4;
ifontsize = 13;
ilinewidth = 2;
iplotBasis = 4; % Best four
iplot= 14;
ibarwidth = 0.8;
strMethod = {'LPPA','BaNPPA-NC','BaNPPA'};
iMethod = 3;
%% get the contents
idataset = 4;
switch idataset
    case 1
        iSample = 2;
        idataU = 200;
    case 2
        iSample = 2;
        idataU = 250;
    case 3
        iSample = 5;
        idataU = 500;
    case 4
        iSample = 6;
        idataU = 600;
end
cNER = cell(1,iMethod);
cMap = cell(1,iMethod);
cSort = cell(iMethod,iplot);
for i = 1:iMethod
    ibest = cPlot{idataset,i}.Best;
    cNER{i} = cPlot{idataset,i}.NER(:,ibest);
    cMap{i} = cPlot{idataset,i}.Map{ibest}(:,1:iSample:idataU);
    [~,I]=max(cMap{i},[],1);
    arrSort = [];
    for c=1:iplot
        cSort{i,c} = find(I==c);
        vWeight = cMap{i}(c,cSort{i,c});
        [~,Imid] = sort(vWeight,'descend');
        cSort{i,c} = cSort{i,c}(Imid);
        arrSort=[arrSort,cSort{i,c}];
    end
    mMapRe = cMap{i}(:,arrSort);
    cMap{i} = mMapRe;
end
vX = 1:iplot;
fig = figure;
for i=1:iMethod
    set(fig,'Position',[83,109,920,380]);
    subplot('Position',[0.08+(i-1)*0.29,0.70,0.25,0.23]);ylim([0,nylim]);xlim([0.5,14.5]);
    hold on;box on;grid on;set(gca,'fontsize',ifontsize);
    if i==1
        ylabel('NER');
    end
    title(strMethod{i});set(gca,'XTick',1:14);set(gca,'XTickLabel',[]);
    bar(vX,cNER{i},ibarwidth,'FaceColor','w');
    subplot('Position',[0.08+(i-1)*0.29,0.15,0.25,0.53]);box on;
    imagesc(cMap{i}');xlabel('Latent function index k');
    if i==1
        ylabel('Sequence ID d');
    end
    set(gca,'fontsize',ifontsize);
end
colorbar('Position',[0.92,0.15,0.03,0.53]);

set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(fig,['.\Result\NERMatrix_Full',num2str(idataset)],'-dpdf','-r0');