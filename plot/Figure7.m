%% plot Matrix and the cluster result
clear;clc;
addpath(genpath('./'));
load('MicroblogMatrix.mat');
dataU = 500;
i=1;
model.Tstart = 3;model.Tend = 15;
iplot = 14;
ikmeans = 0;
ifontsize = 11;
fig = figure;
set(fig,'Position',[184,285,633,277]);
iplotBasis = 6;
% subplot('Position',[0.08,0.2,0.42,0.75]);
% hold on;box on;grid on;set(gca,'fontsize',ifontsize);
if ikmeans == 1
    mMap = cMap{i};
    I = kmeans(mMap',6);
    cSort = cell(1,iplot);
    for c=1:iplotBasis
        cSort{c} = find(I==c);
    end
else
    mMap = cMap{i};
    iSmallU = size(mMap,2);
    [~,I]=max(mMap,[],1);
    cSort = cell(1,iplot);
    arrSort = [];
    for c=1:iplotBasis
        cSort{c} = find(I==c);
        arrSort=[arrSort,cSort{c}];
    end
end
% mMapRe = mMap(:,arrSort);
% imagesc(flipud(kron(mMapRe,ones(10)))); colorbar;
% xlim([0 iSmallU*10]);ylim([0 iplot*10]);
% set(gca,'XTick',100:200:iSmallU*10);
% set(gca,'XTickLabel',10:20:iSmallU);
% set(gca,'YTick',5:20:iplot*10);
% set(gca,'YTickLabel',iplot:-2:1);
% xlabel('Sequence ID d');ylabel('Latent Group Index k');
iLine = 1;

vT = vX;
vColor = [0,0,1];
for ic=1:iplotBasis
    mF = cBasis{i}'*mMap(:,cSort{ic});
    size(mF,2)
    vMeanF = mean(mF,2);
    vMaxF = max(mF,[],2);
    vMinF = min(mF,[],2);
    nSpace = 1/(2+iplotBasis);
    iX = ceil(ic/2);
    iY = mod(ic,2);
    %[0.08+iY*0.48,0.63-(iX-1)*0.27,0.2,0.2]
    subplot('Position',[0.08+(1-iY)*0.49,0.74-(iX-1)*0.28,0.41,0.24]);
    box on;grid on;set(gca,'fontsize',ifontsize);hold on;
    if iLine == 0 % only area
        h = area(vT,[vMinF';vMaxF']');
        h(1).FaceColor = [1,1,1];
        h(2).FaceColor = ([1,1,1]+vColor)/2;
        plot(vT,vMeanF','Color',vColor,'Linewidth',1.5);
    else % plot lines
        for iData = 1:length(cSort{ic})
            plot(vT,mF(:,iData),'Color',([1,1,1]+vColor)/2,'Linewidth',1.5);
        end
        plot(vT,vMeanF','Color',vColor,'Linewidth',1.5);
    end
    xlim([model.Tstart,model.Tend]);ylim([0 70]);
    set(gca,'XTick',model.Tstart:2:model.Tend);
    set(gca,'YTick',[0,35,70]);
    ylabel(['Group ',num2str(ic)]);
    if iX==3
        xlabel('Time t (Hour)');
    else
        set(gca,'XTickLabel',[]);
    end
    
end