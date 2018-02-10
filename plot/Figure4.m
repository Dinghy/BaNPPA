%% Varying Alpha
clear;close all;clc;
addpath(genpath('./'));addpath(genpath('.././'));
arrcolor = [1,0,0;0,0,1;1,0.8,0];
ilinewidth = 2;
fig = figure;
set(fig,'Position',[0,200,1327,353]);
h = cell(1,4);
load('Figure3b.mat');

for idata = 1:length(strTitle)
    subplot('Position',[0.01+0.06*idata+(1-0.1)/(1+length(strTitle))*(idata-1),0.22,0.20,0.68]);hold on;box on;grid on;
    %% plot LPPA
    mLogTest = cPlot{idata,1}.Test;
    vMean = mean(mLogTest,2);
    % best
    [nMax,iMax] = max(vMean);
    vMean = nMax*ones(size(arrAlpha));
    nCov = sqrt(cov(mLogTest(iMax,:))/5);
    
    rectangle('Position',[min(arrAlpha),nMax-1.96*nCov,max(arrAlpha)-min(arrAlpha),1.96*2*nCov],...
        'EdgeColor',arrcolor(1,:)*0.5+[1,1,1]*0.5,'FaceColor',arrcolor(1,:)*0.5+[1,1,1]*0.5);
    h{1} = plot(arrAlpha,vMean,'--','Color',arrcolor(1,:),'LineWidth',ilinewidth);
    % K=14
    nMean = mean(mLogTest(end,:),2);
    vMean = nMean*ones(size(arrAlpha));
    nCov = sqrt(cov(mLogTest(end,:))/5);
    
    rectangle('Position',[min(arrAlpha),nMean-1.96*nCov,max(arrAlpha)-min(arrAlpha),1.96*2*nCov],...
        'EdgeColor',arrcolor(1,:)*0.5+[1,1,1]*0.5,'FaceColor',arrcolor(1,:)*0.5+[1,1,1]*0.5);
    h{2} = plot(arrAlpha,vMean,'Color',arrcolor(1,:),'LineWidth',ilinewidth);
    %% plot BaNPPA and BaNPPA-NC
    for i=1:2
        mLogTest = cPlot{idata,i+1}.Test;
        vMean = mean(mLogTest,2)';
        mCov = cov(mLogTest');
        vCov = sqrt(diag(mCov)/5);
        vMean = vMean(1:length(arrAlpha));
        vCov = vCov(1:length(arrAlpha));
        h{2+i} = errorbar(arrAlpha,vMean,1.96*vCov,'LineWidth',ilinewidth);
        h{2+i}.Color=arrcolor(1+i,:);
    end
    
    
    switch idata
        case 1
            title('Synthetic A');
            xlim([0 9]);
            ylabel('Test Likelihood');
        case 3
            title('Microblog');
            ylim([7.6e4 7.74e4]);
            xlim([0 9]);
        case 4
            title('Citation');
            ylim([1.84e5 1.87e5]);
            xlim([0 9]);
        case 2
            title('Synthetic B');
            xlim([0 9]);
    end
    
    set(gca,'XTick',arrAlpha);
    xlabel('Hyperparameter \alpha');
    set(gca,'fontsize',13);
end
leg1 = legend('LPPA(Best)','LPPA(K=14)','BaNPPA-NC(K=14)','BaNPPA(K=14)','Orientation','horizontal');
set(leg1,'Position',[0,0,1,0.06]);
legend boxoff;
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(fig,'.\Result\Comparison_AlphaVary','-dpdf','-r0');