%% Varying Alpha
clear;close all;clc;
addpath(genpath('./'));
arrcolor = [1,0,0;0,0,1;1,0.8,0];
ilinewidth = 2;
fig = figure;
set(fig,'Position',[3.4,237,1217.6,380.8]);
load('Supp_Figure4.mat');
h = cell(1,4);
for dataset = 1:3
    subplot('Position',[0.08+0.32*(dataset-1),0.25,0.25,0.65]);hold on;box on;grid on;
    mLogTest = cPlot{dataset,1};
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
    
    for i=1:2
        mLogTest = cPlot{dataset,i+1};
        vMean = mean(mLogTest,2)';
        mCov = cov(mLogTest');
        vCov = sqrt(diag(mCov)/5);
        h{2+i} = errorbar(arrAlpha,vMean,1.96*vCov,'LineWidth',ilinewidth);
        h{2+i}.Color=arrcolor(1+i,:);
    end
    
    switch dataset
        case 1
            title('Synthetic C (Actual 6)');
            xlim([1 max(arrAlpha)+0.5]);
        case 2
            title('Synthetic D (Actual 8)');
            % ylim([7.6e4 7.74e4]);
            xlim([1 max(arrAlpha)+0.5]);
        case 3
            title('Synthetic E (Actual 10)');
            % ylim([1.84e5 1.87e5]);
            xlim([1 max(arrAlpha)+0.5]);
    end
    ylabel('Test Likelihood');
    set(gca,'XTick',arrAlpha);
    xlabel('Hyperparameter \alpha');
    set(gca,'fontsize',13);
end
leg1 = legend('LPPA(K=Actual)','LPPA(K=14)','BaNPPA-NC(K=14)','BaNPPA(K=14)','Orientation','horizontal');
set(leg1,'Position',[0,0,1,0.1]);
legend boxoff;
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(fig,'.\Result\Comparison_AlphaVary2','-dpdf','-r0');