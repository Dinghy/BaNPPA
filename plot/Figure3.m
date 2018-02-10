%% plot for one experiment
clear;close all;clc;addpath(genpath('./'));
load('Figure3a.mat');
strTitle = {'Synthetic A','Synthetic B','Microblog','Citation'};
strMethod = {'BaNPPA','BaNPPA-NC','LPPA'};
ifontsize = 13;
ilinewidth = 1.5;
idataset = 3;

strcolor = [1,0.8,0;0,0,1;1,0,0];
strStyle = {'r','b','y'};

mTimeRec = zeros(idataset,4);
%% plot Train and Test separately
cfig = cell(1,2);
for iTrain = 0:1
    fig = figure;
    set(fig,'Position',[0,200,1327,353]);
    for idata = 1:length(strTitle)
        arrPos = [0.01+0.06*idata+(1-0.1)/(1+length(strTitle))*(idata-1),0.22,0.20,0.68];
        if iTrain == 1 && idata == 2
            arrPos(1) = arrPos(1)+0.01;
        end
        h = cell(1,6);
        for i = [2,1,3]
            vX = cPlot{idata,i}.vX;
            mLogTrain = cPlot{idata,i}.Train;
            mLogTest = cPlot{idata,i}.Test;
            options.repeat = 5;
            
            switch i
                case 3
                    mLogTrain = mLogTrain(1:length(vX),:);
                    mLogTest = mLogTest(1:length(vX),:);
                    subplot('Position',arrPos);hold on;box on;grid on;
                    
                    if iTrain == 1
                        vMean = mean(mLogTrain,2)';
                        mCov = cov(mLogTrain');
                        vCov = sqrt(diag(mCov)/options.repeat);
                        h{2*i-1} = errorbar(vX,vMean,1.96*vCov,strStyle{i},'LineWidth',ilinewidth);
                        htmp = h{2*i-1};
                        htmp.Color=strcolor(i,:);
                    else
                        vMean = mean(mLogTest,2)';
                        mCov = cov(mLogTest');
                        vCov = sqrt(diag(mCov)/options.repeat);
                        h{2*i} = errorbar(vX,vMean,1.96*vCov,strStyle{i},'LineWidth',ilinewidth);
                        htmp = h{2*i};
                        htmp.Color=strcolor(i,:);
                    end
                case {1,2}
                    subplot('Position',arrPos);hold on;box on;grid on;
                    if iTrain == 1
                        nMean = mean(mLogTrain(end,:),2)';
                        vMean = nMean*ones(size(vX));
                        nCov = sqrt(cov(mLogTrain(end,:)')/options.repeat);
                        rectangle('Position',[min(vX),nMean-1.96*nCov,max(vX)-min(vX),1.96*2*nCov],...
                            'EdgeColor',strcolor(i,:)*0.4+[1,1,1]*0.6,'FaceColor',strcolor(i,:)*0.4+[1,1,1]*0.6);
                        h{2*i-1} = plot(vX,vMean,'Color',strcolor(i,:),'LineWidth',ilinewidth);
                    else
                        nMean = mean(mLogTest(end,:),2)';
                        vMean = nMean*ones(size(vX));
                        nCov = sqrt(cov(mLogTest(end,:)')/options.repeat);
                        rectangle('Position',[min(vX),nMean-1.96*nCov,max(vX)-min(vX),1.96*2*nCov],...
                            'EdgeColor',strcolor(i,:)*0.4+[1,1,1]*0.6,'FaceColor',strcolor(i,:)*0.4+[1,1,1]*0.6);
                        h{2*i} = plot(vX,vMean,'Color',strcolor(i,:),'LineWidth',ilinewidth);
                    end
            end
        end
        set(gca,'XTick',vX);
        set(gca,'fontsize',13);
        % range modification
        if iTrain == 0 % plot test
            switch idata
                case 1
                    ylim([-7.9e3 -6.6e3]);
                    title('Synthetic A');
                    ylabel('Test Likelihood');
                case 2
                    title('Synthetic B');
                case 3
                    title('Microblog');
                case 4
                    title('Citation');  
            end
        else % plot train
            switch idata
                case 1
                    title('Synthetic A');
                    ylabel('Train Likelihood');
                case 2
                    title('Synthetic B');
                case 3
                    title('Microblog');
                case 4
                    title('Citation');
            end
        end
        xlim([min(vX)-1,max(vX)+1]);
        xlabel('Number of latent functions K');
        if idata == 1
            g = h;
        end
    end
    %%
    if iTrain == 0 % plot test
        leg2 = legend([g{6},g{4},g{2}],{'LPPA','BaNPPA-NC(K=14)','BaNPPA(K=14)'},'Orientation','horizontal');
        set(leg2,'Position',[0,0,1,0.06]);set(leg2,'Box','off');
    else % plot train
        leg2 = legend([h{5},h{3},h{1}],{'LPPA','BaNPPA-NC(K=14)','BaNPPA(K=14)'},'Orientation','horizontal');
        set(leg2,'Position',[0,0,1,0.06]);set(leg2,'Box','off');
    end
    cfig{iTrain+1}=fig;
end
%%
set(cfig{1},'Units','Inches');
pos = get(cfig{1},'Position');
set(cfig{1},'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(cfig{1},'.\Result\Comparison','-dpdf','-r0');
set(cfig{2},'Units','Inches');
pos = get(cfig{2},'Position');
set(cfig{2},'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(cfig{2},'.\Result\ComparisonTrain','-dpdf','-r0');