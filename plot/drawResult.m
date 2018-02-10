function drawResult(fig,fig2,data,model,modelD,options)
%DRAWRESULT PLOT INFERRED INTENSITY AND TRUE INTENSITY
D = size(model.Xm,2);
ifontsize = 13;
vfInt = modelD.vfInt;
nTrain=min(8,model.U);
if model.U <= 8
    arrTrain = 1:model.U;
else
    [~,I] = sort(model.NumData);
    arrTrain = [I(1:4);I(end-3:end)];
end
%% basic calculation
fnRate=zeros(model.K,size(model.plot.Xm,1));
fnLow = zeros(model.K,size(model.plot.Xm,1));
fnHigh = zeros(model.K,size(model.plot.Xm,1));
for c=1:model.K
    Knm=feval(model.cov{2}{1},model.GP.logtheta{c}(1:end-1),model.plot.Xm,model.Xm);
    vMean=Knm*(modelD.Kmm{c}\model.var.Mu{c});
    KnmInvKmm = Knm/modelD.Kmm{c};
    vVar = modelD.diagKnn{c}+sum(KnmInvKmm.*(Knm*(modelD.InvKmmSigma{c}-eye(model.m))),2);
    fnRate(c,:) = vMean.^2+vVar;
%     fnLow(c,:) = ncx2inv(0.05,1,vMean.^2./vVar).*vVar;
%     fnHigh(c,:) = ncx2inv(0.95,1,vMean.^2./vVar).*vVar;
end
% parametrization for amplitude
modelD = varsTransform(model,modelD,options);
%% basis and map result
figure(fig);
drawnow;
mMap = zeros(model.K,model.U);
switch options.inferType
    case 2
        mMap = modelD.mAmp_2;
    case 6
        mMap = modelD.mMap;
end
% mMap = mMap.*repmat(modelD.vfInt,1,model.U);
mMap = mMap./repmat(sum(mMap,1),model.K,1);
% disp(arrTrain(8));
% disp(mMap(:,arrTrain(8)));
if D==1
    subplot(model.K,2,1:2:model.K*2);
else
    if mod(model.K,2)==1
        nRow = model.K+1;
    else
        nRow = model.K;
    end
    subplot(nRow,3,1:3:nRow*3);
    arrplot = sort([2:6:nRow*3,3:6:nRow*3]);
end
imagesc(kron(mMap,ones(10)));colormap(gray);

set(gca,'XTick',round(model.U*10/4):round(model.U*10/4):model.U*10);
set(gca,'XTickLabel',cellstr(arrayfun(@(x)num2str(x),round(model.U/4):round(model.U/4):model.U,'UniformOutput',false)));
set(gca,'YTick',5:10:model.K*10);
set(gca,'YTickLabel',cellstr(arrayfun(@(x)['k',num2str(x)],1:model.K,'UniformOutput',false)));

xlabel('User');
for k=1:model.K 
    if D==1
        subplot(model.K,2,2*k);
        vT = model.plot.Xm;
        plot(vT,fnRate(k,:),'r');hold on;
        plot(vT,fnLow(k,:),'r:');
        plot(vT,fnHigh(k,:),'r:');hold off;
        ylabel(['k',num2str(k),'   '],'rot',0);
        if model.K>10
            set(gca, 'YTick', []);
        end
        if k~= model.K
            set(gca, 'XTick', []);
        end
        
        xlim([model.Tstart,model.Tend]);
        if k==model.K
            xlabel('Time');
        end
    else
        subplot(nRow,3,[arrplot(k),arrplot(k)+3]);
        imagesc(flipud(reshape(fnRate(k,:),model.plot.col,model.plot.row)));colormap(gray);colorbar;
        set(gca, 'YTick', []);
        set(gca, 'XTick', []);
        ylabel(['k',num2str(k),'   '],'rot',0);
    end
end
drawnow;
%% statistic about weight and function
cmap = jet(model.K);
mMap = zeros(model.K,model.U);
switch options.inferType
    case 2
        mMap = modelD.mAmp_2;
    case 6
        mMap = modelD.mMap;
end
mMapBiased = mMap;
mMapBiased = mMapBiased./repmat(sum(mMap,1),model.K,1);
vOccBiased = sum(mMapBiased,2)/model.U;

mMapTrue = mMap.*repmat(modelD.vfInt,1,model.U);
mMapTrue = mMapTrue./repmat(sum(mMapTrue,1),model.K,1);
vOcc = sum(mMapTrue,2)/model.U;

figure(fig2);
subplot(3,1,1);box on;
for c=1:model.K
    bar(c,vOccBiased(c),'FaceColor',cmap(c,:));
    if c == 1
        hold on;
    end
end
ylim([0,0.5]);
xlim([0,model.K+1]);
ylabel('Ratio(Biased)');
set(gca,'fontsize',ifontsize);
set(gca,'XTick',[]);
set(gca,'YTick',0:0.1:1);hold off;

subplot(3,1,2);box on;
for c=1:model.K
    bar(c,modelD.vfInt(c),'FaceColor',cmap(c,:));
    if c == 1
        hold on;
        plot(0:model.K+1,model.A*ones(size(0:model.K+1)),'r');
        plot(0:model.K+1,mean(modelD.vfInt)*ones(size(0:model.K+1)),'b');
    end
end
set(gca,'fontsize',ifontsize);
set(gca,'XTick',1:model.K);
xlim([0,model.K+1]);
set(gca,'XTick',[]);
ylabel('Integration');hold off;

subplot(3,1,3);box on;
for c=1:model.K
    bar(c,vOcc(c),'FaceColor',cmap(c,:));
    if c == 1
        hold on;
    end
end
ylim([0,0.5]);
xlim([0,model.K+1]);
ylabel('Ratio(True)');
set(gca,'fontsize',ifontsize);
set(gca,'XTick',[]);
set(gca,'YTick',0:0.1:1);hold off;


% %% fit result    
% figure(fig2); 
% for i=1:min(data.U,8)
%     u = arrTrain(i);
%     switch options.inferType
%         case 2
%             vAmp = modelD.mAmp_2(:,u);
%             fnRate_u = fnRate'*vAmp;
%         case 6
%             vAmp = mMap(:,u);
%             vAmp = model.amp.s(u)*vAmp;
%             fnRate_u = fnRate'*vAmp;
%     end
%     
%     if D == 1
%         subplot(ceil(nTrain/2),2,i);
%         stem(data.train{u}.X,ones(size(data.train{u}.X)));hold on;
%         plot(model.plot.Xm,fnRate_u,'r');
%         if isfield(data,'flambda')
%             if isfield(data,'test')
%                 vPlot = zeros(size(model.plot.Xm));
%                 for ip = 1:length(vPlot)
%                     vPlot(ip) = data.flambda{u}(model.plot.Xm(ip))/2;
%                 end
%                 plot(model.plot.Xm,vPlot,'m');
%             else
%                 plot(model.plot.Xm,data.flambda{u}(model.plot.Xm),'m');
%             end
%         end
%         xlim([data.XBeg data.XEnd]);title(num2str(u));
%     else
%         subplot(ceil(nTrain/2),4,2*i);
%         plot(model.data{u}.feature(:,1)/5,model.data{u}.feature(:,2)/5,'b+');hold on;
%         xlim([model.Tstart(1),model.Tend(1)]/5);
%         ylim([model.Tstart(2),model.Tend(2)]/5);
%         set(gca, 'YTick', []);
%         set(gca, 'XTick', []);
%         subplot(ceil(nTrain/2),4,2*i-1);
%         fnRate_image = reshape(fnRate_u,model.plot.col,model.plot.row);
%         imagesc(flipud(fnRate_image));title(num2str(u));colormap(gray);
%     end
%     hold off;
% end
drawnow;
end

