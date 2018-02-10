function runServerAlpha(data,options,gz,dgz,ncluster,c)
%RUNSERVER Summary of this function goes here
% Infer Type
% 2, mixture
% 5, mixture+Gamma Process Prior+Regularization;
dataSource = data.dataSource;
options.maxIter = 200;
options.minIter = 5;
options.nEps=1e-3;
options.noise = 2;        % optimize single noise
options.approx = 0;       % 1 approximate 0 do not approximate
options.figure = 1;       % 0 no plotting
options.cluster = ncluster;
warning('off','all');


arrAlpha = [1.1,2,4,6,8,10];
iInfer = options.inferType;
if options.inferType == 6 && isfield(options,'constraint') && options.constraint == 0
    iInfer = iInfer + 1;
end
    
try
    load(['Alpha',num2str(dataSource),'Size',num2str(data.U),'Infer',num2str(iInfer),'.mat']);
catch
    mLogTrain = zeros(length(arrAlpha),options.repeat);
    mLogTest = zeros(length(arrAlpha),options.repeat);
    mTime = zeros(length(arrAlpha),options.repeat);
    model_save = cell(length(arrAlpha),options.repeat);
end

for j = 1:length(arrAlpha)
    options.alpha = arrAlpha(j);
    i = 1;
    while i<= options.repeat
        tstart = tic;
        disp([options.inferType,c,i]);
        model.gz=gz';
        model.dgz=dgz';
        model = varsgpCreateMy(data,model,options);
        if options.inferType >= 5
            nBefore = inf;
            model.la.lambda = ones(options.cluster,1);
            model.la.c = 4;
            options.start = 1;
            while true
                [model,nIter,vTrace,vDiv] = varsgpTrainRec(data,model,options);
                if options.constraint == 0 || isinf(nIter)
                    break
                end
                if sum(abs(nBefore-vTrace(nIter)))<= nIter*1e-3*abs(vTrace(nIter))
                    break
                else
                    options.start = 0;
                    nBefore = vTrace(nIter);
                    model.la.lambda = model.la.lambda + model.la.c*vDiv;
                    model.la.c = model.la.c*4;
                end
            end
        else
            options.start = 1;
            [model,nIter] = varsgpTrainRec(data,model,options);
        end
        if ~isinf(nIter)
            % evaluate model
            mTime(j,i) = toc(tstart);
            model_save{i} = model;
            mLogTrain(j,i) = vargspTestSample(data.train,model,options);
            mLogTest(j,i) = vargspTestSample(data.test,model,options);
            
            fprintf('Test %.4f\t Time %.4f\n',mLogTest(j,i),mTime(j,i));
            modelD.f = 0;
            modelD = calPre(model);
            modelD = calIntegral(model,options,modelD,3);
            modelD = varsTransform(model,modelD,options);
            switch options.inferType
                case 2
                    vamp = modelD.mAmp_2.*repmat(modelD.vfInt,1,model.U);
                case {5,6}
                    vamp = modelD.mMap.*repmat(modelD.vfInt,1,model.U);
            end
            vamp = vamp./repmat(sum(vamp,1),model.K,1);
            i = i+1;
        end
    end
    
    [~,ibest] = max(mLogTest(j,:));
    iInfer = options.inferType;
    if options.inferType == 6 && isfield(options,'constraint') && options.constraint == 0
        iInfer = iInfer + 1;
    end
    save(['Alpha',num2str(j),'Source',num2str(dataSource),'Size',num2str(data.U),'Infer',num2str(iInfer),'Cluster',num2str(options.cluster),'Model.mat'],'model_save','ibest','options');
    save(['Alpha',num2str(dataSource),'Size',num2str(data.U),'Infer',num2str(iInfer),'.mat'],'mLogTest','mLogTrain','mTime');
    
end
end

