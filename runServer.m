function runServer(data,options,gz,dgz,ncluster,c)
%RUNSERVER Summary of this function goes here
% Infer Type
% 2, mixture
% 6, BaNPPA, BaNPPA-NC
dataSource = data.dataSource;
options.maxIter = 200;
options.minIter = 5;
options.nEps=1e-3;
options.noise = 2;         % 2 optimize single noise
% options.approx = 0;      % 1 approximate 0 do not approximate
options.figure = 1;        % 0 no plotting 1 plotting
options.cluster = ncluster;
warning('off','all');

try
    if options.constraint == 0
        load(['Source',num2str(dataSource),'Size',num2str(data.U),'Infer',num2str(options.inferType+1),'.mat']);
    else
        load(['Source',num2str(dataSource),'Size',num2str(data.U),'Infer',num2str(options.inferType),'.mat']);
    end
catch
    mSparse = zeros(options.lenCluster,options.repeat);
    mLogTrain = zeros(options.lenCluster,options.repeat);
    mLogTest = zeros(options.lenCluster,options.repeat);
    mTime = zeros(options.lenCluster,options.repeat);
    model_save = cell(1,options.repeat);
end
i = 1;
while i<= options.repeat
    tstart = tic;
    disp([options.inferType,c,i]);
    model.gz=gz';
    model.dgz=dgz';
    model = varsgpCreateMy(data,model,options);
    if options.inferType == 6
        nBefore = inf;
        model.la.lambda = ones(options.cluster,1);
        model.la.c = 4;
        options.start = 1;
        while true
            if options.detail == 1
                [model,nIter,vTrace,vDiv] = varsgpTrainRec(data,model,options);
            else
                [model,nIter,vTrace,vDiv] = varsgpTrainMy(data,model,options);
            end
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
        if options.detail == 1
            [model,nIter] = varsgpTrainRec(data,model,options);
        else
            [model,nIter] = varsgpTrainMy(data,model,options);
        end
    end
    if ~isinf(nIter)
        % evaluate model
        mTime(c,i) = toc(tstart);
        model_save{i} = model;
        mLogTrain(c,i) = vargspTestSample(data.train,model,options);
        mLogTest(c,i) = vargspTestSample(data.test,model,options);
        if options.inferType == 6
            fprintf('Test %.4f\t Time %.4f\t Alpha %.4f\n',mLogTest(c,i),mTime(c,i),model.prior.alpha);
        else
            fprintf('Test %.4f\t Time %.4f\n',mLogTest(c,i),mTime(c,i));
        end
        modelD.f = 0;
        modelD = calPre(model);
        modelD = calIntegral(model,options,modelD,3);
        modelD = varsTransform(model,modelD,options);
        switch options.inferType
            case 2
                vamp = modelD.mAmp_2.*repmat(modelD.vfInt,1,model.U);
            case 6
                vamp = modelD.mMap.*repmat(modelD.vfInt,1,model.U);
        end
        vamp = vamp./repmat(sum(vamp,1),model.K,1);
        mSparse(c,i) = (sqrt(model.K)-mean(sum(vamp,1)./sqrt(sum(vamp.^2,1))))/(sqrt(model.K)-1);
        i = i+1;
    end
end
[~,ibest] = max(mLogTest(c,:));
if options.detail ~= 1
    if isfield(options,'constraint')
        if options.constraint == 0
            save(['Source',num2str(dataSource),'Size',num2str(data.U),'Infer',num2str(options.inferType+1),'Cluster',num2str(options.cluster),'Model.mat'],'model_save','ibest','options');
            save(['Source',num2str(dataSource),'Size',num2str(data.U),'Infer',num2str(options.inferType+1),'.mat'],'mLogTest','mLogTrain','mSparse','mTime');
        else
            save(['Source',num2str(dataSource),'Size',num2str(data.U),'Infer',num2str(options.inferType),'Cluster',num2str(options.cluster),'Model.mat'],'model_save','ibest','options');
            save(['Source',num2str(dataSource),'Size',num2str(data.U),'Infer',num2str(options.inferType),'.mat'],'mLogTest','mLogTrain','mSparse','mTime');
        end
    else
        save(['Source',num2str(dataSource),'Size',num2str(data.U),'Infer',num2str(options.inferType),'Cluster',num2str(options.cluster),'Model.mat'],'model_save','ibest','options');
        save(['Source',num2str(dataSource),'Size',num2str(data.U),'Infer',num2str(options.inferType),'.mat'],'mLogTest','mLogTrain','mSparse','mTime');
    end
end
end

