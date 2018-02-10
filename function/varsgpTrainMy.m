function [model, nIter, margLogL,vDiv] = varsgpTrainMy(data,model,options)
% VARSGPCONTRAIN Summary of this function goes here
%   Detailed explanation goes here
debug = 0;
nIter = 0;
vEM=zeros(1,options.maxIter);
%% Debug place
if debug==1
    % Extract parameter
    options.debug = 0;
    DebugTest(model,options)
    margLogL = vEM;
    return
end
%% Figures or not
if options.figure == 1
    fig = figure; fig2 = figure;
end
%% Main Iteration
mSave = zeros(options.maxIter,8);
nEMPre = 0;
[LB,UB,LBHyper,UBHyper] = optConstraint(model,options);
fprintf('LB %.4f\tUB %.4f\n',LBHyper(end),UBHyper(end));
for im = 1:options.maxIter
    options.debug=0;
    options.print=0;
    model.derivative = 0;
    modelD = calPre(model);
    modelD = calIntegral(model,options,modelD,3);
    fprintf('Iterartion M:\t%d\n',im);
    %% expectation
    if im > 3 || options.start == 0 % do not optimize amplitude for several rounds
        switch options.inferType
            case 2
                for u = 1:model.U
                    V = extractAmp(model,options,u);
                    % optimize object function
                    myObj = @(VV)varsGpPpBoundAmp(VV,model,modelD,options,u);
                    V = minConf_TMP(myObj,V,zeros(size(V)),inf(size(V)),struct('verbose',0,'maxIter',options.AmpIter));
                    model = returnAmp(V,model,options,u);
                end
            case 6
                % optimize theta
                for u = 1:model.U
                    V = extractAmp(model,options,u);
                    myObj = @(VV)varsGpPpBoundAmp(VV,model,modelD,options,u);
                    V = minConf_TMP(myObj,V,zeros(size(V))+options.nEps,inf(size(V)),struct('verbose',0,'maxIter',options.AmpIter));
                    model = returnAmp(V,model,options,u);
                end
                % Optimize s
                modelD = varsTransform(model,modelD,options);
                for u=1:model.U
                    LD = model.prior.b+sum(modelD.vfInt.*modelD.mMap(:,u));
                    model.amp.s(u) = (model.NumData(u)+model.prior.a-1)/LD;
                end
        end
    end
    if options.figure == 1
        drawResult(fig, fig2, data, model, modelD, options);
    end
    % Optimize the Mu and Sigma
    modelD = varsTransform(model,modelD,options);
    V = extractVariationParamsMy(model,options);
    myObj = @(VV)varsGpPpBoundVarParMy(VV,model,modelD,options);
    V = minConf_TMP(myObj,V,LB,UB,struct('verbose',0,'maxIter',10));
    model = returnVariationParamsMy(model,V,options);
   %% maximization
    % Optimize DP hyperparameters
    if options.inferType == 6  && im > 3
        modelD = calPre(model);
        modelD = calIntegral(model,options,modelD,3);
        modelD = varsTransform(model,modelD,options);
        model.prior.alpha=(model.U*(model.K-1)+1e-6)...
            /(1e-10+sum(sum(modelD.psiTauAB))-sum(sum(modelD.psiTauB)));
%% UnComment this and you would observe that the BaNPPA-NC gets significantly worse.
% Optimize Gamma distribution hyperparameters  
%         V = extractAmpHyper(model,options);
%         myObj = @(VV)varsGpPpAmpHyper(VV,model,modelD,options);
%         V = minConf_TMP(myObj,V,0.1*ones(size(V)),inf(size(V)),struct('verbose',0,'maxIter',10));
%         model = returnAmpHyper(V,model,options);
    end
    % Optimize GP hyperparameters
    W = extractHyperParamsMy(model,options);
    myObj = @(VV)varsGpPpBoundHyperMy(VV,model,options);
    W  = minConf_TMP(myObj,W,LBHyper,UBHyper,struct('verbose',0,'maxIter',10));
    model = returnHyperParamsMy(model,options,W);
    %% plotting
    [vF,vDiv] = calLikelihood(model,options);
    vEM(im) = sum(vF);
    fprintf('\t\t Hyper Likelihood: %.4f\t Hyper GP %.4f\n',vEM(im),model.GP.logtheta{1}(1));
    % get the query range, sometimes LowB and UpB will be infinity, we should put a constraint on sigma
    [LowB,UpB] = varsGpPpQueryRange(model,options);
    if isinf(LowB) || isinf(UpB) || isnan(vEM(im))%% error break
        nIter = inf;
        try
            load('errorCount.mat');
            ierror = ierror+1;
        catch
            ierror = 1;
        end
        save('errorCount.mat','ierror');
        break;
    end 
    % if likelihood change is small, finish iteration
    if options.figure == 1
        mSave(im,:) = recordModel(model,options);
    end
    if abs(vEM(im)-nEMPre)<=1e-3*abs(nEMPre) && (im > options.minIter  || options.start == 0) 
        break
    end
    nEMPre=vEM(im);
end
if options.figure == 1
    strTitle = {'tau0','tau1','mu','Sigma','theta','alpha','a','b'};
    figure;
    for i=1:8
        subplot(2,4,i);
        plot(mSave(1:im,i),'b-+');
        title(strTitle{i});
    end
end
if ~isinf(nIter)
    nIter = im;
end
margLogL = vEM;
end

function vX = recordModel(model,options)
vX = zeros(1,8);
switch options.inferType
    case 2
        vX(1) = sum(sum(abs(model.amp)));
        vX(3) = sum(cellfun(@(x)sum(abs(x)),model.var.Mu));
        vX(4) = sum(cellfun(@(x)sum(sum(abs(x))),model.var.L));
        vX(5) = sum(cellfun(@(x)sum(abs(x)),model.GP.logtheta));
    case 6
        vX(1) = sum(sum(abs(model.amp.a)));
        vX(2) = sum(sum(abs(model.amp.b)));
        vX(3) = sum(cellfun(@(x)sum(abs(x)),model.var.Mu));
        vX(4) = sum(cellfun(@(x)sum(sum(abs(x))),model.var.L));
        vX(5) = sum(cellfun(@(x)sum(abs(x)),model.GP.logtheta));
        vX(6) = model.prior.alpha;
        vX(7) = model.prior.a;
        vX(8) = model.prior.b;
end
end