function model = varsgpCreateMy(data,model,options)
% VARSGPCONCREATE Summary of this function goes here
%   Detailed explanation goes here
model.const.C = 0.5772156649;
model.U = data.U;                  % user number
model.K = options.cluster;         % cluster number
model.S = ceil(model.K/2);         % cluster for each user
model.m = options.m;               % pseudo points numbers
model.Tstart = data.XBeg;          % time start
model.Tend = data.XEnd;            % time end
model.NumData = zeros(model.U,1);  % number of points in each dataset
%% store data
for u=1:model.U
    model.data{u}.feature=unique(data.train{u}.X,'rows');
    model.NumData(u) = size(model.data{u}.feature,1);
end

% prior at each feature point
model.cov = {'covSum',{'covSEard','covNoise'}};

%% amplitude
% 3 parameters to be inferred
switch options.inferType
    case 2
        model.A = mean(model.NumData);
        model.amp = 1+rand(model.K,model.U);
    case 6
        %% we can maximize at the beginning version for a,b in the Gamma prior
        model.A = mean(model.NumData);                  % what integral term you want
        [a,b] = estimateGamma(model.NumData,model.A,0);
        model.prior.a = a;
        model.prior.b = b;
        fprintf('Infer: (Shape a)%.4f\t, (Rate b)%.4f\n',model.prior.a,model.prior.b);
        model.prior.alpha = 1.1;
        if isfield(options,'alpha')
            model.prior.alpha = options.alpha;
        end
        model.amp.s = (model.NumData+model.prior.a-1)/(model.prior.b+model.A);
        model.amp.a = zeros(model.K-1,model.U);
        model.amp.b = zeros(model.K-1,model.U);
        model = varsGpPpTauAmp(model,options);
end
%% pseudo inputs positions
model.Xm = linspace(model.Tstart(1),model.Tend(1),options.m)';
nStep = 0.1;
model.plot.Xm = [model.Tstart(1):nStep:model.Tend(1)]';
model.plot.row = length(model.plot.Xm);

model.m = size(model.Xm,1);
model.prior.gbase = model.A/prod(model.Tend-model.Tstart);
model.prior.g = sqrt(model.prior.gbase)*ones(model.m,1);
%% hyper parameter for Kmm
D = length(model.Tend); % 1 since it is time-sequence
% if data.dataSource == 1 || data.dataSource == 5 
% dd = log(4.3081);                                 % For three additional Synthetic data sets only;
% else                                              % others use log((model.Tend-model.Tstart)')/2;
dd = log((model.Tend-model.Tstart)')/2;
% end

dd(dd==-Inf)=0;
model.GP.logtheta=cell(1,model.K);
model.var.Mu=cell(1,model.K);
model.var.L=cell(1,model.K);
model.var.Sigma=cell(1,model.K);
for c=1:model.K
    %% hyper-parameter in GP
    model.GP.logtheta{c}=zeros(D+2,1);
    model.GP.logtheta{c}(1:D) = dd;
    model.GP.logtheta{c}(D+1) = 0.5*log(model.prior.gbase);
    model.GP.logtheta{c}(D+2) = log(1e-4);
    %% prior for q(m,S)
    model.var.Mu{c} = model.prior.g;
    model.var.Sigma{c} = feval(model.cov{:},model.GP.logtheta{c},model.Xm);
    model.var.L{c} = chol(model.var.Sigma{c},'lower');
end
end
