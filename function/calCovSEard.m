function sRes = calCovSEard(hyp, x, z)
%  partly from GPML
if nargin<2, K = '(D+1)'; return; end              % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode

[n,D] = size(x);
ell = exp(hyp(1:D));                               % characteristic length scale
sf2 = exp(2*hyp(D+1));                                         % signal variance
if D==1 % one dimensional case
    % precompute squared distances
    if dg                                                               % vector kxx
        K = zeros(size(x,1),1);
    else
        if xeqz                                                 % symmetric matrix Kxx
            K = sq_dist(diag(1./ell)*x');
        else                                                   % cross covariances Kxz
            K = sq_dist(diag(1./ell)*x',diag(1./ell)*z');
        end
    end
    
    sRes.Knm = sf2*exp(-K/2);     % covariance
    sRes.DknmDx = cell(1,D+2);
    for i=1:D+2
        if i<=D                                              % length scale parameters
            if dg
                sRes.DknmDx{i} = sRes.Knm*0;
            else
                sRes.DknmDx{i} = sRes.Knm.*K;
            end
        else
            sRes.DknmDx{i} = zeros(size(K));
        end
    end
else
    % precompute squared distances
    if dg                                                               % vector kxx
        K = zeros(size(x,1),1);
    else
        if xeqz                                                 % symmetric matrix Kxx
            K = sq_dist(diag(1./ell)*x');
        else                                                   % cross covariances Kxz
            K = sq_dist(diag(1./ell)*x',diag(1./ell)*z');
        end
    end
    
    sRes.Knm = sf2*exp(-K/2);     % covariance
    sRes.DknmDx = cell(1,D+2);
    for i=1:D+2
        if i<=D                                              % length scale parameters
            if dg
                sRes.DknmDx{i} = sRes.Knm*0;
            else
                if xeqz
                    sRes.DknmDx{i} = sRes.Knm.*sq_dist(x(:,i)'/ell(i));
                else
                    sRes.DknmDx{i} = sRes.Knm.*sq_dist(x(:,i)'/ell(i),z(:,i)'/ell(i));
                end
            end
        else
            sRes.DknmDx{i} = zeros(size(K));
        end
    end
end
