function [mPhi,DmPhidx] = MatrixPhi(V, model)
% MATRIXPHI Summary of this function goes here
%   Detailed explanation goes here
% 1 dimension first
D = size(model.Xm,2); % dimension
M = size(model.Xm,1); % pseudo inputs
% initialization
mPhi = ones(M);
if nargout==2
    DmPhidx = cell(1,length(model.GP.logtheta{1}));
    for d=1:length(model.GP.logtheta{1})
        DmPhidx{d} = zeros(M);
    end
end

mPhi_save = cell(1,D);
mTmp_a = cell(1,D);
mTmp_b = cell(1,D);
valpha_sqrt = exp(V(1:end-2));
valpha = valpha_sqrt.^2;
nGamma = exp(2*V(end-1));
for d=1:D
    nAlphaSqrt = valpha_sqrt(d);
    nAlpha = valpha(d);
    vZ_d = model.Xm(:,d);
    mZbar = repmat(vZ_d,1,M);%+repmat(vZ_d',M,1);
    mQ=mZbar-mZbar';
    mZbar = mZbar+mZbar';
    mZbar = mZbar/(2*nAlphaSqrt);
    mSqQ = mQ.^2/4/nAlpha;
    mExpQ = exp(-mSqQ);
    mZbarMax=mZbar-model.Tend(d)/nAlphaSqrt;
    mZbarMin=mZbar-model.Tstart(d)/nAlphaSqrt;
    mPhi_save{d} = sqrt(pi)*nAlphaSqrt/2*mExpQ.*(-erf(mZbarMax)+erf(mZbarMin));
    mPhi=mPhi.*mPhi_save{d};
    if nargout==2
        mTmp_a{d}=(mSqQ+0.5)*2;
        mTmp_b{d}=nAlphaSqrt*mExpQ.*...
            (mZbarMax.*exp(-(mZbarMax).^2)-mZbarMin.*exp(-(mZbarMin).^2));
    end
end

if nargout==2
    if D==2
        DmPhidx{1}=nGamma^2*(mPhi.*mTmp_a{1}+mTmp_b{1}.*mPhi_save{2});
        DmPhidx{2}=nGamma^2*(mPhi.*mTmp_a{2}+mTmp_b{2}.*mPhi_save{1});
    else
        DmPhidx{1}=nGamma^2*(mPhi_save{1}.*mTmp_a{d}+mTmp_b{d});
    end
end
mPhi=nGamma^2*mPhi;
end

