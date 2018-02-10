function [vData] = sampleIPP(flambda,nBeg,nEnd,nRlim,nRatio)
% Inputs:   # nTlim  : limit on time 0~nTlim
%           # nRlim  : limit on intensity function
%           # flambda: function for lambda
%           # nratio : ratio for the function
% Outputs:  # vData  : return data

nPos = nBeg;
nDebug = 0;
nR=1;

if nargin==4
    nR=nRatio;
end
nRlim=nRlim*nR;

vData = [];
while nPos<=nEnd
    nPos=nPos+exprnd(1/nRlim);
    
    if nDebug
        disp(nPos);
    end
    
    if unifrnd(0,nRlim)<=nR*feval(flambda,nPos)
        vData=[vData,nPos];
    end
end

end

