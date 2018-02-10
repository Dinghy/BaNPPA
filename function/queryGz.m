function [mGG,mDGG] = queryGz(z,gz,dgz)
% QUERYGZ Summary of this function goes here
% linear interpolation
debug = 0;
if debug
    mGG = -z;
    mDGG = ones(size(z));
else
    z = log10(z);
    mPos = 49999*(max(-11,min(11,z))+11)/22+1;
    mPos_l = floor(mPos);
    mPos_u = ceil(mPos);
    mRatio_l = (mPos_u-mPos)./(mPos_u-mPos_l);
    mRatio_u = (mPos-mPos_l)./(mPos_u-mPos_l);
    mBequal = mPos_l==mPos_u;
    mGG = (mBequal).*gz(mPos_l)+(~mBequal).*...
        (mRatio_l.*gz(mPos_l)+mRatio_u.*gz(mPos_u));
    mGG(z<-11) = 0;
    if nargout == 2
        mDGG = (mBequal).*dgz(mPos_l)+(~mBequal).*...
        (mRatio_l.*dgz(mPos_l)+mRatio_u.*dgz(mPos_u));
        mDGG(z<-11) = 2;
    end
    
end
end