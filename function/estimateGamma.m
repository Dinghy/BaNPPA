function [a,b] = estimateGamma(arrCount,T,boutput)
%  input 
%      arrCount -- counting matrix
%      T        -- T window
%      boutput  -- 1 output;0 not output
% output
%      a        -- shape parameter
%      b        -- rate parameter
%% adding a prior
if ~iscolumn(arrCount)
    arrCount = arrCount';
end
gamDis = fitdist(arrCount,'Gamma');
a = gamDis.a;   % shape
b = 1/gamDis.b; % rate
i = 1;
try
    while i<10
        if boutput == 1
            fprintf('\tGamma a:%.4f\t b:%.4f\n',gamDis.a,gamDis.b);
        end
        i = i+1;
        % estimation
        arrRate = (arrCount+a-1)/(b+T);
        % maximization
        gamDis = fitdist(arrRate,'Gamma');
        a = gamDis.a;   % shape
        b = 1/gamDis.b; % rate
    end
catch
    error('Not suitable!');
end
end