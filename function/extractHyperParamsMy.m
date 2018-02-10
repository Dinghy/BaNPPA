function W = extractHyperParamsMy(model,options)
% EXTRACT PARAMETERS FOR KERNEL TO OPTIMIZE
% 0 do not optimize noise 
% 1 optimize noise  
% 2 optimize single noise
W=[];
nNoise = 0;
for c=1:model.K
    if ~iscolumn(model.GP.logtheta{c})
        W=[W;model.GP.logtheta{c}(1:end-2)'];
    else
        W=[W;model.GP.logtheta{c}(1:end-2)];
    end
    nNoise = nNoise+model.GP.logtheta{c}(end);
end
if model.derivative == 1
    W = [W;nNoise];
else
    W = [W;nNoise/model.K];
end

end