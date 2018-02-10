function model = returnHyperParamsMy(model,options,V)
% RETURN PARAMETERS FOR KERNEL TO OPTIMIZE
% 0 do not optimize noise 
% 1 optimize noise  
% 2 optimize single noise, all element share the same noise

nStart=1;
switch options.noise
    case 0
        N=length(model.GP.logtheta{1})-2;
        for c=1:model.K
            model.GP.logtheta{c}(1:end-2)=V(nStart:nStart+N-1);
            nStart=nStart+N;
        end
    case 1
        N=length(model.GP.logtheta{1})-1;
        for c=1:model.K
            model.GP.logtheta{c}(1:end-2)=V(nStart:nStart+N-2);
            model.GP.logtheta{c}(end)=V(nStart+N-1);
            nStart=nStart+N;
        end
    case 2
        N=length(model.GP.logtheta{1})-2;
        for c=1:model.K
            model.GP.logtheta{c}(1:end-2) = V(nStart:nStart+N-1);
            model.GP.logtheta{c}(end) = V(end);
            nStart=nStart+N;
        end
end
end

