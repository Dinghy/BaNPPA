function model = returnAmp(V,model,options,user)
%% RETURN THE AMPLITUDE 
switch options.inferType
    case 2
        model.amp(:,user) = V;
    case 6
        model.amp.a(:,user) = V(1:length(V)/2);
        model.amp.b(:,user) = V(length(V)/2+1:end);
end
end