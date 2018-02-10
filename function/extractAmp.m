function V = extractAmp(model,options,user)
% EXTRACT AMPLITUDE PARAMETER FOR SINGLE USER
switch options.inferType
    case 6
        V = [model.amp.a(:,user);model.amp.b(:,user)];
    case 2
        V = model.amp(:,user);
end
end