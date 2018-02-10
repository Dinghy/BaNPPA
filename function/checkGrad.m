function norm = checkGrad(objFunc,vInput,varargin) 
%  function to check whether gradient caculation is correct
% 
nEps=1e-3;


norm = 0;
[f,df]=feval(objFunc,vInput,varargin{:});

for i=1:length(vInput)      
  vInput(i) = vInput(i) - 0.5*nEps;
  E0 = feval(objFunc,vInput,varargin{:});
  vInput(i) = vInput(i) + nEps;
  E1 =feval(objFunc,vInput,varargin{:});
  vInput(i) = vInput(i) - 0.5*nEps;
  
  diff_grad = (E1 - E0)/nEps;
  
  norm = norm + abs(diff_grad - df(i));
  
  % fprintf('%d: %f vs %f: Diff %f\n',i, diff_grad, df(i),abs(diff_grad - df(i)));
end
% fprintf('Likelihood: %f\n',f);
end 