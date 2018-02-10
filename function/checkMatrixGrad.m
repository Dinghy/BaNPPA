function norm = checkMatrixGrad(objFunc,vInput,model) 
%  function to check whether gradient caculation is correct
% 
nEps=1e-5;
norm = 0;
[~,df]=feval(objFunc,vInput,model);

for i=1:length(vInput)     
  vInput(i) = vInput(i) - 0.5*nEps;
  E0 = feval(objFunc,vInput,model);
  vInput(i) = vInput(i) + nEps;
  E1 =feval(objFunc,vInput,model);
  vInput(i) = vInput(i) - 0.5*nEps;
  diff_grad = (E1 - E0)/nEps;
  nPlus = sum(sum(abs(diff_grad - df{i})));
  norm = norm + nPlus;
  
  fprintf('%d: Diff %f\n',i,nPlus);
end
end 