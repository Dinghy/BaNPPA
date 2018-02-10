function model = returnAmpHyper(V,model,options)
model.prior.a = V(1);
model.prior.b = V(2);
end