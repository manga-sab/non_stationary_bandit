function F = randsample_array(P,X) 
p = cumsum([0; P(1:end-1).'; 1+1e3*eps]); 
[a a] = histc(rand,p); 
F = X(a);