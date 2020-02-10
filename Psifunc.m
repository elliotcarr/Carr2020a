function y = Psifunc(i,j,x,l,lambda)
% Computes Psi_{i,j}(x,s) for j = 1,2 defined in Table 1.

if j == 1
    y = exp(lambda(i,1)*(x-l(i)));
else
    if i == 1
        y = exp(lambda(i,2)*x);
    else
        y = exp(lambda(i,2)*(x-l(i-1)));
    end
end