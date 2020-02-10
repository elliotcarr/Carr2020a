function P = Pfunc(i,s,x,mu,R,gamma,f,l,beta,lambda,Psi,m,inlet,outlet)
% Computes P_{i}(x,s) defined in Eqs (31)-(33) and Table 1.

a0 = inlet{1}; aL = outlet{1};

P = (gamma(i)/s + R(i)*f(i))/(mu(i) + R(i)*s);
if i == 1
    P = (1 + a0/beta(1)*(lambda(1,1)*Psi(1,2,x) - lambda(1,2)*Psi(1,2,l(1))*Psi(1,1,x)))*P;
elseif i == m
    P = (1 + aL/(beta(m))*(lambda(m,2)*Psi(m,1,x) - lambda(m,1)*Psi(m,1,l(m-1))*Psi(m,2,x)))*P;
end