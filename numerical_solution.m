function c = numerical_solution(tspan,R,D,v,mu,gamma,theta,l,f,n,m,inlet,outlet)
% Computes a numerical solution to the multilayer transport model (1)-(6) using a finite volume
% spatial discretisation and ode15s time stepping.

b0 = inlet{2};
bL = outlet{2};

% Node spacing
x = linspace(0,l(end),n)';

% Initial condition
layer = zeros(n,1);
interface = zeros(m-1,1); cnt = 1;
c0 = zeros(n,1);
for k = 1:n
    for i = m:-1:1
        if x(k) == l(i) && ismember(i,1:m-1)
            interface(cnt) = k;
            cnt = cnt + 1;
            c0(k) = (f(i)+f(i+1))/2;
        elseif x(k) < l(i)
            layer(k) = i;
            c0(k) = f(i);
        end
    end
end
h = x(2)-x(1); % uniform spacing

% Solve ODE system
M = eye(n);
if b0 == 0
    M(1,1) = 0;
end
if bL == 0
    M(n,n) = 0;
end
options = odeset('Mass',M,'MassSingular','yes');
F = @(t,c) Ffunc(t,c,R,D,v,mu,gamma,theta,m,h,layer,interface,n,inlet,outlet);
[~,c] = ode15s(F,[0,tspan],c0,options);

% Remove initial solution
c = c(2:end,:);