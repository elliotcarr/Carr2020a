close all, 
clc, clear all

% save_figs = true;
save_figs = false;

font_size = 26;
if save_figs
    path_name = '../../Paper/Figures/';
    addpath('../export_fig-master')
end

Case = '1';
% Case = '2';
% Case = '3';
% Case = '4';
% Case = '5';
% Case = '6';
% Case = '7';
% Case = '8';
% Case = '9';
% Case = '10';
% Case = '11';
% Case = '12';
% Case = '13';
N = 14; [z,w] = cf(N); % poles and residues appearing in Eq (36)

%% Test Cases
c0 = 1;
% Lstep - toggles step inlet boundary condition on or off
if strcmp(Case,'1') || strcmp(Case,'2') || strcmp(Case,'3') || strcmp(Case,'4')
    m = 2;
    L = 20;
    R = [3,3];
    D = [50,50];
    v = [75,75];
    mu = [2,2];
    gamma = [1,1];
    theta = [0.4,0.4];
    l = [10,L];
    tspan = [1e-3,1e-1,0.6,1.0,2,4];
    f = [0,0];
    % Inlet boundary condition
    if strcmp(Case,'1')
        Lstep = false; t0 = 0;
        a0 = v(1); b0 = D(1); g0 = @(t) c0*v(1); G0 = @(s) v(1)*c0/s;
    elseif strcmp(Case,'2')
        Lstep = true; t0 = 0.5;
        a0 = v(1); b0 = D(1); g0 = @(t) c0*v(1); G0 = @(s) v(1)*c0/s;
    elseif strcmp(Case,'3')
        Lstep = false; t0 = 0;
        a0 = 1; b0 = 0; g0 = @(t) c0; G0 = @(s) c0/s;
    elseif strcmp(Case,'4')
        Lstep = true; t0 = 0.5;
        a0 = 1; b0 = 0; g0 = @(t) c0; G0 = @(s) c0/s;
    end
    inlet = {a0, b0, g0, G0};
    % Outlet boundary condition
    aL = 0; bL = 1; gL = @(t) 0; GL = @(s) 0/s;
    outlet = {aL, bL, gL, GL};
    % Analytical solution
    Neigs = 1000; ftol = 1e-20;
    n = 11;
elseif strcmp(Case,'5') || strcmp(Case,'6') || strcmp(Case,'7') || strcmp(Case,'8')   
    m = 2;
    Lstep = false;
    if strcmp(Case,'8')
        L = 30;
        R = [3,2];
        D = [50,20];
        v = [25,40];
        theta = [0.4,0.25];
        mu = [3,4];
        l = [10,L];
        tspan = [0.2,0.4,0.6,0.8,1e3];
        n = 601;
    else
        L = 30; n = 1*L/2+1;
        mu = [0,0];
        R = [1,1];
        l = [10,L];
        tspan = [0.2,0.4,0.6,0.8];
        if strcmp(Case,'5')
            D = [50,20];
            v = [25,40];
            theta = [0.4,0.25];
        elseif strcmp(Case,'6')
            D = [20,50];
            v = [25,40];
            theta = [0.4,0.25];
        elseif strcmp(Case,'7')
            D = [20,50];
            v = [40,25];
            theta = [0.25,0.4];
        end            
    end
    gamma = [0,0];
    f = [0,0];
    % Inlet boundary condition
    a0 = v(1); b0 = D(1); g0 = @(t) c0*v(1); G0 = @(s) v(1)*c0/s; 
    inlet = {a0, b0, g0, G0};
    % Outlet boundary condition
    aL = 0; bL = 1; gL = @(t) 0; GL = @(s) 0/s;
    outlet = {aL, bL, gL, GL};
elseif strcmp(Case,'9') || strcmp(Case,'10') || strcmp(Case,'11') || strcmp(Case,'13')
    L = 30;
    m = 5;
    l = [10,12,20,22,30]; plot_layers = l;
    D = [7,18,7,18,7];
    theta = [0.4,0.5,0.4,0.5,0.4];
    v = [10,8,10,8,10];
    R = [4.25,14,4.25,14,4.25];
    if strcmp(Case,'13')
        mu = [3,2,3,2,3];
        gamma = [1,2,1,2,1]; 
        f = [0,0,0,1,0];
        tspan = [1e-1,1,2,3.5,4,8];
    else
        mu = [0,0,0,0,0];
        gamma = [0,0,0,0,0];
        f = [0,0,0,0,0];
        tspan = [1e-1,1,2,6,10,14,18];
    end
    % Inlet boundary condition
    if strcmp(Case,'9')
        Lstep = false; 
        a0 = v(1); b0 = D(1); g0 = @(t) v(1)*c0; G0 = @(s) v(1)*c0/s; 
        inlet = {a0, b0, g0, G0};
    elseif strcmp(Case,'10') || strcmp(Case,'13')
        Lstep = true; t0 = 3;
        a0 = v(1); b0 = D(1); g0 = @(t) c0*v(1)*heaviside(t0-t); G0 = @(s) v(1)*c0/s; 
        inlet = {a0, b0, g0, G0};
    elseif strcmp(Case,'11')
        Lstep = false;
        t0 = 3; beta = 1/(t0/2); alpha = t0*beta^2;
        a0 = v(1); b0 = D(1); g0 = @(t) c0*v(1)*alpha*t.*exp(-beta*t); G0 = @(s) c0*v(1)*alpha/(s+beta)^2;     
        inlet = {a0, b0, g0, G0};
    end   
    % Outlet boundary condition
    aL = 0; bL = 1; gL = @(t) 0; GL = @(s) 0/s;
    outlet = {aL, bL, gL, GL};
    n = 601;
elseif strcmp(Case,'12')    
    L = 30;
    m = 5;
    l = [10,12,20,22,30]; plot_layers = l;
    D = [7,18,7,18,7];
    theta = [0.4,0.5,0.4,0.5,0.4];
    v = [10,8,10,8,10];
    mu = [0,0,0,0,0];
    R = [4.25,14,4.25,14,4.25];
    gamma = [0,0,0,0,0];
    f = [0,0,0,0,0];
    % Add artificial layer
    m1 = 14; m2 = 18;
    for i = m:-1:1
        if m1 <= l(i)
            layer1 = i;
        end
    end
    l = [l(1:layer1-1),m1,m2,l(layer1:end)];
    D = [D(1:layer1),D(layer1),D(layer1),D(layer1+1:end)];
    theta = [theta(1:layer1),theta(layer1),theta(layer1),theta(layer1+1:end)];
    v = [v(1:layer1),v(layer1),v(layer1),v(layer1+1:end)];
    mu = [mu(1:layer1),mu(layer1),mu(layer1),mu(layer1+1:end)];
    R = [R(1:layer1),R(layer1),R(layer1),R(layer1+1:end)];
    gamma = [gamma(1:layer1),gamma(layer1),gamma(layer1),gamma(layer1+1:end)];
    f = [f(1:layer1),1,f(layer1),f(layer1+1:end)];
    m = length(l);    
    tspan = [0.1,1,2,6,10,14,18]; 
    % Inlet boundary condition
    Lstep = false;
    a0 = 0; b0 = -1; g0 = @(t) 0; G0 = @(s) 0/s; 
    inlet = {a0, b0, g0, G0};
    % Outlet boundary condition
    aL = 0; bL = 1; gL = @(t) 0; GL = @(s) 0/s;
    outlet = {aL, bL, gL, GL};
    n = 601;
end
x = linspace(0,L,n); % plot solution at these points

%% Semi-analytical solution using Laplace transform
c = zeros(length(x),length(tspan));
for j = 1:length(tspan)
    t = tspan(j);
    for k = 1:length(x)
        Us = @(s) Cfunc(s,x(k),v,D,theta,mu,R,gamma,l,f,inlet,outlet);
        c(k,j) = inverse_laplace_transform(Us,t,N,w,z);
        if Lstep && t > t0
            Us = @(s) Cfunc(s,x(k),v,D,theta,mu,R,zeros(m,1),l,zeros(m,1),inlet,outlet);
            c(k,j) = c(k,j) - inverse_laplace_transform(Us,t-t0,N,w,z);
        end
    end
    fprintf('norm(c,inf) = %1.2e | t = %1.2e\n',norm(c(:,j),inf),t)
end

%% Analytical solution (Case 1-4) or Numerical solution (Cases 5--13)
if strcmp(Case,'1') || strcmp(Case,'2') || strcmp(Case,'3') || strcmp(Case,'4')
    % Eigenvalues
    beta = zeros(Neigs,1);
    if b0 == 0
        g = @(beta) beta.*cot(beta) + v*L/(2*D);
    else
        g = @(beta) beta.*cot(beta) - (beta.^2)*D/(v*L) + v*L/(4*D);
    end
    options = optimoptions('fsolve','Display','none','FunctionTolerance',ftol);
    for i = 1:Neigs
        beta(i) = fsolve(g,(2*i-1)*pi/2,options);
    end
    ca = zeros(length(x),length(tspan));
    for j = 1:length(tspan)
        t = tspan(j);
        for k = 1:length(x)
            ca(k,j) = analytical_solution(x(k),t,R(1),D(1),v(1),mu(1),gamma(1),L,f(1),c0,beta,...
                Neigs,Lstep,t0,inlet);
        end
    end
else
    cn = numerical_solution(tspan,R,D,v,mu,gamma,theta,l,f,n,m,inlet,outlet);
    cn = cn';
end

%% Plot solution
colors = [204,235,197; 
168,221,181
123,204,196
78,179,211
43,140,190
8,104,172
8,64,129]/300;

figure;
vec = get(gcf,'Position');
set(gcf,'Color','w','Position',vec.*[1,1,1.25,1])
if strcmp(Case,'12')
    for i = [1:2,5:m-1]
        plot([l(i),l(i)],[0,1],'k--')
        hold on
    end    
else
    for i = 1:m-1
        plot([l(i),l(i)],[0,1],'k--')
        hold on
    end
end
box on
h = zeros(length(tspan),1);
for j = 1:length(tspan)
    h(j) = plot(x,c(:,j)/c0,'Color',colors(j,:),'Linewidth',2.0);
    if strcmp(Case,'1') || strcmp(Case,'2') || strcmp(Case,'3') || strcmp(Case,'4')
        plot(x,ca(:,j)/c0,'o','Color',colors(j,:),'MarkerSize',10,'MarkerIndices',union([1:20:n],n));
    else
        
        plot(x,cn(:,j)/c0,'o','Color',colors(j,:),'MarkerSize',10,'MarkerIndices',union([1:20:n],n));
    end
end
if strcmp(Case,'8')
    % Steady-state solution for case 8
    lambda(1,1) = (v(1)+sqrt(v(1)^2+4*D(1)*mu(1)))/(2*D(1));
    lambda(1,2) = (v(1)-sqrt(v(1)^2+4*D(1)*mu(1)))/(2*D(1));
    lambda(2,1) = (v(2)+sqrt(v(2)^2+4*D(2)*mu(2)))/(2*D(2));
    lambda(2,2) = (v(2)-sqrt(v(2)^2+4*D(2)*mu(2)))/(2*D(2));
    x1 = x(x<=l(1));
    x2 = x(x>=l(1));
    w = sym('w',[4,1]);
    eq = sym('eq',[4,1]);
    us = sym('us',[2,1]);
    syms x
    us(1) = w(1)*exp(lambda(1,1)*x) + w(2)*exp(lambda(1,2)*x);
    us(2) = w(3)*exp(lambda(2,1)*x) + w(4)*exp(lambda(2,2)*x);
    eq(1) = v(1)*subs(us(1),x,0) - D(1)*subs(diff(us(1),x),x,0) == v(1)*c0;
    eq(2) = subs(diff(us(2),x),x,L) == 0;
    eq(3) = subs(us(1),x,l(1)) == subs(us(2),x,l(1));
    eq(4) = theta(1)*D(1)*subs(diff(us(1),x),x,l(1)) == theta(2)*D(2)*subs(diff(us(2),x),x,l(1));
    wt = solve(eq,w);
    wt = struct2cell(wt);
    plot(x1,eval(wt{1})*exp(lambda(1,1)*x1) + eval(wt{2})*exp(lambda(1,2)*x1),'kx',...
        'MarkerSize',10,'MarkerIndices',union([1:20:length(x1)],length(x1)))
    h(length(tspan)+1) = plot(x2,eval(wt{3})*exp(lambda(2,1)*x2) + eval(wt{4})*exp(lambda(2,2)*x2),...
        'kx','MarkerSize',10,'MarkerIndices',union([1:20:length(x2)],length(x2)));
end
if strcmp(Case,'5') || strcmp(Case,'6') || strcmp(Case,'7') || strcmp(Case,'8')
    axis([0 20 0 1]) % truncate domain for Cases 5-8
else
    axis([0 L 0 1])
end
set(gca,'Fontsize',font_size-2,'TickLabelInterpreter','latex','Xtick',0:5:L,'XMinorTick','on',...
    'YTick',[0,0.5,1],'YMinorTick','on','TickDir','out')
xlabel('$x$ $[\textrm{cm}]$','Interpreter','LaTeX','Fontsize',font_size)
ylabel('$c(x,t)/c_{0}$ $[-]$','Interpreter','LaTeX','Fontsize',font_size)
title(['Case ',Case],'Interpreter','LaTeX','Fontsize',font_size);
if strcmp(Case,'8')
    leg_str = cell(length(tspan)+1,1);
    for j = 1:length(tspan)-1
        leg_str{j} = ['$\,t = ',num2str(tspan(j)),'$'];
    end
    leg_str{j+1} = ['$\,t = 10^{3}$'];
    leg_str{j+2} = ['$\,t\rightarrow\infty$'];
else
    leg_str = cell(length(tspan),1);
    for j = 1:length(tspan)
        leg_str{j} = ['$t = ',num2str(tspan(j)),'$'];
    end
end
leg = legend(h,leg_str,'Interpreter','LaTeX','Fontsize',font_size-6,'Location','EastOutside');
leg.Title.String = '$t$ [days]';
if save_figs
    feval('export_fig',[path_name,['Case',Case]],'-pdf')
end

%% Comparison to previous published results

% Concentration values from Guerrero (2013), Liu et al (1998) and Leij and van Genuchten (1995).
if strcmp(Case,'5')
    % Guerrero (2013)
    CITT(:,1) = [0.884 0.742 0.561 0.375 0.222 0.142 0.063 0.021 0.005 0.001 0.000]';
    CITT(:,2) = [0.963 0.915 0.841 0.746 0.645 0.579 0.480 0.372 0.264 0.168 0.094]';
    CITT(:,3) = [0.987 0.969 0.940 0.901 0.858 0.829 0.781 0.722 0.651 0.567 0.473]';
    CITT(:,4) = [0.995 0.988 0.977 0.962 0.945 0.933 0.914 0.889 0.858 0.819 0.770]';
    % Liu et al (1998)
    GITT(:,1) = [0.884 0.742 0.561 0.374 0.222 0.142 0.063 0.021 0.005 0.001 0.000]';
    GITT(:,2) = [0.963 0.915 0.841 0.746 0.645 0.579 0.480 0.372 0.265 0.169 0.094]';
    GITT(:,3) = [0.987 0.969 0.940 0.901 0.858 0.829 0.781 0.722 0.651 0.567 0.473]';
    GITT(:,4) = [0.995 0.988 0.977 0.962 0.945 0.933 0.914 0.889 0.858 0.819 0.770]';
    % Leij and van Genuchten (1995)
    Linv(:,1) = [0.884 0.742 0.561 0.375 0.222 0.142 0.063 0.021 0.005 0.001 0.000]';
    Linv(:,2) = [0.963 0.915 0.841 0.746 0.645 0.579 0.480 0.372 0.264 0.168 0.094]';
    Linv(:,3) = [0.987 0.969 0.940 0.901 0.858 0.829 0.781 0.722 0.651 0.567 0.473]';
    Linv(:,4) = [0.995 0.988 0.977 0.962 0.945 0.933 0.914 0.889 0.858 0.819 0.770]';
elseif strcmp(Case,'6')
    % Guerrero (2013)
    CITT(:,1) = [0.978 0.868 0.634 0.345 0.131 0.033 0.011 0.003 0.001 0.000 0.000]';
    CITT(:,2) = [0.998 0.984 0.942 0.849 0.693 0.496 0.370 0.257 0.166 0.098 0.054]';
    CITT(:,3) = [1.000 0.998 0.991 0.972 0.930 0.853 0.784 0.699 0.601 0.498 0.395]';
    CITT(:,4) = [1.000 1.000 0.999 0.995 0.986 0.966 0.944 0.913 0.871 0.817 0.751]';
    % Liu et al (1998)
    GITT(:,1) = [0.977 0.867 0.633 0.345 0.131 0.033 0.011 0.003 0.001 0.000 0.000]';
    GITT(:,2) = [0.998 0.984 0.942 0.849 0.693 0.496 0.370 0.257 0.166 0.099 0.054]';
    GITT(:,3) = [1.000 0.998 0.991 0.972 0.929 0.853 0.783 0.698 0.601 0.498 0.395]';
    GITT(:,4) = [1.000 1.000 0.999 0.995 0.986 0.966 0.944 0.913 0.871 0.817 0.750]';
    % Leij and van Genuchten (1995)
    Linv(:,1) = [0.978 0.868 0.634 0.345 0.131 0.033 0.011 0.003 0.001 0.000 0.000]';
    Linv(:,2) = [0.998 0.984 0.942 0.849 0.693 0.496 0.370 0.257 0.166 0.098 0.054]';
    Linv(:,3) = [1.000 0.998 0.991 0.972 0.930 0.853 0.784 0.699 0.601 0.498 0.395]';
    Linv(:,4) = [1.000 1.000 0.999 0.995 0.986 0.966 0.944 0.913 0.871 0.817 0.751]';
elseif strcmp(Case,'7')
    % Guerrero (2013)
    CITT(:,1) = [0.999 0.988 0.928 0.764 0.496 0.152 0.049 0.013 0.003 0.000 0.000]';
    CITT(:,2) = [1.000 1.000 0.999 0.995 0.976 0.780 0.600 0.418 0.262 0.148 0.075]';
    CITT(:,3) = [1.000 1.000 1.000 1.000 0.998 0.940 0.870 0.773 0.653 0.522 0.393]';
    CITT(:,4) = [1.000 1.000 1.000 1.000 0.999 0.979 0.952 0.911 0.851 0.774 0.681]';    
    % Liu et al (1998)
    GITT(:,1) = [0.999 0.987 0.928 0.763 0.495 0.152 0.050 0.013 0.003 0.000 0.000]';
    GITT(:,2) = [1.000 1.000 0.999 0.995 0.976 0.779 0.600 0.418 0.262 0.148 0.075]';
    GITT(:,3) = [1.000 1.000 1.000 1.000 0.998 0.939 0.870 0.773 0.653 0.522 0.393]';    
    GITT(:,4) = [1.000 1.000 1.000 1.000 0.999 0.978 0.952 0.910 0.851 0.774 0.681]';    
    % Leij and van Genuchten (1995)
    Linv(:,1) = [0.999 0.988 0.928 0.764 0.496 0.152 0.049 0.013 0.003 0.000 0.000]';
    Linv(:,2) = [1.000 1.000 0.999 0.995 0.976 0.780 0.600 0.417 0.262 0.148 0.075]';
    Linv(:,3) = [1.000 1.000 1.000 1.000 0.998 0.940 0.870 0.773 0.653 0.522 0.393]';
    Linv(:,4) = [1.000 1.000 1.000 1.000 0.999 0.979 0.952 0.911 0.851 0.774 0.681]';
end

%% Tables 4 and 5
if strcmp(Case,'1') || strcmp(Case,'2') || strcmp(Case,'3') || strcmp(Case,'4')
    fprintf('%s',Case)
    for j = 1:length(tspan)
        fprintf(' & \\num{%1.2e}',norm(c(:,j)-ca(:,j),2))
    end
    fprintf('\\\\\n');
elseif strcmp(Case,'5') || strcmp(Case,'6') || strcmp(Case,'7')
    fprintf('%s',Case)
    xvec = zeros(11,1);
    for k = 1:11
        xvec(k) = find(x == 2*(k-1));
    end
    for k = 1:11
        fprintf(' & %g',x(xvec(k)))
        for j = 1:length(tspan)
            vec = [round(c(xvec(k),j),3),CITT(k,j),GITT(k,j),Linv(k,j)];
            ind = find(vec(1:4) ~= vec(1));
            fprintf(' & %1.3f',vec(1))
            for i = 2:4
                if ismember(i,ind)
                    fprintf(' & \\cellcolor{tableshade}{%1.3f}',vec(i))
                else
                    fprintf(' & %1.3f',vec(i))
                end
            end
        end
        fprintf('\\\\\n');
    end
end