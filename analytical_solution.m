function c = analytical_solution(x,t,R,D,v,mu,gamma,L,f,C0,beta,N,Lstep,ts,Lbnd)
% Analytical solutions for single-layer (homogeneous) media [see van Genuchten and Alves (1982)]

b0 = Lbnd{2};
u = v*sqrt(1+(4*mu*D/(v^2)));

if (v*L/D) > min(5 + 40*v*t/(R*L),100) || (v*L/D) > min(5 + 40*v*(t-ts)/(R*L),100)
    
    % Use approximate solutions
    if b0 == 0
        A = exp(-mu*t/R)*(1 - 0.5*erfc((R*x-v*t)/(2*sqrt(D*R*t))) - ...
            0.5*exp(v*x/D)*erfc((R*x+v*t)/(2*sqrt(D*R*t))) - ...
            0.5*(2 + (v*(2*L-x))/D + v^2*t/(D*R))*exp(v*L/D)*erfc((R*(2*L-x)+v*t)/(2*sqrt(D*R*t))) + ...
            sqrt(v^2*t/(pi*D*R))*exp(v*L/D - (R/(4*D*t))*(2*L-x + v*t/R)^2));
        B3 = 0.5*exp((v-u)*x/(2*D))*erfc((R*x-u*t)/(2*sqrt(D*R*t))) + ...
            0.5*exp((v+u)*x/(2*D))*erfc((R*x+u*t)/(2*sqrt(D*R*t))) + ...
            (u-v)/(2*(u+v))*exp(((v+u)*x-2*u*L)/(2*D))*erfc((R*(2*L-x)-u*t)/(2*sqrt(D*R*t))) + ...
            (u+v)/(2*(u-v))*exp(((v-u)*x+2*u*L)/(2*D))*erfc((R*(2*L-x)+u*t)/(2*sqrt(D*R*t))) - ...
            (v^2/(2*mu*D))*exp(v*L/D - mu*t/R)*erfc((R*(2*L-x)+v*t)/(2*sqrt(D*R*t)));
        B4 = 1 + (u-v)/(u+v)*exp(-u*L/D);
        B = B3/B4;
        c = gamma/mu + (f-gamma/mu)*A + (C0 - gamma/mu)*B;
        if Lstep && t > ts
            t = t - ts;
            Bs3 = 0.5*exp((v-u)*x/(2*D))*erfc((R*x-u*t)/(2*sqrt(D*R*t))) + ...
                0.5*exp((v+u)*x/(2*D))*erfc((R*x+u*t)/(2*sqrt(D*R*t))) + ...
                (u-v)/(2*(u+v))*exp(((v+u)*x-2*u*L)/(2*D))*erfc((R*(2*L-x)-u*t)/(2*sqrt(D*R*t))) + ...
                (u+v)/(2*(u-v))*exp(((v-u)*x+2*u*L)/(2*D))*erfc((R*(2*L-x)+u*t)/(2*sqrt(D*R*t))) - ...
                (v^2/(2*mu*D))*exp(v*L/D - mu*t/R)*erfc((R*(2*L-x)+v*t)/(2*sqrt(D*R*t)));
            Bs = Bs3/B4;
            c = c - C0*Bs;
        end
        
    else
        
        A = exp(-mu*t/R)*(1 - 0.5*erfc((R*x-v*t)/(2*sqrt(D*R*t))) - ...
            sqrt((v^2*t)/(pi*D*R))*exp(-((R*x-v*t)^2)/(4*D*R*t)) + ...
            0.5*(1 + v*x/D + v^2*t/(D*R))*exp(v*x/D)*erfc((R*x+v*t)/(2*sqrt(D*R*t))) - ...
            sqrt(4*v^2*t/(pi*D*R))*(1 + v/(4*D)*(2*L-x+(v*t)/R))*exp(v*L/D - (R/(4*D*t))*(2*L-x+(v*t)/R)^2) + ...
            v/D*(2*L-x + 3*v*t/(2*R) + (v/(4*D))*(2*L-x+(v*t)/R)^2)*exp(v*L/D)*erfc((R*(2*L-x)+v*t)/(2*sqrt(D*R*t))));
        B3 = v/(v+u)*exp((v-u)*x/(2*D))*erfc((R*x-u*t)/(2*sqrt(D*R*t))) + ...
            v/(v-u)*exp((v+u)*x/(2*D))*erfc((R*x+u*t)/(2*sqrt(D*R*t))) + ...
            v^2/(2*mu*D)*exp(v*x/D - mu*t/R)*erfc((R*x+v*t)/(2*sqrt(D*R*t))) + ...
            v^2/(2*mu*D)*(v*(2*L-x)/D + v^2*t/(D*R) + 3 + v^2/(mu*D))*exp(v*L/D - mu*t/R)*erfc((R*(2*L-x)+v*t)/(2*sqrt(D*R*t))) - ...
            v^3/(mu*D)*sqrt(t/(pi*D*R))*exp(v*L/D - mu*t/R - R/(4*D*t)*(2*L-x+v*t/R)^2) + ...
            v*(u-v)/((u+v)^2)*exp(((v+u)*x-2*u*L)/(2*D))*erfc((R*(2*L-x)-u*t)/(2*sqrt(D*R*t))) - ...
            v*(u+v)/((u-v)^2)*exp(((v-u)*x+2*u*L)/(2*D))*erfc((R*(2*L-x)+u*t)/(2*sqrt(D*R*t)));
        B4 = 1 - (u-v)^2/(u+v)^2*exp(-u*L/D);
        B = B3/B4;
        c = gamma/mu + (f-gamma/mu)*A + (C0 - gamma/mu)*B;
        
        if Lstep && t > ts
            t = t - ts;
            Bs3 = v/(v+u)*exp((v-u)*x/(2*D))*erfc((R*x-u*t)/(2*sqrt(D*R*t))) + ...
                v/(v-u)*exp((v+u)*x/(2*D))*erfc((R*x+u*t)/(2*sqrt(D*R*t))) + ...
                v^2/(2*mu*D)*exp(v*x/D - mu*t/R)*erfc((R*x+v*t)/(2*sqrt(D*R*t))) + ...
                v^2/(2*mu*D)*(v*(2*L-x)/D + v^2*t/(D*R) + 3 + v^2/(mu*D))*exp(v*L/D - mu*t/R)*erfc((R*(2*L-x)+v*t)/(2*sqrt(D*R*t))) - ...
                v^3/(mu*D)*sqrt(t/(pi*D*R))*exp(v*L/D - mu*t/R - R/(4*D*t)*(2*L-x+v*t/R)^2) + ...
                v*(u-v)/((u+v)^2)*exp(((v+u)*x-2*u*L)/(2*D))*erfc((R*(2*L-x)-u*t)/(2*sqrt(D*R*t))) - ...
                v*(u+v)/((u-v)^2)*exp(((v-u)*x+2*u*L)/(2*D))*erfc((R*(2*L-x)+u*t)/(2*sqrt(D*R*t)));
            Bs = Bs3/B4;
            c = c - C0*Bs;
        end
        
    end
    
else
    
    % Use full expansion solutions    
    A = 0.0;
    B2 = 0.0;
    for m = 1:N
        if b0 == 0
            E = 2*beta(m)*sin(beta(m)*x/L)/(beta(m)^2 + (v*L/(2*D))^2 + v*L/(2*D));
        else
            E = (2*v*L/D * beta(m)*(beta(m)*cos(beta(m)*x/L) + v*L/(2*D)*sin(beta(m)*x/L))) / ...
                ((beta(m)^2 + (v*L/(2*D))^2 + v*L/D)*(beta(m)^2 + (v*L/(2*D))^2));
        end
        A = A + E * exp(v*x/(2*D) - mu*t/R - v^2*t/(4*D*R) - beta(m)^2*D*t/(L^2 * R));
        B2 = B2 + (E * (beta(m)^2 + (v*L/(2*D))^2) * exp(v*x/(2*D) - mu*t/R - v^2*t/(4*D*R) - beta(m)^2*D*t/(L^2*R))) / ...
            (beta(m)^2 + (v*L/(2*D))^2 + mu*L^2/D);
    end
    
    if b0 == 0
        B1 = (exp((v-u)*x/(2*D)) + (u-v)/(u+v)*exp((v+u)*x/(2*D)-u*L/D)) / ...
            (1 + (u-v)/(u+v)*exp(-u*L/D));
    else
        B1 = (exp((v-u)*x/(2*D)) + (u-v)/(u+v)*exp(((v+u)*x-2*u*L)/(2*D))) / ...
            ((u+v)/(2*v) - (u-v)^2/(2*v*(u+v))*exp(-u*L/D));
    end
    B = B1-B2;
    c = gamma/mu + (f-gamma/mu)*A + (C0 - gamma/mu)*B;
    
    if Lstep && t > ts
        B2 = 0.0; t = t - ts;
        for m = 1:N
            if b0 == 0
                E = 2*beta(m)*sin(beta(m)*x/L)/(beta(m)^2 + (v*L/(2*D))^2 + v*L/(2*D));
            else
                E = (2*v*L/D * beta(m)*(beta(m)*cos(beta(m)*x/L) + v*L/(2*D)*sin(beta(m)*x/L))) / ...
                    ((beta(m)^2 + (v*L/(2*D))^2 + v*L/D)*(beta(m)^2 + (v*L/(2*D))^2));
            end
            B2 = B2 + (E * (beta(m)^2 + (v*L/(2*D))^2) * exp(v*x/(2*D) - mu*t/R - v^2*t/(4*D*R) - beta(m)^2*D*t/(L^2*R))) / ...
                (beta(m)^2 + (v*L/(2*D))^2 + mu*L^2/D);
        end
        B = B1-B2;
        c = c - C0*B;
    end
    
end