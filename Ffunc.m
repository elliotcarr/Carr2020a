function F = Ffunc(t,c,R,D,v,mu,gamma,theta,m,h,layer,interface,n,inlet,outlet)
% Computes F(c) defined in Appendix A.

a0 = inlet{1};
b0 = inlet{2};
g0 = inlet{3};
aL = outlet{1};
bL = outlet{2};
gL = outlet{3};

F = zeros(n,1);

% Treatment of inlet boundary condition
if b0 == 0
    F(1) = a0*c(1) - g0(t);
else
    F(1) = D(1)*(c(2)-c(1))/h - v(1)*(c(1)+c(2))/2 + D(1)/b0*g0(t) + ...
        (v(1)-D(1)*a0/b0)*c(1) + h/2*(-mu(1)*c(1)+gamma(1));
    F(1) = F(1)/(h/2*R(1));
end

% Treatment of outlet boundary condition
if bL == 0
    F(n) = aL*c(n) - gL(t);
else
    F(n) = D(m)/bL*gL(t) - (v(m) + D(m)*aL/bL)*c(n) - D(m)*(c(n)-c(n-1))/h + ...
        v(m)*(c(n-1)+c(n))/2 + h/2*(-mu(m)*c(n)+gamma(m));
    F(n) = F(n)/(h/2*R(m));
end


for k = 2:n-1
    
    [mem,indx] = ismember(k,interface);
    
    if mem
        
        i = indx; % node located at interface between layers i and i+1
        F(k) = theta(i+1)*D(i+1)*(c(k+1)-c(k))/h - theta(i+1)*v(i+1)*(c(k)+c(k+1))/2 - ...
            theta(i)*D(i)*(c(k)-c(k-1))/h + v(i)*theta(i)*(c(k-1)+c(k))/2 + ...
            h/2*(theta(i)*(-mu(i)*c(k)+gamma(i))+theta(i+1)*(-mu(i+1)*c(k)+gamma(i+1)));
        F(k) = F(k)/(h/2*(theta(i)*R(i)+theta(i+1)*R(i+1)));
        
    else
        
        i = layer(k); % node located in interior of layer i
        F(k) = D(i)*(c(k+1)-c(k))/h - v(i)*(c(k)+c(k+1))/2 - D(i)*(c(k)-c(k-1))/h ...
            + v(i)*(c(k-1)+c(k))/2 + h*(-mu(i)*c(k)+gamma(i));
        F(k) = F(k)/(h*R(i));
        
    end
    
end