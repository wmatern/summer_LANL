function flux = ACMy(u,a,Cy)

u_half_p1 = u(3) - u(2)+ eps;
u_half_m1 = u( - u_m + eps;

% alpha (j,i) = 2*abs((abs(u_half).^2.5-abs(u_half_m).^2.5)./...
%     (abs(u_half).^2.5+abs(u_half_m)^2.5+eps));

alpha = 66*(abs(u_p-2*u+u_m)./...
    (abs(u_p-u)+abs(u-u_m)+eps)).^2;

L_half   = 0.5*(abs(a)-Cy*a.^2).*...
            (u_half-minmod(u_half_m,minmod(u_half,...
            u_half(j+1,i))));
    
L_half_m = 0.5*(abs(a)-Cy*a.^2).*...
            (u_half-minmod(u_half_m,minmod(u_half,...
            u_half(j+1,i))));

L     (j,i) = sign(L_half).*(max(0,max(sign(L_half).*...
    minmod(alpha.*L_half(j-1,i),L_half),sign(L_half(j-1,i)).*...
    minmod(L_half(j-1,i), alpha.*L_half))));

if u_half == 0
    mu2   (j,i) = 0;
else
    mu2   (j,i) = (L(j+1,i)-L(j,i))./u_half;
end


g1    (j,i) = 0.5*(abs(a+mu2(j,i))-Cy*(a+mu2(j,i)).^2).*u_half;
g     (j,i) = minmod(g1(j,i),g1(j-1,i));

if u_half == 0
    gamma(j,i) = 0;
else
    gamma(j,i) = (g(j+1,i)-g(j,i))./u_half;
end

flux = 0.5*(a.*(u+u_p)+(L(j,i)+L(j+1,i))...
    + (g(j,i)+g(j+1,i))-abs(a+mu2(j,i)+gamma(j,i)).*(u_half));

end