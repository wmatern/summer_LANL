function flux = ACMy(u,a,Cy,q)

u_half_p3 = u(7) - u(6); + eps;
u_half_p2 = u(6) - u(5); + eps;
u_half_p1 = u(5) - u(4); + eps;
u_half_m1 = u(4) - u(3); + eps;
u_half_m2 = u(3) - u(2); + eps;
u_half_m3 = u(2) - u(1); + eps;

alpha = 2*abs((abs(u_half_p1)^2.5-abs(u_half_m1)^2.5)/...
    (abs(u_half_p1)^2.5+abs(u_half_m1)^2.5+eps));

% alpha = q*(abs(u(5)-2*u(4)+u(3))./...
%     (abs(u(5)-u(4))+abs(u(4)-u(3))+eps)).^2;

L_half_p2 = 0.5*(abs(a(2))-Cy*a(2).^2)*...
            (u_half_p2-minmod(u_half_p1,minmod(u_half_p2,...
            u_half_p3)));
        
L_half_p1 = 0.5*(abs(a(2))-Cy*a(2).^2)*...
            (u_half_p1-minmod(u_half_m1,minmod(u_half_p1,...
            u_half_p2)));
    
L_half_m1 = 0.5*(abs(a(2))-Cy*a(2).^2)*...
            (u_half_m1-minmod(u_half_m2,minmod(u_half_m1,...
            u_half_p1)));

L_half_m2 = 0.5*(abs(a(2))-Cy*a(2).^2)*...
            (u_half_m2-minmod(u_half_m3,minmod(u_half_m2,...
            u_half_m1)));

L_p      = sign(L_half_p2)*(max(0,       max(sign(L_half_p2)*...
    minmod(alpha*L_half_p1,       L_half_p2),sign(L_half_p2)*...
    minmod(      L_half_p1, alpha*L_half_p2))));

L        = sign(L_half_p1)*(max(0,       max(sign(L_half_p1)*...
    minmod(alpha*L_half_m1,       L_half_p1),sign(L_half_p1)*...
    minmod(      L_half_m1, alpha*L_half_p1))));

L_m      = sign(L_half_m1)*(max(0,       max(sign(L_half_m1)*...
    minmod(alpha*L_half_m2,       L_half_m1),sign(L_half_m1)*...
    minmod(      L_half_m2, alpha*L_half_m1))));

if u_half_p1 == 0
    mu_p = 0;
else
    mu_p = (L_p-L)/u_half_p1;
end
if u_half_m1 == 0
    mu_m = 0;
else
    mu_m = (L-L_m)/u_half_m1;
end

g_half_p = 0.5*(abs(a(2)+mu_p )-Cy*(a(2)+mu_p )^2)*u_half_p1;
g_half_m = 0.5*(abs(a(2)+mu_m )-Cy*(a(2)+mu_m )^2)*u_half_m1;
g_p1     = minmod(g_half_p,g_half_m);

if u_half_p2 == 0
    mu_p2 = 0;
else
    mu_p2 = (L_p-L)/u_half_p2;
end


g_half_p = 0.5*(abs(a(2)+mu_p2)-Cy*(a(2)+mu_p2)^2)*u_half_p2;
g_half_m = 0.5*(abs(a(2)+mu_p )-Cy*(a(2)+mu_p )^2)*u_half_p1;
g_p2     = minmod(g_half_p,g_half_m);

if u_half_p1 == 0
    gamma_p1 = 0;
else
    gamma_p1 = (g_p2-g_p1)./u_half_p1;
end

flux = 0.5*(a(2).*(u(4)+u(5))+(L+L_p)...
    + (g_p1+g_p2)-abs(a(2)+mu_p+gamma_p1)*(u_half_p1));

end