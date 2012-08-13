function [rho Psr] = analytical_sol(t, gamma, rho_minus, rho_plus, nx)

%Initial Conditions
P_minus   = rho_minus; P_plus   = rho_plus;
v_minus   = 0        ; v_plus   = 0;
I0 = 1;

%Constants
lambda = P_minus/P_plus;

rho                = zeros(nx+4,1);
rho(1:(nx+4)/2   ) = rho_minus;
rho((nx+4)/2:nx+4) = rho_plus;

Psr                = zeros(nx+4,1);
Psr(1:(nx+4)/2   ) = P_minus;
Psr((nx+4)/2:nx+4) = P_plus;

P    = fzero(@(x) (1-x)^2/(gamma*(1+x)-1+x)-...
    2*gamma/(gamma-1)^2*(1-(x/gamma)^(.5*(1-1/gamma)))^2,1);
Pc   = P_minus/P;
rhoR = rho_plus*(gamma*(1+P)+(P-1))/(gamma*(1+P)-(P-1));
rhoL = rho_plus*(P*lambda^(gamma-1))^(1/gamma);
IR   = I0*P*rho_plus/rhoR;
IL   = I0*P*rho_plus/rhoL;
uc   = 2*sqrt(gamma*I0/(gamma-1))*(1-sqrt(IL/I0));
vs   = uc*rhoR/(rhoR-rho_plus);
va   = uc - sqrt(gamma*(gamma-1)*I0);
vb   = -sqrt(gamma*(gamma-1)*I0);

dista = floor((nx+4)/2+t*va)+1;
distb = floor((nx+4)/2+t*vb)+1;
distc = ceil ((nx+4)/2+t*uc)+1;
dists = ceil ((nx+4)/2+t*vs)+1;

rho(1:distb)     = rho_minus;
rho(distb:dista) = (rho_minus-rhoL)/(distb-dista)*...
    (0:dista-distb)+rho_minus;
rho(dista:distc) = rhoL;
rho(distc:dists) = rhoR;

Psr(1:dista)     = P_minus;
Psr(dists:nx)    = P_plus;
Psr(dista:dists) = Pc;