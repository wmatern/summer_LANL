clear all; clc;% close all;
% Lax-Wendroff Scheme for solving the compressible Euler Equations

% Domain initialization
nx   =  64 ; ny   =  64 ;
xmax =  2.0; ymax =  2.0;
xmin = -2.0; ymin = -2.0;
dx = (xmax-xmin)/(nx-1);
dy = (ymax-ymin)/(ny-1);
[X,Y] = meshgrid(xmin:dx:xmax,ymin:dy:ymax);

% Time initialization
tmax = 100;
CFL  = 0.05;
dt   = CFL*min(dx,dy);
Cx   = dt/dx;
Cy   = dt/dy;

% Switches
TVD   = 2;                  % 1=MinMod, 2=Superbee, 3=Nonlinear Superbee
BC    = 2;                  % 1=Reflective, 2=Periodic
SCH   = 3;                  % 1=Lax-Wendroff, 2=MPDATA, 3=ACM
alpha = 2;                  % Constatnt for ACM

% Physical Constants
gamma = 1.4;

% Initial Conditions
RHO = ones (ny+4,nx+4);
MX  = zeros(ny+4,nx+4);
MY  = ones (ny+4,nx+4);
E   = ones (ny+4,nx+4);
P   = ones (ny+4,nx+4);

% RHO(3:ny+2,3:nx+2) = exp(-2*(X.^2+Y.^2)) + 1;
% RHO(ny/2-ny/4:ny/2+ny/4,:) = 2;
RHO(ny/2+2:ny+4,:) = 0.5;
% MY (ny/2+2:ny+4,:) = 0.5;
MY = RHO;
P  (ny/2+2:ny+4,:) = 0.5;

RHO_New1 = RHO;
RHO_New  = RHO;
RHO_New3 = RHO;

figure(1),surf(X,Y,RHO(3:ny+2,3:nx+2)), title('RHO')
figure(2),surf(X,Y,MY (3:ny+2,3:nx+2)), title('MY' )
figure(3),plot(MY(:,34),'.-'         ), title('MY' )

RHO_New = RHO;


for t=0:dt:tmax
MY = RHO;
% Boundary Conditions
if BC == 1 % Reflective BC
    RHO(2,:) = RHO(3,:); P(2,:) = P(3,:); E(2,:) = E(3,:);
    RHO(1,:) = RHO(2,:); P(1,:) = P(2,:); E(1,:) = E(2,:);
    RHO(:,2) = RHO(:,3); P(:,2) = P(:,3); E(:,2) = E(:,3);
    RHO(:,1) = RHO(:,2); P(:,1) = P(:,2); E(:,1) = E(:,2);
    
    MX(2,:) = -MX(3,:); MY(2,:) = -MY(3,:);
    MX(1,:) =  MX(2,:); MY(1,:) =  MY(2,:);
    MX(:,2) = -MX(:,3); MY(:,2) = -MY(:,3);
    MX(:,1) =  MX(:,2); MY(:,1) =  MY(:,2);
    
    RHO(ny+3,:) = RHO(ny+2,:); P(ny+3,:) = P(ny+2,:); E(ny+3,:) = E(ny+2,:);
    RHO(ny+4,:) = RHO(ny+2,:); P(ny+4,:) = P(ny+2,:); E(ny+4,:) = E(ny+2,:);
    RHO(:,nx+3) = RHO(:,nx+2); P(:,nx+3) = P(:,nx+2); E(:,nx+3) = E(:,nx+2);
    RHO(:,nx+4) = RHO(:,nx+2); P(:,nx+4) = P(:,nx+2); E(:,nx+4) = E(:,nx+2);
    
    MX(ny+3,:) = -MX(ny+2,:); MY(ny+3,:) = -MY(ny+2,:);
    MX(ny+4,:) =  MX(ny+3,:); MY(ny+4,:) =  MY(ny+3,:);
    MX(:,nx+3) = -MX(:,nx+2); MY(:,nx+3) = -MY(:,nx+2);
    MX(:,nx+4) =  MX(:,nx+3); MY(:,nx+4) =  MY(:,nx+3);
    
elseif BC == 2 % Periodic BC
    RHO(2,:) = RHO(ny+2,:); P(2,:) = P(ny+2,:); E(2,:) = E(ny+2,:);
    RHO(1,:) = RHO(ny+1,:); P(1,:) = P(ny+1,:); E(1,:) = E(ny+1,:);
    RHO(:,2) = RHO(:,nx+2); P(:,2) = P(:,nx+2); E(:,2) = E(:,nx+2);
    RHO(:,1) = RHO(:,nx+1); P(:,1) = P(:,nx+1); E(:,1) = E(:,nx+1);
    
    MX(2,:) = MX(ny+2,:); MY(2,:) = MY(ny+2,:);
    MX(1,:) = MX(ny+1,:); MY(1,:) = MY(ny+1,:);
    MX(:,2) = MX(:,nx+2); MY(:,2) = MY(:,nx+2);
    MX(:,1) = MX(:,nx+1); MY(:,1) = MY(:,nx+1);
    
    RHO(ny+3,:) = RHO(3,:); P(ny+3,:) = P(3,:); E(ny+3,:) = E(3,:);
    RHO(ny+4,:) = RHO(4,:); P(ny+4,:) = P(4,:); E(ny+4,:) = E(4,:);
    RHO(:,nx+3) = RHO(:,3); P(:,nx+3) = P(:,3); E(:,nx+3) = E(:,3);
    RHO(:,nx+4) = RHO(:,4); P(:,nx+4) = P(:,4); E(:,nx+4) = E(:,4);
    
    MX(ny+3,:) = MX(3,:); MY(ny+3,:) = MY(3,:);
    MX(ny+4,:) = MX(4,:); MY(ny+4,:) = MY(4,:);
    MX(:,nx+3) = MX(:,3); MY(:,nx+3) = MY(:,3);
    MX(:,nx+4) = MX(:,4); MY(:,nx+4) = MY(:,4);
end

    RHO_New(2,:) = RHO_New(ny+2,:);
    RHO_New(1,:) = RHO_New(ny+1,:);
    RHO_New(:,2) = RHO_New(:,nx+2);
    RHO_New(:,1) = RHO_New(:,nx+1);
    RHO_New(ny+3,:) = RHO_New(3,:);
    RHO_New(ny+4,:) = RHO_New(4,:);
    RHO_New(:,nx+3) = RHO_New(:,3);
    RHO_New(:,nx+4) = RHO_New(:,4);
    
    RHO_New1(2,:) = RHO_New1(ny+2,:);
    RHO_New1(1,:) = RHO_New1(ny+1,:);
    RHO_New1(:,2) = RHO_New1(:,nx+2);
    RHO_New1(:,1) = RHO_New1(:,nx+1);
    RHO_New1(ny+3,:) = RHO_New1(3,:);
    RHO_New1(ny+4,:) = RHO_New1(4,:);
    RHO_New1(:,nx+3) = RHO_New1(:,3);
    RHO_New1(:,nx+4) = RHO_New1(:,4);
    
    RHO_New3(2,:) = RHO_New3(ny+2,:);
    RHO_New3(1,:) = RHO_New3(ny+1,:);
    RHO_New3(:,2) = RHO_New3(:,nx+2);
    RHO_New3(:,1) = RHO_New3(:,nx+1);
    RHO_New3(ny+3,:) = RHO_New3(3,:);
    RHO_New3(ny+4,:) = RHO_New3(4,:);
    RHO_New3(:,nx+3) = RHO_New3(:,3);
    RHO_New3(:,nx+4) = RHO_New3(:,4);



% First Pass
j = 1:ny+2;
i = 1:nx+2;
% X Direction
% Density
RHOx (j,i) = 0.5*(RHO(j+2,i+2)+RHO(j+2,i+1))...
    - Cx/2*(MX(j+2,i+2)-MX(j+2,i+1));

% X Momentum
MXx  (j,i) = 0.5*(MX(j+2,i+2)+MX(j+2,i+1))...
    - Cx/2*((MX(j+2,i+2).^2./RHO(j+2,i+2)+P(j+2,i+2))...
          - (MX(j+2,i+1).^2./RHO(j+2,i+1)+P(j+2,i+1)));
% Y Momentum
MYx  (j,i) = 0.5*(MY(j+2,i+2)+MY(j+2,i+1))...
    - Cx/2*(MX(j+2,i+2).*MY(j+2,i+2)./RHO(j+2,i+2)...
          - MX(j+2,i+1).*MY(j+2,i+1)./RHO(j+2,i+1));
% Energy Flux
Ex   (j,i) = 0.5*(E(j+2,i+2)+E(j+2,i+1))...
    - Cx/2*((E(j+2,i+2)+P(j+2,i+2)).*MX(j+2,i+2)./RHO(j+2,i+2)...
          - (E(j+2,i+1)+P(j+2,i+1)).*MX(j+2,i+1)./RHO(j+2,i+1));
% Pressure
Px   (j,i) = (gamma-1)*Ex(j,i)-0.5*((MXx(j,i).^2+MYx(j,i).^2)./RHOx(j,i));

% Y Direction
RHOy (j,i) = 0.5*(RHO(j+2,i+2)+RHO(j+1,i+2))...
    - Cy/2*(MY(j+2,i+2)-MY(j+1,i+2));
% X Momentum
MXy  (j,i) = 0.5*(MX(j+2,i+2)+MX(j+1,i+2))...
    - Cy/2*(MY(j+2,i+2).*MX(j+2,i+2)./RHO(j+2,i+2)...
          - MY(j+1,i+2).*MX(j+1,i+2)./RHO(j+1,i+2));
% Y Momentum
MYy  (j,i) = 0.5*(MY(j+2,i+2)+MY(j+1,i+2))...
    - Cy/2*((MY(j+2,i+2).^2./RHO(j+2,i+2)+P(j+2,i+2))...
          - (MY(j+1,i+2).^2./RHO(j+1,i+2)+P(j+1,i+2)));
% Energy Flux
Ey   (j,i) = 0.5*(E(j+2,i+2)+E(j+1,i+2))...
    - Cy/2*((E(j+2,i+2)+P(j+2,i+2)).*MY(j+2,i+2)./RHO(j+2,i+2)...
          - (E(j+1,i+2)+P(j+1,i+2)).*MY(j+1,i+2)./RHO(j+1,i+2));
% Pressure
Py   (j,i) = (gamma-1)*Ey(j,i)-0.5*((MXy(j,i).^2+MYy(j,i).^2)./RHOy(j,i));

% Second Pass
j = 3:ny+2;
i = 3:nx+2;

nuplusy  = abs(MYy(j-1,i-2)./RHOy(j-1,i-2)...
    +sqrt(gamma*Py(j-1,i-2)./RHOy(j-1,i-2)))*Cy;
nuplusx  = abs(MXx(j-2,i-1)./RHOx(j-2,i-1)...
    +sqrt(gamma*Px(j-2,i-1)./RHOx(j-2,i-1)))*Cx;
numinusy = abs(MYy(j-2,i-2)./RHOy(j-2,i-2)...
    +sqrt(gamma*Py(j-2,i-2)./RHOy(j-2,i-2)))*Cy;
numinusx = abs(MXx(j-2,i-2)./RHOx(j-2,i-2)...
    +sqrt(gamma*Px(j-2,i-2)./RHOx(j-2,i-2)))*Cx;

% Mass Conservation
% if SCH == 1 % Use Lax-Wendroff
RHO = RHO_New1;       
RHO_New1 (j,i) = RHO(j,i) - Cx*(MXx(j-2,i-1)-MXx(j-2,i-2))...
                         - Cy*(MYy(j-1,i-2)-MYy(j-2,i-2))...
    - wminusx(RHO,j,i,numinusx,TVD).*(RHO(j  ,i  )-RHO(j  ,i-1))...
    + wplusx (RHO,j,i,nuplusx ,TVD).*(RHO(j  ,i+1)-RHO(j  ,i  ))...
    - wminusy(RHO,j,i,numinusy,TVD).*(RHO(j  ,i  )-RHO(j-1,i  ))...
    + wplusy (RHO,j,i,nuplusy ,TVD).*(RHO(j+1,i  )-RHO(j  ,i  ));

% elseif SCH == 2 % Use MPDATA
RHO = RHO_New;
U_halfp (j,i) = (MX(j,i+1)./RHO(j,i+1) + MX(j,i  )./RHO(j,i  ))*Cx/2;
U_halfm (j,i) = (MX(j,i  )./RHO(j,i  ) + MX(j,i-1)./RHO(j,i-1))*Cx/2;
V_halfp (j,i) = (MY(j+1,i)./RHO(j+1,i) + MY(j,i  )./RHO(j,i  ))*Cy/2;
V_halfm (j,i) = (MY(j,i  )./RHO(j,i  ) + MY(j-1,i)./RHO(j-1,i))*Cy/2;

RHO_New (j,i) = RHO(j,i) - (MPDATA(RHO(j,i  ),RHO(j,i+1),U_halfp(j,i))...
                          - MPDATA(RHO(j,i-1),RHO(j,i  ),U_halfm(j,i)))...
                         - (MPDATA(RHO(j,i  ),RHO(j+1,i),V_halfp(j,i))...
                          - MPDATA(RHO(j-1,i),RHO(j,i  ),V_halfm(j,i)));
                      
U_d_p   (j,i) = dt/dx*((abs(U_halfp(j,i))-U_halfp(j,i).^2).*...
    d_px(j,i  ,RHO_New) - U_halfp(j,i).*V_halfp(j,i).*d_py(j,i  ,RHO_New));
U_d_m   (j,i) = dt/dx*((abs(U_halfm(j,i))-U_halfm(j,i).^2).*...
    d_px(j,i-1,RHO_New) - U_halfm(j,i).*V_halfm(j,i).*d_py(j,i-1,RHO_New));
V_d_p   (j,i) = dt/dy*((abs(V_halfp(j,i))-V_halfp(j,i).^2).*...
    d_py(j  ,i,RHO_New) - U_halfp(j,i).*V_halfp(j,i).*d_px(j  ,i,RHO_New));
V_d_m   (j,i) = dt/dy*((abs(V_halfm(j,i))-V_halfm(j,i).^2).*...
    d_py(j-1,i,RHO_New) - U_halfm(j,i).*V_halfm(j,i).*d_px(j-1,i,RHO_New));

RHO_New (j,i) = RHO_New(j,i) - (MPDATA(RHO_New(j,i  ),RHO_New(j,i+1),U_d_p(j,i))...
                              - MPDATA(RHO_New(j,i-1),RHO_New(j,i  ),U_d_m(j,i)))...
                             - (MPDATA(RHO_New(j,i  ),RHO_New(j+1,i),V_d_p(j,i))...
                              - MPDATA(RHO_New(j-1,i),RHO_New(j,i  ),V_d_m(j,i)));

% elseif SCH == 3 % Use ACM
RHO = RHO_New3;
nuplusy  = abs(MYy(j-1,i-2)./RHOy(j-1,i-2)...
    -sqrt(gamma*Py(j-1,i-2)./RHOy(j-1,i-2)))*Cy;
nuplusx  = abs(MXx(j-2,i-1)./RHOx(j-2,i-1)...
    -sqrt(gamma*Px(j-2,i-1)./RHOx(j-2,i-1)))*Cx;
numinusy = abs(MYy(j-2,i-2)./RHOy(j-2,i-2)...
    -sqrt(gamma*Py(j-2,i-2)./RHOy(j-2,i-2)))*Cy;
numinusx = abs(MXx(j-2,i-2)./RHOx(j-2,i-2)...
    -sqrt(gamma*Px(j-2,i-2)./RHOx(j-2,i-2)))*Cx;
RHO_New3 (j,i) = RHO(j,i) - Cx*(MXx(j-2,i-1)-MXx(j-2,i-2))...
                         - Cy*(MYy(j-1,i-2)-MYy(j-2,i-2))...
    - ACMx(RHO,j  ,i-1,alpha).*(RHO(j  ,i  )-RHO(j  ,i-1))...
    + ACMx(RHO,j  ,i  ,alpha).*(RHO(j  ,i+1)-RHO(j  ,i  ))...
    - ACMy(RHO,j-1,i  ,alpha).*(RHO(j  ,i  )-RHO(j-1,i  ))...
    + ACMy(RHO,j  ,i  ,alpha).*(RHO(j+1,i  )-RHO(j  ,i  ))...
    - wminusx(RHO,j,i,numinusx,TVD).*(RHO(j  ,i  )-RHO(j  ,i-1))...
    + wplusx (RHO,j,i,nuplusx ,TVD).*(RHO(j  ,i+1)-RHO(j  ,i  ))...
    - wminusy(RHO,j,i,numinusy,TVD).*(RHO(j  ,i  )-RHO(j-1,i  ))...
    + wplusy (RHO,j,i,nuplusy ,TVD).*(RHO(j+1,i  )-RHO(j  ,i  ));


% end

RHO = RHO_New;
MY  = RHO_New;
plot1(1:length(j)  ,1) = RHO_New1(j,34);
if BC == 2, plot1(length(j)+j-2,1) = RHO_New1(j,34); end
plot2(1:length(j)  ,1) = RHO_New(j,34);
if BC == 2, plot2(length(j)+j-2,1) = RHO_New(j,34); end
plot3(1:length(j)  ,1) = RHO_New3(j,34);
if BC == 2, plot3(length(j)+j-2,1) = RHO_New3(j,34); end

if mod(t,10*dt) == 0
    %     figure(1),surf(X,Y,RHO(j,i)), title('RHO')
    %     figure(2),surf(X,Y,MY (j,i)), title('MY' )
    figure(3),plot(plot1,'.-'), title('RHO')
    axis([1 128 0 1.1])
    figure(4),plot(plot2,'.-'), title('RHO')
    axis([1 128 0 1.1])
    figure(5),plot(plot3,'.-'), title('RHO')
    axis([1 128 0 1.1])
end
end