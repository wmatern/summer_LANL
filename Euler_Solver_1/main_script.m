clear all; clc; %close all;
% Lax-Wendroff Scheme for solving the compressible Euler Equations

% Domain initialization
nx   =  32 ; ny   =  32 ;
xmax =  2.0; ymax =  2.0;
xmin = -2.0; ymin = -2.0;
dx = (xmax-xmin)/(nx-1);
dy = (ymax-ymin)/(ny-1);
[X,Y] = meshgrid(xmin:dx:xmax,ymin:dy:ymax);

% Time initialization
tmax = 100;
CFL  = 0.2;
dt   = CFL*min(dx,dy);
Cx   = dt/dx;
Cy   = dt/dy;

% Switches
TVD = 2;                  % 1=MinMod, 2=Superbee, 3=Nonlinear Superbee
BC  = 2;                  % 1=Reflective, 2=Periodic
SCH = 4;                  % 1=Lax-Wendroff, 2=MPDATA,4=FEM

% Physical Constants
gamma = 1.4;

% Initial Conditions
RHO = ones (ny+4,nx+4);
MX  = zeros(ny+4,nx+4);
MY  = ones(ny+4,nx+4);
E   = ones (ny+4,nx+4)*20;
P   = ones (ny+4,nx+4);

% RHO(3:ny+2,3:nx+2) = exp(-2*(X.^2+Y.^2)) + 1;
% RHO(ny/2-ny/4:ny/2+ny/4,:) = 2;
RHO(ny/4:3*ny/4,:) = .5;
%RHO(3:nx+2,3:ny+2) = 1+.25*cos((2*pi/4)*Y);
%     RHO(2,:) = RHO(3,:); 
%     RHO(1,:) = RHO(2,:); 
%     RHO(:,2) = RHO(:,3); 
%     RHO(:,1) = RHO(:,2);
%     RHO(ny+3,:) = RHO(ny+2,:);
%     RHO(ny+4,:) = RHO(ny+2,:);
%     RHO(:,nx+3) = RHO(:,nx+2);
%     RHO(:,nx+4) = RHO(:,nx+2);
MY = RHO;
E = RHO;
P = (gamma-1)*(E-0.5*((MX.^2 + MY.^2)./RHO));
% MY (ny/2+2:ny+4,:) = 0.5;
% P  (ny/2+2:ny+4,:) = 0.5;

figure(1),surf(X,Y,RHO(3:ny+2,3:nx+2)), title('RHO')
figure(2),surf(X,Y,MY(3:ny+2,3:nx+2) ), title('MY' )
%figure(3),plot(MY(:,34),'.-'         ), title('MY' )

RHO_New = RHO;
err = eps*(randi([-1,1],nx,ny));
RHO(3:ny+2,3:nx+2) = RHO(3:ny+2,3:nx+2) + 0*eps*(randi([-1,1],nx,ny));


for t=0:dt:tmax
    t1 = tic;
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
Px   (j,i) = (gamma-1)*(Ex(j,i)-0.5*((MXx(j,i).^2+MYx(j,i).^2)./RHOx(j,i)));

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
Py   (j,i) = (gamma-1)*(Ey(j,i)-0.5*((MXy(j,i).^2+MYy(j,i).^2)./RHOy(j,i)));

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
if SCH == 1 % Use Lax-Wendroff
% nuplusy  = abs(MYy(j-1,i-2)./RHOy(j-1,i-2)...
%     - sqrt(gamma*Py(j-1,i-2)./RHOy(j-1,i-2)))*Cy;
% nuplusx  = abs(MXx(j-2,i-1)./RHOx(j-2,i-1)...
%     - sqrt(gamma*Px(j-2,i-1)./RHOx(j-2,i-1)))*Cx;
% numinusy = abs(MYy(j-2,i-2)./RHOy(j-2,i-2)...
%     - sqrt(gamma*Py(j-2,i-2)./RHOy(j-2,i-2)))*Cy;
% numinusx = abs(MXx(j-2,i-2)./RHOx(j-2,i-2)...
%     - sqrt(gamma*Px(j-2,i-2)./RHOx(j-2,i-2)))*Cx;
       
RHO_New (j,i) = RHO(j,i) - Cx*(MXx(j-2,i-1)-MXx(j-2,i-2))...
                         - Cy*(MYy(j-1,i-2)-MYy(j-2,i-2))...
    - wminusx(RHO,j,i,numinusx,TVD).*(RHO(j  ,i  )-RHO(j  ,i-1))...
    + wplusx (RHO,j,i,nuplusx ,TVD).*(RHO(j  ,i+1)-RHO(j  ,i  ))...
    - wminusy(RHO,j,i,numinusy,TVD).*(RHO(j  ,i  )-RHO(j-1,i  ))...
    + wplusy (RHO,j,i,nuplusy ,TVD).*(RHO(j+1,i  )-RHO(j  ,i  ));

elseif SCH == 2 % Use MPDATA
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

elseif BC == 2 && SCH == 4
    %Remove extra points used for BCs
    RHO_New = RHO(3:ny+2,3:nx+2);
    MX_t = MX(3:ny+2,3:nx+2);
    MY_t = MY(3:ny+2,3:nx+2);
    a = 1;
    U = reshape(MX_t./RHO_New,nx*ny,1);
    V = reshape(MY_t./RHO_New,nx*ny,1);
    X = reshape(X,nx*ny,1);
    Y = reshape(Y,nx*ny,1);
    %tri = delaunay(X,Y);
    n = 1;
    for k = 1:nx*ny
        if rem(k,ny) ~= 0 && k < (nx-1)*ny
            tri(n,1:3) = [k, k+1, k+ny];
            n = n+1;
            tri(n,1:3) = [k+1, k+1+ny, k+ny];
            n = n+1;
        end
    end
    
    alias = [(1:ny)',(nx*ny-(ny-1):nx*ny)';(1+ny:ny:nx*ny-(ny-1))',(2*ny:ny:nx*ny)'];
    alias = fliplr(alias);
    nalias = 1:nx*ny; %any node that is not on the right side of alias
    nalias(alias(:,2)) = [];
    
    %Assemble the Global stiffness matrix from local matrices
    M = zeros(nx*ny,nx*ny);
    K = zeros(nx*ny,nx*ny);
    n = 1:length(tri);
    Me1 = zeros(3,3);
    Me2 = zeros(3,3);
    Me3 = zeros(3,3);
    h(n) = .707*(xmax-xmin)/(nx-1);
    A(n) = abs(.5*(X(tri(n,1)).*Y(tri(n,2)) - X(tri(n,1)).*Y(tri(n,3))...
                 + X(tri(n,2)).*Y(tri(n,3)) - X(tri(n,2)).*Y(tri(n,1))...
                 + X(tri(n,3)).*Y(tri(n,1)) - X(tri(n,3)).*Y(tri(n,2))))';
    tic
    for n = 1:length(tri)
        if rem(n,2) == 1
            u = -(U(tri(n,3))+U(tri(n,2))+U(tri(n,1))+U(tri(n+1,2)))/4;
            v = -(V(tri(n,3))+V(tri(n,2))+V(tri(n,1))+V(tri(n+1,2)))/4;
            x = (X(tri(n,3))+X(tri(n,2))+X(tri(n,1))+X(tri(n+1,2)))/4;
            y = (Y(tri(n,3))+Y(tri(n,2))+Y(tri(n,1))+Y(tri(n+1,2)))/4;
        else
            u = -(U(tri(n,3))+U(tri(n,2))+U(tri(n,1))+U(tri(n-1,1)))/4;
            v = -(V(tri(n,3))+V(tri(n,2))+V(tri(n,1))+V(tri(n-1,1)))/4;
            x = (X(tri(n,3))+X(tri(n,2))+X(tri(n,1))+X(tri(n-1,1)))/4;
            y = (Y(tri(n,3))+Y(tri(n,2))+Y(tri(n,1))+Y(tri(n-1,1)))/4;
        end
        t_u(n) = u;
        t_v(n) = v;
        t_x(n) = x;
        t_y(n) = y;
        %v = -1;
        %u = 0;
        uv_mag = sqrt(u^2+v^2);
        Me1(1:3,1:3) = A(n)*[1/6,1/12,1/12;1/12,1/6,1/12;1/12,1/12,1/6];
        Me2(1:3,1:3) = (a*h(n)*u/(12*uv_mag))*...
            [Y(tri(n,2))-Y(tri(n,3)),Y(tri(n,2))-Y(tri(n,3)),Y(tri(n,2))-Y(tri(n,3));... 
             Y(tri(n,3))-Y(tri(n,1)),Y(tri(n,3))-Y(tri(n,1)),Y(tri(n,3))-Y(tri(n,1));...
             Y(tri(n,1))-Y(tri(n,2)),Y(tri(n,1))-Y(tri(n,2)),Y(tri(n,1))-Y(tri(n,2))];
        Me3(1:3,1:3) = (a*h(n)*v/(12*uv_mag))*...
            [X(tri(n,3))-X(tri(n,2)),X(tri(n,3))-X(tri(n,2)),X(tri(n,3))-X(tri(n,2));... 
             X(tri(n,1))-X(tri(n,3)),X(tri(n,1))-X(tri(n,3)),X(tri(n,1))-X(tri(n,3));...
             X(tri(n,2))-X(tri(n,1)),X(tri(n,2))-X(tri(n,1)),X(tri(n,2))-X(tri(n,1))];
        M(tri(n,1),tri(n,:)) = M(tri(n,1),tri(n,:)) + Me1(1,:) + Me2(1,:) + Me3(1,:);
        M(tri(n,2),tri(n,:)) = M(tri(n,2),tri(n,:)) + Me1(2,:) + Me2(2,:) + Me3(2,:);
        M(tri(n,3),tri(n,:)) = M(tri(n,3),tri(n,:)) + Me1(3,:) + Me2(3,:) + Me3(3,:);
        Ke2r1(1,1:3) = ((u/6) + (u*a*h(n)/(8*A(n))) * ((u/uv_mag)*(Y(tri(n,2))-Y(tri(n,3))) + (v/uv_mag)*(X(tri(n,3))-X(tri(n,2))))) * [Y(tri(n,2))-Y(tri(n,3)), Y(tri(n,3))-Y(tri(n,1)), Y(tri(n,1))-Y(tri(n,2))]... 
                     + ((v/6) + (v*a*h(n)/(8*A(n))) * ((u/uv_mag)*(Y(tri(n,2))-Y(tri(n,3))) + (v/uv_mag)*(X(tri(n,3))-X(tri(n,2))))) * [X(tri(n,3))-X(tri(n,2)), X(tri(n,1))-X(tri(n,3)), X(tri(n,2))-X(tri(n,1))];
        Ke2r2(1,1:3) = ((u/6) + (u*a*h(n)/(8*A(n))) * ((u/uv_mag)*(Y(tri(n,3))-Y(tri(n,1))) + (v/uv_mag)*(X(tri(n,1))-X(tri(n,3))))) * [Y(tri(n,2))-Y(tri(n,3)), Y(tri(n,3))-Y(tri(n,1)), Y(tri(n,1))-Y(tri(n,2))]... 
                     + ((v/6) + (v*a*h(n)/(8*A(n))) * ((u/uv_mag)*(Y(tri(n,3))-Y(tri(n,1))) + (v/uv_mag)*(X(tri(n,1))-X(tri(n,3))))) * [X(tri(n,3))-X(tri(n,2)), X(tri(n,1))-X(tri(n,3)), X(tri(n,2))-X(tri(n,1))];
        Ke2r3(1,1:3) = ((u/6) + (u*a*h(n)/(8*A(n))) * ((u/uv_mag)*(Y(tri(n,1))-Y(tri(n,2))) + (v/uv_mag)*(X(tri(n,2))-X(tri(n,1))))) * [Y(tri(n,2))-Y(tri(n,3)), Y(tri(n,3))-Y(tri(n,1)), Y(tri(n,1))-Y(tri(n,2))]... 
                     + ((v/6) + (v*a*h(n)/(8*A(n))) * ((u/uv_mag)*(Y(tri(n,1))-Y(tri(n,2))) + (v/uv_mag)*(X(tri(n,2))-X(tri(n,1))))) * [X(tri(n,3))-X(tri(n,2)), X(tri(n,1))-X(tri(n,3)), X(tri(n,2))-X(tri(n,1))];
        K(tri(n,1),tri(n,:)) = K(tri(n,1),tri(n,:)) + Ke2r1(1,1:3);
        K(tri(n,2),tri(n,:)) = K(tri(n,2),tri(n,:)) + Ke2r2(1,1:3);
        K(tri(n,3),tri(n,:)) = K(tri(n,3),tri(n,:)) + Ke2r3(1,1:3);
    end
    toc
    figure(6);
    plot3(t_x,t_y,t_v,'.');
    
    %Enforce periodic boundary conditions.
    for m = 1:length(alias)
        K(:,alias(m,1)) = K(:,alias(m,1)) + K(:,alias(m,2));
        M(:,alias(m,1)) = M(:,alias(m,1)) + M(:,alias(m,2));
    end
    for m = 1:length(alias)
        K(alias(m,1),:) = K(alias(m,1),:) + K(alias(m,2),:);
        M(alias(m,1),:) = M(alias(m,1),:) + M(alias(m,2),:);
    end
    
    K(alias(:,2),:) = [];
    K(:,alias(:,2)) = [];
    M(alias(:,2),:) = [];
    M(:,alias(:,2)) = [];
    RHO_New(alias(:,2)) = [];
    RHO_New = RHO_New';
    
    M = sparse(M);
    K = sparse(K);
    
    options = odeset('Mass',M);
    f = @(t, RHO_New)(-K * RHO_New);
    
    tic
    [T,RHO_New] = ode45(f,[0 dt],RHO_New,options);
    toc
    RHO_New = RHO_New(end,1:end)';
    %theta = 0; %1 - forward Euler
    %A = (M+(1-theta)*dt*K);
    %RHO_New = (M+(1-theta)*dt*K) \ ((M - theta*dt*K)*RHO_New);
    
    %Restore the normal form of the variables
    temp = zeros((ny)*(nx),1);
    temp(nalias) = RHO_New;
    
    temp(alias(:,2)) = temp(alias(:,1)); %Having trouble with repeats
    temp(alias(:,2)) = temp(alias(:,1));
    
    RHO_New = reshape(temp,ny,nx);
    n = 3:ny+2; m = 3:nx+2;
    temp = zeros(ny+4,nx+4);
    temp(n,m) = RHO_New(n-2,m-2);
    RHO_New = temp;
    X = reshape(X,ny,nx);
    Y = reshape(Y,ny,nx);
end

% X Momentum
% nuplusy  = abs(MYy(j-1,i-2)./RHOy(j-1,i-2))*Cy;
% nuplusx  = abs(MXx(j-2,i-1)./RHOx(j-2,i-1))*Cx;
% numinusy = abs(MYy(j-2,i-2)./RHOy(j-2,i-2))*Cy;
% numinusx = abs(MXx(j-2,i-2)./RHOx(j-2,i-2))*Cx;

MX_New  (j,i) = MX(j,i) - Cx*((MXx(j-2,i-1).^2./RHOx(j-2,i-1)+Px(j-2,i-1))...
                            - (MXx(j-2,i-2).^2./RHOx(j-2,i-2)+Px(j-2,i-2)))...
    - Cy*(MYy(j-1,i-2).*MXy(j-1,i-2)./RHOy(j-1,i-2)...
        - MYy(j-2,i-2).*MXy(j-2,i-2)./RHOy(j-2,i-2))...
    - wminusx(MX,j,i,numinusx,TVD).*(MX(j  ,i  )-MX(j  ,i-1))...
    + wplusx (MX,j,i,nuplusx ,TVD).*(MX(j  ,i+1)-MX(j  ,i  ))...
    - wminusy(MX,j,i,numinusy,TVD).*(MX(j  ,i  )-MX(j-1,i  ))...
    + wplusy (MX,j,i,nuplusy ,TVD).*(MX(j+1,i  )-MX(j  ,i  ));

% Y Momentum
MY_New  (j,i) = MY(j,i) - Cy*((MYy(j-1,i-2).^2./RHOy(j-1,i-2)+Py(j-1,i-2))...
                            - (MYy(j-2,i-2).^2./RHOy(j-2,i-2)+Py(j-2,i-2)))...
    - Cx*(MXx(j-2,i-1).*MYx(j-2,i-1)./RHOx(j-2,i-1)...
        - MXx(j-2,i-2).*MYx(j-2,i-2)./RHOx(j-2,i-2))...
    - wminusx(MY,j,i,numinusx,TVD).*(MY(j  ,i  )-MY(j  ,i-1))...
    + wplusx (MY,j,i,nuplusx ,TVD).*(MY(j  ,i+1)-MY(j  ,i  ))...
    - wminusy(MY,j,i,numinusy,TVD).*(MY(j  ,i  )-MY(j-1,i  ))...
    + wplusy (MY,j,i,nuplusy ,TVD).*(MY(j+1,i  )-MY(j  ,i  ));

% Energy Equation
% nuplusy  = abs(MYy(j-1,i-2)./RHOy(j-1,i-2)...
%     +sqrt(gamma*Py(j-1,i-2)./RHOy(j-1,i-2)))*Cy;
% nuplusx  = abs(MXx(j-2,i-1)./RHOx(j-2,i-1)...
%     +sqrt(gamma*Px(j-2,i-1)./RHOx(j-2,i-1)))*Cx;
% numinusy = abs(MYy(j-2,i-2)./RHOy(j-2,i-2)...
%     +sqrt(gamma*Py(j-2,i-2)./RHOy(j-2,i-2)))*Cy;
% numinusx = abs(MXx(j-2,i-2)./RHOx(j-2,i-2)...
%     +sqrt(gamma*Px(j-2,i-2)./RHOx(j-2,i-2)))*Cx;

E_New   (j,i) = E(j,i) - Cx*(((Ex(j-2,i-1)+Px(j-2,i-1)).*MXx(j-2,i-1)./RHOx(j-2,i-1)) ...
                           - ((Ex(j-2,i-2)+Px(j-2,i-2)).*MXx(j-2,i-2)./RHOx(j-2,i-2)))...
                       - Cy*(((Ey(j-1,i-2)+Py(j-1,i-2)).*MYy(j-1,i-2)./RHOy(j-1,i-2)) ...
                           - ((Ey(j-2,i-2)+Py(j-2,i-2)).*MYy(j-2,i-2)./RHOy(j-2,i-2)))...
    - wminusx(E,j,i,numinusx,TVD).*(E(j  ,i  )-E(j  ,i-1))...
    + wplusx (E,j,i,nuplusx ,TVD).*(E(j  ,i+1)-E(j  ,i  ))...
    - wminusy(E,j,i,numinusy,TVD).*(E(j  ,i  )-E(j-1,i  ))...
    + wplusy (E,j,i,nuplusy ,TVD).*(E(j+1,i  )-E(j  ,i  ));

% Pressure
P       (j,i) = (gamma-1)*(E_New(j,i)-0.5*((MX_New(j,i).^2 ...
    + MY_New(j,i).^2)./RHO_New(j,i)));


RHO(j,i) = RHO_New(j,i);
MY (j,i) = MY_New (j,i);
MX (j,i) = MX_New (j,i);
E  (j,i) = E_New  (j,i);

plot1(1:length(j)  ,1) = RHO(j,nx/2);
if BC == 2, plot1(length(j)+j-2,1) = RHO(j,nx/2); end

if mod(t,1*dt) == 0
    figure(1),surf(X,Y,RHO(j,i))  , title('RHO')
    figure(2),surf(X,Y,MY (j,i))  , title('MY' )
    figure(3),plot(plot1,'.-') , title('RHO')
    figure(4),surf(X,Y,MY(j,i)./RHO(j,i))  , title('V')
    figure(5),surf(X,Y,MX(j,i)./RHO(j,i))  , title('X')
end

% if t>2.54, break; end
telapsed = toc(t1)
end