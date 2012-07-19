clear all;
%Mesh Generation
nx = 64;
ny = 64;
Lx = 1;
Ly = 1;
xp = linspace(0,Lx,nx);
yp = linspace(0,Ly,ny);
[X,Y] = meshgrid(xp,yp);
X = reshape(X,nx*ny,1);
Y = reshape(Y,nx*ny,1);
%tri = delaunay(X,Y);
n = 1;
for i = 1:nx*ny
    if rem(i,ny) ~= 0 && i < (nx-1)*ny
        tri(n,1:3) = [i, i+1, i+ny];
        n = n+1;
        tri(n,1:3) = [i+1, i+1+ny, i+ny];
        n = n+1;
    end
end
%trimesh(tri,X,Y,zeros(length(X),1));
%Boundary Conditions
BC = 1; %0 - stationary boundary %1 - periodic

%Variable initialization
u = -1;
v = 0;
a = 1;
phip = zeros(nx,ny);

i = floor(1+0*nx/8):floor(1+7.9*nx/8);
j = floor(1+3*nx/8):floor(1+5*nx/8);
phip(i,j) = 1;

phi = reshape(phip,nx*ny,1);
%trimesh(tri,X,Y,phi);
phi_plot = zeros(nx*ny,1);

%Assemble the Global stiffness matrix from local matrices
n = 1:length(tri);
Me1 = zeros(3,3);
Me2 = zeros(3,3);
Me3 = zeros(3,3);
h(n) = .707*Lx/(nx-1);
A(n) = abs(.5*(X(tri(n,1)).*Y(tri(n,2)) - X(tri(n,1)).*Y(tri(n,3))...
             + X(tri(n,2)).*Y(tri(n,3)) - X(tri(n,2)).*Y(tri(n,1))...
             + X(tri(n,3)).*Y(tri(n,1)) - X(tri(n,3)).*Y(tri(n,2))))';
u = -u;
v = -v;
%This is where the matrix calculation should go
%
%END matrix calc

dt = .1;
nt = 1000;
for t=0:dt:dt*nt
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    M = zeros(nx*ny,nx*ny);
K = zeros(nx*ny,nx*ny);
    U = 0*exp(-1000*(X-.2).^2) + 1;
V = 0*eps*(randi([-1,1],length(X),1));
%figure(2);
%plot3(X,Y,V,'.');
for n = 1:length(tri)
    uv_mag = sqrt(u^2+v^2);
    if rem(n,2) == 1
        u = (U(tri(n,3))+U(tri(n,2))+U(tri(n,1))+U(tri(n+1,2)))/4;
        v = (V(tri(n,3))+V(tri(n,2))+V(tri(n,1))+V(tri(n+1,2)))/4;
    else
        u = (U(tri(n,3))+U(tri(n,2))+U(tri(n,1))+U(tri(n-1,1)))/4;
        v = (V(tri(n,3))+V(tri(n,2))+V(tri(n,1))+V(tri(n-1,1)))/4;
    end
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

%Set to zero at the boundaries
if BC == 0
    b = [1:ny,ny+1:ny:nx*ny-(2*ny-1),nx*ny-(ny-1):nx*ny,2*ny:ny:(nx-1)*ny];
    nb = 1:length(X); nb(b) = [];
    %plot(X(nb),Y(nb),'.');
    K(b,:) = [];
    K(:,b) = [];
    M(b,:) = [];
    M(:,b) = [];
    phi(b) = [];
    
    %Periodic BCs
elseif BC == 1
    alias = [(1:ny)',(nx*ny-(ny-1):nx*ny)';(1:ny:nx*ny-(2*ny-1))',(ny:ny:(nx-1)*ny)'];
    nalias = 1:nx*ny;
    nalias(alias(:,2)) = [];
    
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
    phi(alias(:,2)) = [];
end

M = sparse(M);
K = sparse(K);

options = odeset('Mass',M);
f = @(t, phi)(-K * phi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[T,phinew] = ode45(f,[0 dt],phi,options);
phi = phinew(end,1:end)';
if BC == 0
    phi_plot(nb) = phi;
elseif BC == 1
    temp = zeros((ny)*(nx),1);
    temp(nalias) = phi;
    
    temp(alias(:,2)) = temp(alias(:,1)); %Having trouble with repeats
    temp(alias(:,2)) = temp(alias(:,1));
    phi_plot = temp;
    phi = phi_plot;
end
figure(1);
%plot3(X,Y,phi_plot,'.');
%trimesh(tri,X,Y,phi_plot);
surf(reshape(X,nx,ny),reshape(Y,nx,ny),reshape(phi_plot,nx,ny));
%plot(Y(2*ny:3*ny),phi_plot(2*ny:3*ny),'.');
drawnow;
end