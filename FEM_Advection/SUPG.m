clear all;
%Mesh Generation
nx = 10;
ny = 10;
Lx = 1;
Ly = 1;
xp = linspace(0,Lx,nx);
yp = linspace(0,Ly,ny);
[X,Y] = meshgrid(xp,yp);
X = reshape(X,nx*ny,1);
Y = reshape(Y,nx*ny,1);
tri = delaunay(X,Y);
b = [1:ny,ny+1:ny:nx*ny-(2*ny-1),nx*ny-(ny-1):nx*ny,2*ny:ny:(nx-1)*ny];
nb = 1:length(X); nb(b) = [];
plot(X(nb),Y(nb),'.');

%Variable initialization
u = 0;
v = .2;
phip = zeros(nx,ny);
i = floor(1+0*nx/8):floor(1+3*nx/8);
j = floor(1+0*nx/8):floor(1+7.9*nx/8);
phip(i,j) = 1;
a = 1;
phi = reshape(phip,nx*ny,1);
trimesh(tri,X,Y,phi);
M = zeros(nx*ny,nx*ny);
K = zeros(nx*ny,nx*ny);
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
for n = 1:length(tri)
    Me1(1:3,1:3) = A(n)*[1/6,1/12,1/12;1/12,1/6,1/12;1/12,1/12,1/6];
    Me2(1:3,1:3) = (a*h(n)*sign(u)/12)*...
        [Y(tri(n,2))-Y(tri(n,3)),Y(tri(n,2))-Y(tri(n,3)),Y(tri(n,2))-Y(tri(n,3));... 
         Y(tri(n,3))-Y(tri(n,1)),Y(tri(n,3))-Y(tri(n,1)),Y(tri(n,3))-Y(tri(n,1));...
         Y(tri(n,1))-Y(tri(n,2)),Y(tri(n,1))-Y(tri(n,2)),Y(tri(n,1))-Y(tri(n,2))];
    Me3(1:3,1:3) = (a*h(n)*sign(v)/12)*...
        [X(tri(n,3))-X(tri(n,2)),X(tri(n,3))-X(tri(n,2)),X(tri(n,3))-X(tri(n,2));... 
         X(tri(n,1))-X(tri(n,3)),X(tri(n,1))-X(tri(n,3)),X(tri(n,1))-X(tri(n,3));...
         X(tri(n,2))-X(tri(n,1)),X(tri(n,2))-X(tri(n,1)),X(tri(n,2))-X(tri(n,1))];
    M(tri(n,1),tri(n,:)) = M(tri(n,1),tri(n,:)) + Me1(1,:) + Me2(1,:) + Me3(1,:);
    M(tri(n,2),tri(n,:)) = M(tri(n,2),tri(n,:)) + Me1(2,:) + Me2(2,:) + Me3(2,:);
    M(tri(n,3),tri(n,:)) = M(tri(n,3),tri(n,:)) + Me1(3,:) + Me2(3,:) + Me3(3,:);
    Ke2r1(1,1:3) = ((u/6) + (u*a*h(n)/(8*A(n))) * (sign(u)*(Y(tri(n,2))-Y(tri(n,3))) + sign(v)*(X(tri(n,3))-X(tri(n,2))))) * [Y(tri(n,2))-Y(tri(n,3)), Y(tri(n,3))-Y(tri(n,1)), Y(tri(n,1))-Y(tri(n,2))]... 
                 + ((v/6) + (v*a*h(n)/(8*A(n))) * (sign(u)*(Y(tri(n,2))-Y(tri(n,3))) + sign(v)*(X(tri(n,3))-X(tri(n,2))))) * [X(tri(n,3))-X(tri(n,2)), X(tri(n,1))-X(tri(n,3)), X(tri(n,2))-X(tri(n,1))];
    Ke2r2(1,1:3) = ((u/6) + (u*a*h(n)/(8*A(n))) * (sign(u)*(Y(tri(n,3))-Y(tri(n,1))) + sign(v)*(X(tri(n,1))-X(tri(n,3))))) * [Y(tri(n,2))-Y(tri(n,3)), Y(tri(n,3))-Y(tri(n,1)), Y(tri(n,1))-Y(tri(n,2))]... 
                 + ((v/6) + (v*a*h(n)/(8*A(n))) * (sign(u)*(Y(tri(n,3))-Y(tri(n,1))) + sign(v)*(X(tri(n,1))-X(tri(n,3))))) * [X(tri(n,3))-X(tri(n,2)), X(tri(n,1))-X(tri(n,3)), X(tri(n,2))-X(tri(n,1))];
    Ke2r3(1,1:3) = ((u/6) + (u*a*h(n)/(8*A(n))) * (sign(u)*(Y(tri(n,1))-Y(tri(n,2))) + sign(v)*(X(tri(n,2))-X(tri(n,1))))) * [Y(tri(n,2))-Y(tri(n,3)), Y(tri(n,3))-Y(tri(n,1)), Y(tri(n,1))-Y(tri(n,2))]... 
                 + ((v/6) + (v*a*h(n)/(8*A(n))) * (sign(u)*(Y(tri(n,1))-Y(tri(n,2))) + sign(v)*(X(tri(n,2))-X(tri(n,1))))) * [X(tri(n,3))-X(tri(n,2)), X(tri(n,1))-X(tri(n,3)), X(tri(n,2))-X(tri(n,1))];
    K(tri(n,1),tri(n,:)) = K(tri(n,1),tri(n,:)) + Ke2r1(1,1:3);
    K(tri(n,2),tri(n,:)) = K(tri(n,2),tri(n,:)) + Ke2r2(1,1:3);
    K(tri(n,3),tri(n,:)) = K(tri(n,3),tri(n,:)) + Ke2r3(1,1:3);
end

%Set to zero at the boundaries
K(b,:) = [];
K(:,b) = [];
M(b,:) = [];
M(:,b) = [];
phi(b) = [];

M = sparse(M);
K = sparse(K);

options = odeset('Mass',M);
f = @(t, phi)(-K * phi);
dt = .1;
nt = 30;
for t=0:dt:dt*nt
    [T,phinew] = ode45(f,[0 dt],phi,options);
    phi = phinew(end,1:end)';
    phi_plot(nb) = phi;
    %plot3(X,Y,phi_plot,'.');
    %trimesh(tri,X,Y,phi_plot);
    %plot(Y(2*ny:3*ny),phi_plot(2*ny:3*ny),'.');
    drawnow;
end