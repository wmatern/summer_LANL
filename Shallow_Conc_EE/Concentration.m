clear all; close all; clc;

% Create Mesh
xmin = -2; dx = 0.1; xmax = 2;
ymin = -2; dy = 0.1; ymax = 2;

[X, Y] = meshgrid(xmin:dx:ymax,ymin:dx:ymax);

% Timestep and maximum time
tmax = 10; dt = 0.001;

% Set Gravity
g = 9.8;

h = zeros(length(X(:,1))+4,length(Y(:,1))+4);
for i=1:length(X(:,1))
    for j=1:length(Y(:,1))
        h(i+2,j+2) = exp(-10*(X(i,j)^2+Y(i,j)^2))+1;
    end
end

phi = zeros(length(X(:,1))+4,length(Y(:,1))+4);
for i=1:length(X(:,1))
    for j=1:length(Y(:,1))
        phi(i+2,j+2) = exp(-10*(X(i,j)^2+Y(i,j)^2));
    end
end

u = zeros(length(X(:,1))+4,length(Y(:,1))+4);
v = zeros(length(X(:,1))+4,length(Y(:,1))+4);

% figure(1),  contour(X,Y,phi(3:length(X(:,1))+2,3:length(Y(1,:))+2));
k=0;
% Solve ODE
for t = 0:dt:tmax
    % BCs for h
    h(1,:) = h(length(X(:,1))+1,:               );
    h(2,:) = h(length(X(:,1))+2,:               );
    h(:,1) = h(:               ,length(Y(:,1))+1);
    h(:,2) = h(:               ,length(Y(:,1))+2);
    
    h(length(X(:,1))+3,:               ) = h(3,:);
    h(length(X(:,1))+4,:               ) = h(4,:);
    h(:               ,length(Y(:,1))+3) = h(:,3);
    h(:               ,length(Y(:,1))+4) = h(:,4);
    
    % BCs for u
    u(1,:) = u(length(X(:,1))+1,:               );
    u(2,:) = u(length(X(:,1))+2,:               );
    u(:,1) = u(:               ,length(Y(:,1))+1);
    u(:,2) = u(:               ,length(Y(:,1))+2);
    
    u(length(X(:,1))+3,:               ) = u(3,:);
    u(length(X(:,1))+4,:               ) = u(4,:);
    u(:               ,length(Y(:,1))+3) = u(:,3);
    u(:               ,length(Y(:,1))+4) = u(:,4);
    
    % BCs for v
    v(1,:) = v(length(X(:,1))+1,:               );
    v(2,:) = v(length(X(:,1))+2,:               );
    v(:,1) = v(:               ,length(Y(:,1))+1);
    v(:,2) = v(:               ,length(Y(:,1))+2);
    
    v(length(X(:,1))+3,:               ) = v(3,:);
    v(length(X(:,1))+4,:               ) = v(4,:);
    v(:               ,length(Y(:,1))+3) = v(:,3);
    v(:               ,length(Y(:,1))+4) = v(:,4);
    
    % BCs for phi
    phi(1,:) = phi(length(X(:,1))+1,:               );
    phi(2,:) = phi(length(X(:,1))+2,:               );
    phi(:,1) = phi(:               ,length(Y(:,1))+1);
    phi(:,2) = phi(:               ,length(Y(:,1))+2);
    
    phi(length(X(:,1))+3,:               ) = phi(3,:);
    phi(length(X(:,1))+4,:               ) = phi(4,:);
    phi(:               ,length(Y(:,1))+3) = phi(:,3);
    phi(:               ,length(Y(:,1))+4) = phi(:,4);

    hu = h.*u; hv = h.*v;
    hu2 = h.*u.^2+0.5*g*h.^2; hv2 = h.*v.^2+0.5*g*h.^2;
    huv = h.*u.*v;
    for i = 3:length(X(:,1))+2
        for j = 3:length(Y(1,:))+2
            h_new (i,j) = h(i,j) +dt*(d_dx(hu,i,j)/dx+d_dy(hv,i,j)/dy);
            hu_new(i,j) = hu(i,j)+dt*(d_dx(hu2,i,j)/dx+d_dy(huv,i,j)/dy);
            hv_new(i,j) = hv(i,j)+dt*(d_dx(huv,i,j)/dx+d_dy(hv2,i,j)/dy);
            phi   (i,j) = phi(i,j)+dt*(u(i,j)/dx*d_dx(phi,i,j)...
                            +v(i,j)/dy*d_dy(phi,i,j));
        end
    end
    h = h_new;
    u = hu_new./h;
    v = hv_new./h;
    k=k+1;
    if mod(k,5)==0
 %       figure(2), contour(X,Y,phi(3:length(X(:,1))+2,3:length(Y(1,:))+2));
        figure(3), surf   (X,Y,phi(3:length(X(:,1))+2,3:length(Y(1,:))+2));
        axis([-2 2 -2 2 -2 2])
        figure(4), surf   (X,Y,h(3:length(X(:,1))+2,3:length(Y(1,:))+2));
        axis([-2 2 -2 2 0 2])
    end
end