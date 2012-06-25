clear all; clc; close all
% SHALLOW_WATER
% MATLAB version of: WAVE -- 2D Shallow Water Model
% New Mexico Supercomputing Challenge, Glorieta Kickoff 2007
%
% Lax-Wendroff finite difference method.
% Reflective boundary conditions.
% Random water drops initiate gravity waves.
% Surface plot displays height colored by momentum.
% Plot title shows t = simulated time and tv = a measure of total variation.
% An exact solution to the conservation law would have constant tv.
% Lax-Wendroff produces nonphysical oscillations and increasing tv.
%
% Cleve Moler, The MathWorks, Inc.
% Derived from C programs by
%    Bob Robey, Los Alamos National Laboratory.
%    Joseph Koby, Sarah Armstrong, Juan-Antonio Vigil, Vanessa Trujillo, McCurdy School.
%    Jonathan Robey, Dov Shlachter, Los Alamos High School.
% See:
%    http://en.wikipedia.org/wiki/Shallow_water_equations
%    http://www.amath.washington.edu/~rjl/research/tsunamis
%    http://www.amath.washington.edu/~dgeorge/tsunamimodeling.html
%    http://www.amath.washington.edu/~claw/applications/shallow/www

% Parameters

n  = 64;                  % grid size
g  = 9.8;                 % gravitational constant
dt = 0.02;                % hardwired timestep
dx = 1.0;
dy = 1.0;
m  = 3;                   % number of materials
Cx = dt/dx;
Cy = dt/dy;

nplotstep = 8;            % plot interval
ndrops = 5;               % maximum number of drops
dropstep = 500;           % drop interval
D = droplet (2.0,11);     % simulate a water drop

% Initialize graphics

[surfplot,phiplot,heightplot,top,start,stop] = initgraphics (n,m);

% Outer loop, restarts.

while get (stop,'value') == 0
   set (start,'value',0)
   
   H = ones (n+4,n+4);   U = zeros (n+4,n+4);  V = zeros (n+4,n+4);
   Hx = zeros (n+3,n+2); Ux = zeros (n+3,n+2); Vx = zeros (n+3,n+2);
   Hy = zeros (n+2,n+3); Uy = zeros (n+2,n+3); Vy = zeros (n+2,n+3);
   phi = zeros(n+4,n+4,m);
   phix = zeros(n+3,n+2,m);
   phiy = zeros(n+2,n+3,m);
   ndrop = ceil (rand*ndrops);
   nstep = 0;
   
   phi(:,:,1) = 1;

   % Inner loop, time steps.

   while get (start,'value')==0 && get (stop,'value')==0
       nstep = nstep + 1;
       
       % Initial Dam Break
       if nstep == 1
           w = size (D,1);
           hconc = 1;
           
           % Concentration of materials added
           % Material 1 added
           i = ceil (.5 *(n-w))+(1:ceil(w/2));
           j = ceil (.5 *(n-w))+(1:w);
           S = ones(ceil(length(D)/2),length(D));
           phi(i,j,2) = hconc*S;
           
           % Material 2 added
           i = ceil (.5 *(n-w))+(ceil(w/2)+1:w+1);
           j = ceil (.5 *(n-w))+(1:w);
           S = ones(ceil(length(D)/2),length(D));
           phi(i,j,3) = hconc*S;
           
           phi(:,:,1) = phi(:,:,1) - phi(:,:,2) - phi(:,:,3);
           
           i = ceil (.5 *(n-w))+(1:w);
           j = ceil (.5 *(n-w))+(1:w);
           S = ones(size(D));
           hdrop = 3;
           H (i,j) = H (i,j) + hdrop*S;
       end

       
%        % Random water drops
%        if mod (nstep,dropstep) == 0 && nstep <= ndrop*dropstep
%            w = size (D,1);
%            i = ceil (rand*(n-w))+(1:w);
%            j = ceil (rand*(n-w))+(1:w);
%            hdrop = rand;
%            phi (i,j) = (hdrop*D + phi(i,j) .* H(i,j))./(hdrop*D + H(i,j));
%            S = ones(size(D));
%            H (i,j) = H (i,j) + hdrop*S;
%        end


       % Reflective boundary conditions
       H (:,2) = H (:,3);      U (:,2) = U (:,3);       V (:,2) = -V(:,3);
       H (:,n+3) = H (:,n+2);  U (:,n+3) = U (:,n+2);   V (:,n+3) = -V(:,n+2);
       H (2,:) = H (3,:);      U (2,:) = -U(3,:);      V (2,:) = V (3,:);
       H (n+3,:) = H (n+2,:);  U (n+3,:) = -U(n+2,:);  V (n+3,:) = V (n+2,:);
       phi (:,2,:) = phi (:,3,:);      
       phi (:,n+3,:) = phi (:,n+2,:);  
       phi (2,:,:) = phi (3,:,:);     
       phi (n+3,:,:) = phi (n+2,:,:); 
       
       H (:,1) = H (:,2);      U (:,1) = U (:,2);       V (:,1) = V(:,2);
       H (:,n+4) = H (:,n+3);  U (:,n+4) = U (:,n+3);   V (:,n+4) = V(:,n+3);
       H (1,:) = H (2,:);      U (1,:) = U(2,:);      V (1,:) = V (2,:);
       H (n+4,:) = H (n+3,:);  U (n+4,:) = U(n+3,:);  V (n+4,:) = V (n+3,:);
       phi (:,1,:) = phi (:,2,:);      
       phi (:,n+4,:) = phi (:,n+3,:);  
       phi (1,:,:) = phi (2,:,:);     
       phi (n+4,:,:) = phi (n+3,:,:); 
       
       % First half step
       % x direction
       i = 2:n+2;
       j = 2:n+1;
       
       % phi
       for k = 1:m
       phix(i,j,k) = (phi(i+1,j+1,k) + phi(i,j+1,k))/2 - ...
           dt/(4*dx) * (U(i+1,j+1)./H(i+1,j+1)+U(i,j+1)./H(i,j+1)) .*...
           (phi(i+1,j+1,k)-phi(i,j+1,k));
       end
       
       % height
       Hx (i,j) = (H (i+1,j+1)+H (i,j+1))/2 ...
           - dt/(2*dx)*(U (i+1,j+1)-U (i,j+1));
   
       % x momentum
       Ux (i,j) = (U (i+1,j+1)+U (i,j+1))/2 -  ...
           dt/(2*dx)*((U (i+1,j+1).^2./H(i+1,j+1) + g/2*H (i+1,j+1).^2) - ...
           (U (i,j+1).^2./H(i,j+1) + g/2*H (i,j+1).^2));
       
       % y momentum
       Vx (i,j) = (V (i+1,j+1)+V (i,j+1))/2 - ...
           dt/(2*dx)*((U (i+1,j+1).*V (i+1,j+1)./H(i+1,j+1)) - ...
           (U (i,j+1).*V (i,j+1)./H(i,j+1)));
       
       % y direction
       i = 2:n+1;
       j = 2:n+2;
      
       % phi
       for k=1:m
       phiy(i,j,k) = (phi(i+1,j+1,k) + phi(i+1,j,k))/2 -...
           dt/(4*dy)*(V(i+1,j+1)./H(i+1,j+1)+V(i+1,j)./H(i+1,j)) .*...
           (phi(i+1,j+1,k)-phi(i+1,j,k));
       end
       
       % height
       Hy (i,j) = (H (i+1,j+1)+H (i+1,j))/2 - dt/(2*dy)*(V (i+1,j+1)-V (i+1,j));
       
       % x momentum
       Uy (i,j) = (U (i+1,j+1)+U (i+1,j))/2 - ...
           dt/(2*dy)*((V (i+1,j+1).*U (i+1,j+1)./H(i+1,j+1)) - ...
           (V (i+1,j).*U (i+1,j)./H(i+1,j)));
       % y momentum
       Vy (i,j) = (V (i+1,j+1)+V (i+1,j))/2 - ...
           dt/(2*dy)*((V (i+1,j+1).^2./H(i+1,j+1) + g/2*H (i+1,j+1).^2) - ...
           (V (i+1,j).^2./H(i+1,j) + g/2*H (i+1,j).^2));
   
       % Second half step
       i = 3:n+2;
       j = 3:n+2;
       
       
       % concentration
       nuplusy  = abs(Vy(i-1,j  )./Hy(i-1,j  ))*Cy;
       nuplusx  = abs(Ux(i  ,j-1)./Hx(i  ,j-1))*Cx;
       numinusy = abs(Vy(i-1,j-1)./Hy(i-1,j-1))*Cy;
       numinusx = abs(Ux(i-1,j-1)./Hx(i-1,j-1))*Cx;

       for k=1:m
       phi (i,j,k) = phi (i,j,k) - (dt/(2*dx))*(Ux (i,j-1)./Hx(i,j-1)...
           + Ux (i-1,j-1)./Hx(i-1,j-1)) .* (phix(i,j-1,k)-phix(i-1,j-1,k))...
           - (dt/(2*dy))*(Vy (i-1,j)./Hy(i-1,j)...
           + Vy (i-1,j-1)./Hy(i-1,j-1)) .* (phiy(i-1,j,k)-phiy(i-1,j-1,k))...
           - wminusx(phi(:,:,k),i,j,numinusx,2).*(phi(i  ,j  ,k)-phi(i-1,j  ,k))...
           + wplusx (phi(:,:,k),i,j,nuplusx ,2).*(phi(i+1,j  ,k)-phi(i  ,j  ,k))...
           - wminusy(phi(:,:,k),i,j,numinusy,2).*(phi(i  ,j  ,k)-phi(i  ,j-1,k))...
           + wplusy (phi(:,:,k),i,j,nuplusy ,2).*(phi(i  ,j+1,k)-phi(i  ,j  ,k));
       end
       
       % height
       nuplusy  = abs(abs(Vy(i-1,j  )./Hy(i-1,j  )) - sqrt(g*Hy(i-1,j  )))*Cy;
       nuplusx  = abs(abs(Ux(i  ,j-1)./Hx(i  ,j-1)) - sqrt(g*Hx(i  ,j-1)))*Cx;
       numinusy = abs(abs(Vy(i-1,j-1)./Hy(i-1,j-1)) - sqrt(g*Hy(i-1,j-1)))*Cy;
       numinusx = abs(abs(Ux(i-1,j-1)./Hx(i-1,j-1)) - sqrt(g*Hx(i-1,j-1)))*Cx;
       
       H (i,j) = H (i,j)...
           - (dt/dx)*(Ux (i,j-1)-Ux (i-1,j-1))...
           - (dt/dy)*(Vy (i-1,j)-Vy (i-1,j-1))...
           - wminusx(H,i,j,numinusx,1).*(H(i  ,j  )-H(i-1,j  ))...
           + wplusx (H,i,j,nuplusx ,1).*(H(i+1,j  )-H(i  ,j  ))...
           - wminusy(H,i,j,numinusy,1).*(H(i  ,j  )-H(i  ,j-1))...
           + wplusy (H,i,j,nuplusy ,1).*(H(i  ,j+1)-H(i  ,j  ));
       
       % x momentum
       nuplusy  = (abs(Vy(i-1,j  )./Hy(i-1,j  )) + sqrt(g*Hy(i-1,j  )))*Cy;
       nuplusx  = (abs(Ux(i  ,j-1)./Hx(i  ,j-1)) + sqrt(g*Hx(i  ,j-1)))*Cx;
       numinusy = (abs(Vy(i-1,j-1)./Hy(i-1,j-1)) + sqrt(g*Hy(i-1,j-1)))*Cy;
       numinusx = (abs(Ux(i-1,j-1)./Hx(i-1,j-1)) + sqrt(g*Hx(i-1,j-1)))*Cx;
       
       U (i,j) = U (i,j)...
           - (dt/dx)*((Ux(i,j-1).^2./Hx(i,j-1) + g/2*Hx(i,j-1).^2)...
           - (Ux (i-1,j-1).^2./Hx(i-1,j-1) + g/2*Hx (i-1,j-1).^2))...
           - (dt/dy)*(Vy (i-1,j).*Uy (i-1,j)./Hy(i-1,j)...
           - Vy (i-1,j-1).*Uy (i-1,j-1)./Hy(i-1,j-1))...
           - wminusx(U,i,j,numinusx,1).*(U(i  ,j)-U(i-1,j  ))...
           + wplusx (U,i,j,nuplusx ,1).*(U(i+1,j)-U(i  ,j  ))...
           - wminusy(U,i,j,numinusy,1).*(U(i  ,j)-U(i  ,j-1))...
           + wplusy (U,i,j,nuplusy ,1).*(U(i,j+1)-U(i  ,j  ));
       
       % y momentum
       nuplusy  = (abs(Vy(i-1,j  )./Hy(i-1,j  )) + sqrt(g*Hy(i-1,j  )))*Cy;
       nuplusx  = (abs(Ux(i  ,j-1)./Hx(i  ,j-1)) + sqrt(g*Hx(i  ,j-1)))*Cx;
       numinusy = (abs(Vy(i-1,j-1)./Hy(i-1,j-1)) + sqrt(g*Hy(i-1,j-1)))*Cy;
       numinusx = (abs(Ux(i-1,j-1)./Hx(i-1,j-1)) + sqrt(g*Hx(i-1,j-1)))*Cx;
       
       V (i,j) = V (i,j)...
           - (dt/dx)*(Ux(i,j-1).*Vx(i,j-1)./Hx(i,j-1)...
           - Ux(i-1,j-1).*Vx(i-1,j-1)./Hx(i-1,j-1))...
           - (dt/dy)*((Vy(i-1,j).^2./Hy(i-1,j)+g/2*Hy (i-1,j).^2)...
           - (Vy(i-1,j-1).^2./Hy(i-1,j-1)+g/2*Hy(i-1,j-1).^2))...
           - wminusx(V,i,j,numinusx,1).*(V(i  ,j)-V(i-1,j  ))...
           + wplusx (V,i,j,nuplusx ,1).*(V(i+1,j)-V(i  ,j  ))...
           - wminusy(V,i,j,numinusy,1).*(V(i  ,j)-V(i  ,j-1))...
           + wplusy (V,i,j,nuplusy ,1).*(V(i,j+1)-V(i  ,j  ));
   
       % Update plot
       if mod (nstep,nplotstep) == 0
          figure(1);
          C = abs (U (i,j)) + abs (V (i,j));  % Color shows momemtum
          t = nstep*dt;
          tv = norm (C,'fro');
          set (surfplot,'zdata',H (i,j),'cdata',C);
          set (top,'string',sprintf ('t = %6 .2f,  tv = %6 .2f',t,tv))
          drawnow
          
          figure(2);
          set (phiplot,'zdata',H (i,j),'cdata',phi(i,j,1));
          drawnow
          
          for fig = 3:m+2
              figure(fig);
              set(heightplot(fig), 'zdata', phi(i,j,fig-2).*H(i,j),'cdata',...
                  phi(i,j,fig-2).*H(i,j));
              drawnow
          end
          
%           figure(5)
%           plot(phi(i,length(phi)/2))
       end
      
       if any (any (isnan (H))), break, end  % Unstable, restart
       if any (any (isinf (H))), break, end  % Unstable, restart
   end
   if any (any (isnan (H))), break, end  % Unstable, restart
   if any (any (isinf (H))), break, end  % Unstable, restart
end
% close (figure(1)), close(figure(2))
