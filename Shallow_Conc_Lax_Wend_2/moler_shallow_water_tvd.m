function shallow_water
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

n = 64;                  % grid size
g = 9.8;                 % gravitational constant
dt = 0.02;               % hardwired timestep
dx = 1.0;
dy = 1.0;
nplotstep = 8;           % plot interval
ndrops = 5;              % maximum number of drops
dropstep = 500;          % drop interval
D = droplet (2.0,11);     % simulate a water drop

% Initialize graphics

[surfplot,phiplot,uplot,vplot,top,start,stop] = initgraphics (n);

% Outer loop, restarts.

while get (stop,'value') == 0
   set (start,'value',0)
   
   H = ones (n+4,n+4);   U = zeros (n+4,n+4);  V = zeros (n+4,n+4);
   Hx = zeros (n+3,n+2); Ux = zeros (n+3,n+2); Vx = zeros (n+3,n+2);
   Hy = zeros (n+2,n+3); Uy = zeros (n+2,n+3); Vy = zeros (n+2,n+3);
   phi = zeros(n+4,n+4);
   phix = zeros(n+3,n+2);
   phiy = zeros(n+2,n+3);
   ndrop = ceil (rand*ndrops);
   nstep = 0;

   % Inner loop, time steps.

   while get (start,'value')==0 && get (stop,'value')==0
       nstep = nstep + 1;

       % Random water drops
       if mod (nstep,dropstep) == 0 && nstep <= ndrop*dropstep
           w = size (D,1);
           i = ceil (rand*(n-w))+(1:w);
           j = ceil (rand*(n-w))+(1:w);
           hdrop = rand;
           phi (i,j) = (hdrop*D + phi(i,j) .* H(i,j))./(hdrop*D + H(i,j));
           S = ones(size(D));
           H (i,j) = H (i,j) + hdrop*S;
       end
     
       % Reflective boundary conditions
       H (:,2) = H (:,3);      U (:,2) = U (:,3);       V (:,2) = -V(:,3);
       H (:,n+3) = H (:,n+2);  U (:,n+3) = U (:,n+2);   V (:,n+3) = -V(:,n+2);
       H (2,:) = H (3,:);      U (2,:) = -U(3,:);      V (2,:) = V (3,:);
       H (n+3,:) = H (n+2,:);  U (n+3,:) = -U(n+2,:);  V (n+3,:) = V (n+2,:);
       phi (:,2) = phi (:,3);      
       phi (:,n+3) = phi (:,n+2);  
       phi (2,:) = phi (3,:);     
       phi (n+3,:) = phi (n+2,:); 
       
       H (:,1) = H (:,2);      U (:,1) = U (:,2);       V (:,1) = V(:,2);
       H (:,n+4) = H (:,n+3);  U (:,n+4) = U (:,n+3);   V (:,n+4) = V(:,n+3);
       H (1,:) = H (2,:);      U (1,:) = U(2,:);      V (1,:) = V (2,:);
       H (n+4,:) = H (n+3,:);  U (n+4,:) = U(n+3,:);  V (n+4,:) = V (n+3,:);
       phi (:,1) = phi (:,2);      
       phi (:,n+4) = phi (:,n+3);  
       phi (1,:) = phi (2,:);     
       phi (n+4,:) = phi (n+3,:); 
       
       % First half step
   
       % x direction
       i = 2:n+2;
       j = 2:n+1;
   
       % phi
       phix(i,j) = (phi(i+1,j+1) + phi(i,j+1))/2 - dt/(4*dx) * (U(i+1,j+1)./H(i+1,j+1)+U(i,j+1)./H(i,j+1)) .* (phi(i+1,j+1)-phi(i,j+1));
       
       % height
       Hx (i,j) = (H (i+1,j+1)+H (i,j+1))/2 - dt/(2*dx)*(U (i+1,j+1)-U (i,j+1));
   
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
       phiy(i,j) = (phi(i+1,j+1) + phi(i+1,j))/2 - dt/(4*dy)*(V(i+1,j+1)./H(i+1,j+1)+V(i+1,j)./H(i+1,j)) .* (phi(i+1,j+1)-phi(i+1,j));
       
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
   
       %TVD
       duminus1 = H(i-1,j) - H(i-2,j);
       duminus2 = U(i-1,j) - U(i-2,j);
       duplus1 = H(i+1,j) - H(i,j);
       duplus2 = U(i+1,j) - U(i,j);
       duhalf1 = H(i,j) - H(i-1,j);
       duhalf2 = U(i,j) - U(i-1,j);
       rdenom = max(duhalf1.^2 + duhalf2.^2, eps);
       rnumplus = duplus1.*duhalf1 + duplus2.*duhalf2;
       rnumminus = duminus1.*duhalf1 + duminus2.*duhalf2;
       rplus = rnumplus./rdenom;
       rminus = rnumminus./rdenom;
       q = max(min(1.0, min(rminus, rplus)),0.0);
       nu = abs(Ux(i-2,j-2)) + sqrt(g*Hx(i-2,j-2))*dt/dx;
       cv = nu.*(1-nu);
       wminusx = .5*cv.*(1-q);
       
       duminus1 = H(i,j) - H(i-1,j);
       duminus2 = U(i,j) - U(i-1,j);
       duplus1 = H(i+2,j) - H(i+1,j);
       duplus2 = U(i+2,j) - U(i+1,j);
       duhalf1 = H(i+1,j) - H(i,j);
       duhalf2 = U(i+1,j) - U(i,j);
       rdenom = max(duhalf1.^2 + duhalf2.^2, eps);
       rnumplus = duplus1.*duhalf1 + duplus2.*duhalf2;
       rnumminus = duminus1.*duhalf1 + duminus2.*duhalf2;
       rplus = rnumplus./rdenom;
       rminus = rnumminus./rdenom;
       q = max(min(1.0, min(rminus, rplus)),0.0);
       nu = abs(Ux(i-1,j-2)) + sqrt(g*Hx(i-1,j-2))*dt/dx;
       cv = nu.*(1-nu);
       wplusx = .5*cv.*(1-q);       
       
       duminus1 = H(i,j-1) - H(i,j-2);
       duminus2 = V(i,j-1) - V(i,j-2);
       duplus1 = H(i,j+1) - H(i,j);
       duplus2 = V(i,j+1) - V(i,j);
       duhalf1 = H(i,j) - H(i,j-1);
       duhalf2 = V(i,j) - V(i,j-1);
       rdenom = max(duhalf1.^2 + duhalf2.^2, eps);
       rnumplus = duplus1.*duhalf1 + duplus2.*duhalf2;
       rnumminus = duminus1.*duhalf1 + duminus2.*duhalf2;
       rplus = rnumplus./rdenom;
       rminus = rnumminus./rdenom;
       q = max(min(1.0, min(rminus, rplus)),0.0);
       nu = abs(Uy(i-2,j-2)) + sqrt(g*Hy(i-2,j-2))*dt/dy;
       cv = nu.*(1-nu);
       wminusy = .5*cv.*(1-q);
       
       duminus1 = H(i,j) - H(i,j-1);
       duminus2 = V(i,j) - V(i,j-1);
       duplus1 = H(i,j+2) - H(i,j+1);
       duplus2 = V(i,j+2) - V(i,j+1);
       duhalf1 = H(i,j+1) - H(i,j);
       duhalf2 = V(i,j+1) - V(i,j);
       rdenom = max(duhalf1.^2 + duhalf2.^2, eps);
       rnumplus = duplus1.*duhalf1 + duplus2.*duhalf2;
       rnumminus = duminus1.*duhalf1 + duminus2.*duhalf2;
       rplus = rnumplus./rdenom;
       rminus = rnumminus./rdenom;
       q = max(min(1.0, min(rminus, rplus)),0.0);
       nu = abs(Uy(i-2,j-1)) + sqrt(g*Hy(i-2,j-1))*dt/dy;
       cv = nu.*(1-nu);
       wplusy = .5*cv.*(1-q);
       
       if (max(max(abs(wplusx))) > .5) || (max(max(abs(wplusy))) > .5) || (max(max(abs(wminusx))) > .5) || (max(max(abs(wminusy))) > .5)
           max(max(abs(wplusx)))
           max(max(abs(wplusy)))
           max(max(abs(wminusx)))
           max(max(abs(wminusy)))
           max(max(abs(imag(wplusx))))
           max(max(abs(imag(wplusy))))
           max(max(abs(imag(wminusx))))
           max(max(abs(imag(wminusy))))
       end
       
       % height
       phi (i,j) = phi (i,j) - (dt/(2*dx))*(Ux (i,j-1)./Hx(i,j-1)+Ux (i-1,j-1)./Hx(i-1,j-1)) .* (phix(i,j-1)-phix(i-1,j-1)) - ...
                         (dt/(2*dy))*(Vy (i-1,j)./Hy(i-1,j)+Vy (i-1,j-1)./(Hy(i-1,j-1))) .* (phiy(i-1,j)-phiy(i-1,j-1));
                     
       % height
       H (i,j) = H (i,j) - (dt/dx)*(Ux (i,j-1)-Ux (i-1,j-1)) - wminusx.*(H(i,j) - H(i-1,j)) + wplusx.*(H(i+1,j)-H(i,j))...
                         - (dt/dy)*(Vy (i-1,j)-Vy (i-1,j-1)) - wminusy.*(H(i,j) - H(i,j-1)) + wplusy.*(H(i,j+1)-H(i,j));
       % x momentum
       U (i,j) = U (i,j) - (dt/dx)*((Ux (i,j-1).^2./Hx(i,j-1) + g/2*Hx (i,j-1).^2) - wminusx.*(U(i,j) - U(i-1,j)) + wplusx.*(U(i+1,j)-U(i,j)) - ...
                         (Ux (i-1,j-1).^2./Hx(i-1,j-1) + g/2*Hx (i-1,j-1).^2)) ...
                       - (dt/dy)*((Vy (i-1,j).*Uy (i-1,j)./Hy(i-1,j)) - wminusy.*(U(i,j) - U(i,j-1)) + wplusy.*(U(i,j+1)-U(i,j)) - ...
                         (Vy (i-1,j-1).*Uy (i-1,j-1)./Hy(i-1,j-1)));
       % y momentum
       V (i,j) = V (i,j) - (dt/dx)*((Ux (i,j-1).*Vx (i,j-1)./Hx(i,j-1)) + wminusx.*(V(i,j) - V(i-1,j)) - wplusx.*(V(i+1,j)-V(i,j)) - ...
                         (Ux (i-1,j-1).*Vx (i-1,j-1)./Hx(i-1,j-1))) ...
                       - (dt/dy)*((Vy (i-1,j).^2./Hy(i-1,j) + g/2*Hy (i-1,j).^2) + wminusy.*(V(i,j) - V(i,j-1)) - wplusy.*(V(i,j+1)-V(i,j)) - ...
                         (Vy (i-1,j-1).^2./Hy(i-1,j-1) + g/2*Hy (i-1,j-1).^2));
   
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
          set (phiplot,'zdata',H (i,j),'cdata',phi(i,j));
          drawnow
          
          figure(3);
          set (uplot,'zdata',U(i,j),'cdata',U(i,j));
          drawnow
          
          figure(4);
          set (vplot,'zdata',V (i,j),'cdata',V(i,j));
          drawnow
       end
      
       if all (all (isnan (H))), break, end  % Unstable, restart
   end
end
close (figure(1)), close(figure(2))

% ------------------------------------

function D = droplet (height,width)
% DROPLET  2D Gaussian
% D = droplet (height,width)
   [x,y] = ndgrid (-1:(2/(width-1)):1);
   D = height*exp (-5*(x.^2+y.^2));

% ------------------------------------

function [surfplot,phiplot,uplot,vplot,top,start,stop] = initgraphics (n);
% INITGRAPHICS  Initialize graphics for shallow_water.
% [surfplot,top,start,stop] = initgraphics (n)
% returns handles to a surface plot, its title, and two uicontrol toggles.

   if exist ('bigscreen','file')
      bigscreen
   end
   clf
   shg
   figure(1);
   set (gcf,'numbertitle','off','name','Shallow_water')
   x = (0:n-1)/(n-1);
   surfplot = surf (x,x,ones (n,n),zeros (n,n));
   grid off
   axis ([0 1 0 1 -1 3])
   caxis ([-1 1])
   shading faceted
   c = (1:64)'/64;
   cyan = [0*c c c];
   colormap (cyan)
   top = title ('Click start');
   start = uicontrol ('position',[20 20 80 20],'style','toggle','string','start');
   stop = uicontrol ('position',[120 20 80 20],'style','toggle','string','stop');
   
   figure(2);
   phiplot = surf (x,x,ones (n,n),zeros (n,n));
   grid off
   axis([0 1 0 1 -1 3])
   caxis  ([-1 1])
   shading faceted
   r = (1:64)'/64;
   red = [r r 0*r];
   colormap(red)
   
   figure(3);
   uplot = surf (x,x,ones (n,n),zeros (n,n));
   grid off
   axis([0 1 0 1 -1 3])
   caxis  ([-1 1])
   shading faceted
   r = (1:64)'/64;
   red = [r r 0*r];
   colormap(red)
   
   figure(4);
   vplot = surf (x,x,ones (n,n),zeros (n,n));
   grid off
   axis([0 1 0 1 -1 3])
   caxis  ([-1 1])
   shading faceted
   r = (1:64)'/64;
   red = [r r 0*r];
   colormap(red)