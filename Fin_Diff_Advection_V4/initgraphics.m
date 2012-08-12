function [surfplot,phiplot,heightplot,top] = initgraphics (n,nx,ny)
% INITGRAPHICS  Initialize graphics for shallow_water.
% [surfplot,top,start,stop] = initgraphics (n)
% returns handles to a surface plot, its title, and two uicontrol toggles.

%    if exist ('bigscreen','file')
%       bigscreen
%    end
%    clf
%    shg
   figure(1);
%   set (gcf,'numbertitle','off','name','Shallow_water')
   x = (0:n-1)/(n-1);
   surfplot = surf (1:nx,1:ny,ones (n,n),zeros (n,n));
   grid off
   caxis ([-1 1])
   shading faceted
   c = (1:64)'/64;
   cyan = [0*c c c];
   colormap (cyan)
   xlabel('X'); ylabel('Y');
   top = title ('Click start');
%    start = uicontrol ('position',[20 20 80 20],'style','toggle','string','start');
%    stop = uicontrol ('position',[120 20 80 20],'style','toggle','string','stop');
   
   figure(2);
   phiplot = surf (1:nx,1:ny,ones (n,n),zeros (n,n));
   grid off
   caxis  ([-1 1])
   shading faceted
   r = (1:64)'/64;
   red = [r r 0*r];
   colormap(red)
   
   
   figure(3)
   grid off
   heightplot(k) = surf (1:nx,1:ny,ones (n,n),zeros (n,n));
   axis ([X(1,1) X(n,n) Y(1,1) Y(n,n) 0 3])
   caxis  ([-1 1])
   shading faceted
   r = (1:64)'/64;
   red = [r r 0*r];
   colormap(red)
   