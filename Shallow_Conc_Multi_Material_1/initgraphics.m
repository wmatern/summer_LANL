function [surfplot,phiplot,heightplot,top,start,stop] = initgraphics (n,m)
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
   xlabel('X'); ylabel('Y');
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
   
   for k=3:m+2
       figure(k)
       grid off
       heightplot(k) = surf (x,x,ones (n,n),zeros (n,n));
       axis([0 1 0 1 0 3])
       caxis  ([-1 1])
       shading faceted
       r = (1:64)'/64;
       red = [r r 0*r];
       colormap(red)
   end
   