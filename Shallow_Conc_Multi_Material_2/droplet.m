function D = droplet (height,width)
% DROPLET  2D Gaussian
% D = droplet (height,width)
   [x,y] = ndgrid (-1:(2/(width-1)):1);
   D = height*exp (-5*(x.^2+y.^2));