function [phi U V] = init_cond(X,Y,nx,ny,cond)
if cond == 1
    xnorm = X/nx;
    phi = ceil(xnorm-1/12) - ceil(xnorm-3/12);
    phi = phi + ((xnorm-5/12)/(1/12) - ceil(xnorm-6/12).*(2*(xnorm - 6/12)/(1/12))).*(ceil(xnorm-5/12) - ceil(xnorm-7/12));
    phi = phi + exp(-(150*(xnorm-10/12).^2));
    
    U  = ones (nx+4,ny+4);
    V  = zeros(nx+4,ny+4);
elseif cond == 2
    ynorm = Y/ny;
    phi = ceil(ynorm-1/12) - ceil(ynorm-3/12);
    phi = phi + ((ynorm-5/12)/(1/12) - ceil(ynorm-6/12).*(2*(ynorm - 6/12)/(1/12))).*(ceil(ynorm-5/12) - ceil(ynorm-7/12));
    phi = phi + exp(-(150*(ynorm-10/12).^2));
    
    V  = -ones (nx+4,ny+4);
    U  = zeros(nx+4,ny+4);
elseif cond == 3
    xnorm = X/nx;
    ynorm = Y/ny;
    phi = (ceil(xnorm-1/12) - ceil(xnorm-2/12)).*(ceil(ynorm-3/12) - ceil(ynorm-4/12));
    phi = phi + (ceil(xnorm-1/12) - ceil(xnorm-2/12)).*(ceil(ynorm-5/12) - ceil(ynorm-6/12));
    phi = phi + (ceil(xnorm-1/12) - ceil(xnorm-2/12)).*(ceil(ynorm-7/12) - ceil(ynorm-8/12));
    
    U = U_init(X,Y,nx,ny);
    V = V_init(X,Y,nx,ny);
elseif cond == 4
    phi = ones(nx,ny);
    phi(:,ny/2:ny) = 0.5;
    
    V  = ones (nx+4,ny+4);
    U  = zeros(nx+4,ny+4);
end
temp = zeros(nx+4,ny+4);
temp(3:nx+2,3:ny+2) = reshape(phi,nx,ny);
phi = temp;

end

