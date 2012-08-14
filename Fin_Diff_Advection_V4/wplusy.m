function w = wplusy (H,U,V,i,j,nu,method)
duminus1  = H(i  ,j  ) - H(i  ,j-1);
duplus1   = H(i  ,j+2) - H(i  ,j+1);
duhalf1   = H(i  ,j+1) - H(i  ,j  );
rdenom    = max(duhalf1.^2, eps);
rnumplus  = duplus1 .*duhalf1;
rnumminus = duminus1.*duhalf1;
rplus     = rnumplus ./rdenom;
rminus    = rnumminus./rdenom;

if V(i,j) > 0
    theta = rminus;
else
    theta = rplus;
end

if     method == 1 % use minmod
    q         = max(min(1.0, min(rminus, rplus)),0.0);
elseif method == 2 % use superbee
    q          = max(0.0, max(min(1.0, 2*theta),...
        min(2.0, theta)));
elseif method == 3 % upwind
    q          = 0;
else
    q          = 1;
end

cv = nu.*(1-nu);
w  = 0.5*cv.*(1-q);
end