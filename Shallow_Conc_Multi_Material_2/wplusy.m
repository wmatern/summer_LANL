function w = wplusy (H,i,j,nu,method)
duminus1  = H(i  ,j  ) - H(i  ,j-1);
duplus1   = H(i  ,j+2) - H(i  ,j+1);
duhalf1   = H(i  ,j+1) - H(i  ,j  );
rdenom    = max(duhalf1.^2, eps);
rnumplus  = duplus1 .*duhalf1;
rnumminus = duminus1.*duhalf1;
rplus     = rnumplus ./rdenom;
rminus    = rnumminus./rdenom;

if     method == 1 % use minmod
    q         = max(min(1.0, min(rminus, rplus)),0.0);
elseif method == 2 % use superbee
    q          = max(0.0, max(min(1.0, 2*min(rminus, rplus)),...
                    min(2.0, min(rminus, rplus))));
end

cv = nu.*(1-nu);
w  = 0.5*cv.*(1-q);
end