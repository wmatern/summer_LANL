function w = wminusx(H,U,V,i,j,nu,method)
duminus   = H(i-1,j  ) - H(i-2,j  );
duplus    = H(i+1,j  ) - H(i  ,j  );
duhalf    = H(i  ,j  ) - H(i-1,j  );
rdenom    = max(duhalf.^2, eps);
rnumplus  = duplus .*duhalf;
rnumminus = duminus.*duhalf;
rplus     = rnumplus ./rdenom;
rminus    = rnumminus./rdenom;

if U(i,j) > 0
    theta = rminus;
else
    theta = rplus;
end

if     method == 1 % use minmod
    q         = max(min(1.0, min(rminus, rplus)),0.0);
elseif method == 2 % use superbee
    q          = max(0.0, max(min(1.0, 2*theta),...
        min(2.0, theta)));
else
    q          = 1;
end

cv = nu.*(1-nu);
w  = 0.5*cv.*(1-q);

end