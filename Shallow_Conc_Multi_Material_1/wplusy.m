function w = wplusy(H,i,j,Vy,Hy,g,C)
duminus1  = H(i,j)   - H(i,j-1);
duplus1   = H(i,j+2) - H(i,j+1);
duhalf1   = H(i,j+1) - H(i,j  );
rdenom    = max(duhalf1.^2, eps);
rnumplus  = duplus1 .*duhalf1;
rnumminus = duminus1.*duhalf1;
rplus     = rnumplus ./rdenom;
rminus    = rnumminus./rdenom;
q         = max(min(1.0, min(rminus, rplus)),0.0);
nu        = (abs(Vy(i-1,j)./Hy(i-1,j)) + sqrt(g*Hy(i-1,j)))*C;
cv        = nu.*(1-nu);
w	      = .5*cv.*(1-q);
end