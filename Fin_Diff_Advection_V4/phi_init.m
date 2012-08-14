function init = phi_init(X,Y,nx,ny,cond)
if cond == 0
    xnorm = X/nx;
    init = ceil(xnorm-1/12) - ceil(xnorm-3/12);
    init = init + ((xnorm-5/12)/(1/12) - ceil(xnorm-6/12).*(2*(xnorm - 6/12)/(1/12))).*(ceil(xnorm-5/12) - ceil(xnorm-7/12));
    init = init + exp(-(150*(xnorm-10/12).^2));
elseif cond == 2
    ynorm = Y/ny;
    init = ceil(ynorm-1/12) - ceil(ynorm-3/12);
    init = init + ((ynorm-5/12)/(1/12) - ceil(ynorm-6/12).*(2*(ynorm - 6/12)/(1/12))).*(ceil(ynorm-5/12) - ceil(ynorm-7/12));
    init = init + exp(-(150*(ynorm-10/12).^2));
elseif cond == 3
    xnorm = X/nx;
    ynorm = Y/ny;
    init = (ceil(xnorm-1/12) - ceil(xnorm-2/12)).*(ceil(ynorm-3/12) - ceil(ynorm-4/12));
    init = init + (ceil(xnorm-1/12) - ceil(xnorm-2/12)).*(ceil(ynorm-5/12) - ceil(ynorm-6/12));
    init = init + (ceil(xnorm-1/12) - ceil(xnorm-2/12)).*(ceil(ynorm-7/12) - ceil(ynorm-8/12));
elseif cond == 4
    ynorm = Y/ny;
    init = ceil(ynorm-0) - ceil(ynorm-0.5);
    
else
    init = zeros(nx,ny);
end
end

