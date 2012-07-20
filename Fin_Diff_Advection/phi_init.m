function init = phi_init(X,Y,nx,ny)
    xnorm = X/nx;
    init = ceil(xnorm-1/12) - ceil(xnorm-3/12);
    init = init + ((xnorm-5/12)/(1/12) - ceil(xnorm-6/12).*(2*(xnorm - 6/12)/(1/12))).*(ceil(xnorm-5/12) - ceil(xnorm-7/12));
    init = init + exp(-(150*(xnorm-10/12).^2));
end

