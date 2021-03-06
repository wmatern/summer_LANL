function b = beta_dn_y(p,p_new,v,i,j,dt,dx)

b = (p(i,j)-psi_min(p,i,j,p_new))./...
    (dt/dx*abs(my_plus(v(i,j  )).*p(i,j  )-my_minus(v(i,j+1)).*p(i,j+1))+eps);
end

function p = psi_min(p,i,j,p_new)
p = min(min(min(p    (i,j-1),p    (i,j+1)),p    (i,j)),...
    min(min(    p_new(i,j-1),p_new(i,j+1)),p_new(i,j)));
end