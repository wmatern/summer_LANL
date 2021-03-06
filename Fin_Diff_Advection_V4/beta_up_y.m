function b = beta_up_y(p,p_new,v,i,j,dt,dx)

b = (psi_max(p,i,j,p_new)-p(i,j))./...
    (dt/dx*abs(my_plus(v(i,j-1)).*p(i,j-1)-my_minus(v(i,j)).*p(i,j  ))+eps);
end

function p = psi_max(p,i,j,p_new)
p = max(max(max(p(i,j-1),p(i,j+1)),p(i,j)),...
    max(max(p_new(i,j-1),p_new(i,j+1)),p_new(i,j)));
end