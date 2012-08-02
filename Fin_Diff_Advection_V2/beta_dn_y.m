function b = beta_dn_y(p,i,j,A_out)

b = (p(i,j) - min(p(i,j-1),min(p(i,j),p(i,j+1))))./(A_out+eps);