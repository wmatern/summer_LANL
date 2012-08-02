function b = beta_up_y(p,i,j,A_in)

b = (max(p(i,j-1),max(p(i,j),p(i,j+1))) - p(i,j))./(A_in+eps);